#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "TissueSimulation.hpp"
#include "Exception.hpp"
#include "CancerParameters.hpp"
#include "RandomNumberGenerator.hpp"
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <set>
#include "TrianglesMeshWriter.cpp"
#include "TrianglesMeshReader.hpp"
#include "SimulationTime.hpp"
#include "ColumnDataWriter.hpp"
#include "WntCellCycleModel.hpp"
#include "WntGradient.hpp"
#include "OutputFileHandler.hpp"
#include "LogFile.cpp"

template<unsigned DIM> 
TissueSimulation<DIM>::TissueSimulation(Crypt<DIM>& rCrypt, bool deleteCrypt)
  :  mrCrypt(rCrypt)
{
    srandom(0);
    mDeleteCrypt = deleteCrypt;
    mpParams = CancerParameters::Instance();
    mpRandomGenerator = RandomNumberGenerator::Instance();
    
    mDt = 1.0/120.0; // Timestep of 30 seconds (as per Meineke)
    mEndTime = 0.0; // hours - this is set later on.
    
    // defaults
    mOutputDirectory = "";
    mReMesh = true;
    mOutputCellTypes = false ;
    mNoBirth = false;
    mMaxCells = 10*mrCrypt.rGetMesh().GetNumNodes();
    mMaxElements = 10*mrCrypt.rGetMesh().GetNumElements();
    mWntIncluded = false;
    mNumBirths = 0;
    mNumDeaths = 0;
    mUseNonFlatBottomSurface = false;
    
    // whether to use a cutoff point, ie specify zero force if two 
    // cells are greater than a certain distance apart. By default 
    // we don't have a cutoff point
    mUseCutoffPoint = false;
    mCutoffPoint = 1e10;
    
    
    assert(SimulationTime::Instance()->IsStartTimeSetUp());
    // start time must have been set to create crypt which includes cell cycle models
    
    mrCrypt.SetMaxCells(mMaxCells);
    mrCrypt.SetMaxElements(mMaxElements);
}

/**
 * Free any memory allocated by the constructor.
 * This frees the crypt and cell killers, if they were created by de-serialization.
 */
template<unsigned DIM> 
TissueSimulation<DIM>::~TissueSimulation()
{
    if (mDeleteCrypt)
    {
        delete &mrCrypt;
        for (typename std::vector<AbstractCellKiller<DIM>*>::iterator it=mCellKillers.begin();
             it != mCellKillers.end();
             ++it)
        {
            delete *it;
        }
    }
}

template<unsigned DIM> 
void TissueSimulation<DIM>::WriteVisualizerSetupFile(std::ofstream& rSetupFile)
{
    assert(DIM==2); // this is 2d specific
    rSetupFile << "MeshWidth\t" << mrCrypt.rGetMesh().GetWidth(0u);// get furthest distance between nodes in the x-direction
    rSetupFile.close();
}

template<unsigned DIM>  
unsigned TissueSimulation<DIM>::DoCellBirth()
{
    //assert (!mNoBirth);
    if (mNoBirth)
    {
        return 0;
    }

    unsigned num_births_this_step = 0;

    // Iterate over all cells, seeing if each one can be divided
    for (typename Crypt<DIM>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        MeinekeCryptCell& cell = *cell_iter;
        Node<DIM>* p_our_node = cell_iter.GetNode();
        
        // Check for this cell dividing
        // Construct any influences for the cell cycle...

        std::vector<double> cell_cycle_influences;
        if (mWntIncluded)
        {
            #define COVERAGE_IGNORE
            assert(DIM==2);
            #undef COVERAGE_IGNORE
            double y = p_our_node->rGetLocation()[1];
            double wnt_stimulus = mWntGradient.GetWntLevel(y);
            cell_cycle_influences.push_back(wnt_stimulus);
        }
        
        // check if this cell is ready to divide - if so create a new cell etc.
        if (cell.ReadyToDivide(cell_cycle_influences))
        {
            // Create new cell
            MeinekeCryptCell new_cell = cell.Divide();
            // std::cout << "Cell division at node " << cell.GetNodeIndex() << "\n";
        
            // Add a new node to the mesh
            c_vector<double, DIM> new_location = CalculateDividingCellCentreLocations(cell_iter);
            
            MeinekeCryptCell *p_new_cell=mrCrypt.AddCell(new_cell, new_location);
            std::set<MeinekeCryptCell*> new_cell_pair;
            new_cell_pair.insert(&cell); //Parent cell
            new_cell_pair.insert(p_new_cell); //New cell (the clue's in the name)
            
            mDivisionPairs.insert(new_cell_pair);
            num_births_this_step++;
        } // if (ready to divide)
    } // cell iteration loop
   
    return num_births_this_step;
}


template<unsigned DIM> 
unsigned TissueSimulation<DIM>::DoCellRemoval()
{
    unsigned num_deaths_this_step=0;
        
    // this labels cells as dead or apoptosing. It does not actually remove the cells, 
    // crypt.RemoveDeadCells() needs to be called for this.
    for(unsigned killer_index = 0; killer_index<mCellKillers.size(); killer_index++)
    {
        mCellKillers[killer_index]->TestAndLabelCellsForApoptosisOrDeath();
    }
    
    num_deaths_this_step += mrCrypt.RemoveDeadCells(); 
    
    return num_deaths_this_step;
}



template<unsigned DIM> 
c_vector<double, DIM> TissueSimulation<DIM>::CalculateDividingCellCentreLocations(typename Crypt<DIM>::Iterator parentCell)
{
    double separation = CancerParameters::Instance()->GetDivisionSeparation();
    c_vector<double, DIM> parent_coords = parentCell.rGetLocation();
    c_vector<double, DIM> daughter_coords;
        
    // pick a random direction and move the parent cell backwards by 0.5*sep in that
    // direction and return the position of the daughter cell (0.5*sep forwards in the
    // random vector direction

    // Make a random direction vector of the required length
    c_vector<double, DIM> random_vector;
    
    if(DIM==1)
    {
        random_vector(0) = 0.5*separation;

        daughter_coords = parent_coords+random_vector;
        parent_coords = parent_coords-random_vector;
    }   
    else if(DIM==2)
    {
        double random_angle = RandomNumberGenerator::Instance()->ranf();
        random_angle *= 2.0*M_PI;
        
        random_vector(0) = 0.5*separation*cos(random_angle);
        random_vector(1) = 0.5*separation*sin(random_angle);
        
        c_vector<double, 2> proposed_new_parent_coords = parent_coords-random_vector;
        c_vector<double, 2> proposed_new_daughter_coords = parent_coords+random_vector;
        
        if (   (proposed_new_parent_coords(1) >= BottomSurfaceProfile(proposed_new_parent_coords(0)))
            && (proposed_new_daughter_coords(1) >= BottomSurfaceProfile(proposed_new_daughter_coords(0))))
        {
             // We are not too close to the bottom of the crypt
            // move parent
            parent_coords = proposed_new_parent_coords;
            daughter_coords = proposed_new_daughter_coords;
        }
        else
        {
            proposed_new_daughter_coords = parent_coords+2.0*random_vector;
            while (proposed_new_daughter_coords(1) < BottomSurfaceProfile(proposed_new_daughter_coords(0)))
            {
                random_angle = RandomNumberGenerator::Instance()->ranf();
                random_angle *= 2.0*M_PI;
                
                random_vector(0) = separation*cos(random_angle);
                random_vector(1) = separation*sin(random_angle);
                proposed_new_daughter_coords = parent_coords+random_vector;
            }
            daughter_coords = proposed_new_daughter_coords;
        }
        
//        if (   (proposed_new_parent_coords(1) <= BottomSurfaceProfile(proposed_new_parent_coords(0)))
//            && (proposed_new_daughter_coords(1) <= BottomSurfaceProfile(proposed_new_daughter_coords(0))))
//        {
//            // the worst case, both parent and daughter would move to below the surface
//            // move them across and then up onto the surface. This is the only case
//            // when the cells won't be 0.1 away from each other
//            parent_coords = proposed_new_parent_coords;
//            parent_coords(1) = BottomSurfaceProfile(parent_coords(0));
//
//            daughter_coords = proposed_new_daughter_coords;
//            daughter_coords(1) = BottomSurfaceProfile(daughter_coords(0));
//        }
//        else if (proposed_new_parent_coords(1) <= BottomSurfaceProfile(proposed_new_parent_coords(0)))
//        {
//            // Leave parent, move daughter twice as far
//            daughter_coords = parent_coords + random_vector;
//        }
//        else if (proposed_new_daughter_coords(1) <= BottomSurfaceProfile(proposed_new_daughter_coords(0)))
//        {
//            // don't move parent, move daughter in opposite direction (twice as far)
//            daughter_coords = parent_coords - random_vector;
//        }
//        else        
//        {   
//            // We are not too close to the bottom of the crypt
//            // move parent
//            parent_coords = proposed_new_parent_coords;
//            daughter_coords = proposed_new_daughter_coords;
//        }
        
        
        

        assert(daughter_coords(1)>=0.0);// to make sure dividing cells stay in the crypt
        assert(parent_coords(1)>=0.0);// to make sure dividing cells stay in the crypt
    }
    else if(DIM==3)
    {
        double random_zenith_angle = RandomNumberGenerator::Instance()->ranf();// phi 
        random_zenith_angle *= M_PI;
        double random_azimuth_angle = RandomNumberGenerator::Instance()->ranf();// theta
        random_azimuth_angle *= 2*M_PI;
        
        random_vector(0) = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
        random_vector(1) = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
        random_vector(2) = 0.5*separation*cos(random_zenith_angle);

        daughter_coords = parent_coords+random_vector;
        parent_coords = parent_coords-random_vector;
    }
    
    
    // set the parent to use this location
    ChastePoint<DIM> parent_coords_point(parent_coords);
    mrCrypt.MoveCell(parentCell, parent_coords_point);
    return daughter_coords;
}



template<unsigned DIM>  
std::vector<c_vector<double, DIM> > TissueSimulation<DIM>::CalculateVelocitiesOfEachNode()
{
    std::vector<c_vector<double, DIM> > drdt(mrCrypt.rGetMesh().GetNumAllNodes());
    for (unsigned i=0; i<drdt.size(); i++)
    {
        drdt[i]=zero_vector<double>(DIM);
    }

    for(typename Crypt<DIM>::SpringIterator spring_iterator=mrCrypt.SpringsBegin();
        spring_iterator!=mrCrypt.SpringsEnd();
        ++spring_iterator)
    {
        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
         
        double damping_constantA = mpParams->GetDampingConstantNormal();
        double damping_constantB = mpParams->GetDampingConstantNormal();
        
        
        if(   (spring_iterator.rGetCellA().GetMutationState()==HEALTHY)
           || (spring_iterator.rGetCellA().GetMutationState()==APC_ONE_HIT))
        {
            damping_constantA = mpParams->GetDampingConstantNormal();
        }
        else
        {
            damping_constantA = mpParams->GetDampingConstantMutant();
        }
        
        if(   (spring_iterator.rGetCellB().GetMutationState()==HEALTHY)
           || (spring_iterator.rGetCellB().GetMutationState()==APC_ONE_HIT))
        {
            damping_constantB = mpParams->GetDampingConstantNormal();
        }
        else
        {
            damping_constantB = mpParams->GetDampingConstantMutant();
        }    
       
        // these cannot be ghost nodes anymore..
        // the both apply forces on each other
        drdt[nodeB_global_index] -= force / damping_constantB;
        drdt[nodeA_global_index] += force / damping_constantA;
    }
    
    return drdt;
}



template<unsigned DIM> 
c_vector<double, DIM> TissueSimulation<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex)
{
    assert(nodeAGlobalIndex!=nodeBGlobalIndex);
    c_vector<double, DIM> unit_difference;
    c_vector<double, DIM> node_a_location = mrCrypt.rGetMesh().GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = mrCrypt.rGetMesh().GetNode(nodeBGlobalIndex)->rGetLocation();
    
    // there is reason not to substract one position from the other (cyclidrical meshes). clever gary
    unit_difference = mrCrypt.rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);   
    
    double distance_between_nodes = norm_2(unit_difference);
    
    unit_difference /= distance_between_nodes;
    
    if(mUseCutoffPoint)
    {
        if( distance_between_nodes >= mCutoffPoint )
        {
            return zero_vector<double>(DIM);
            //c_vector<double,DIM>() is not guaranteed to be fresh memory
            // ie return zero force;
        }
    }
    
    double rest_length = 1.0;
        
    double ageA = mrCrypt.rGetCellAtNodeIndex(nodeAGlobalIndex).GetAge();
    double ageB = mrCrypt.rGetCellAtNodeIndex(nodeBGlobalIndex).GetAge();
    
    if (ageA<1.0 && ageB<1.0 )
    {
        // Spring Rest Length Increases to normal rest length from ???? to normal rest length, 1.0, over 1 hour
        std::set<MeinekeCryptCell *> cell_pair;
        cell_pair.insert(&(mrCrypt.rGetCellAtNodeIndex(nodeAGlobalIndex)));
        cell_pair.insert(&(mrCrypt.rGetCellAtNodeIndex(nodeBGlobalIndex)));
        unsigned count=mDivisionPairs.count(cell_pair);
        if (count==1)
        {   
            double lambda=CancerParameters::Instance()->GetDivisionRestingSpringLength();
            rest_length=(lambda+(1.0-lambda)*ageA);           
        }
        
        if (ageA+ SimulationTime::Instance()->GetTimeStep() >= 1.0)
        {
            //This spring is about to go out of scope
            mDivisionPairs.erase(cell_pair);
            //\todo some of the pairs aren't accessible at this stage(if we lose them due to a remesh
        }
       
    }
    assert(rest_length<=1.0);
    return mpParams->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
}



template<unsigned DIM> 
void TissueSimulation<DIM>::UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rDrDt)
{
    // update ghost positions first because they do not affect the real cells
    mrCrypt.UpdateGhostPositions(mDt);

    // Iterate over all cells to update their positions.
    for (typename Crypt<DIM>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        MeinekeCryptCell& cell = *cell_iter;
        unsigned index = cell.GetNodeIndex();
        
        ChastePoint<DIM> new_point(mrCrypt.rGetMesh().GetNode(index)->rGetLocation() + mDt*rDrDt[index]);
        
        if(DIM==2)
        {
               
            // stem cells are fixed if no wnt, so reset the x-value to the old x-value           
            if((cell.GetCellType()==STEM) && (!mWntIncluded))
            {
                new_point.rGetLocation()[0] = mrCrypt.rGetMesh().GetNode(index)->rGetLocation()[0];
                new_point.rGetLocation()[1] = mrCrypt.rGetMesh().GetNode(index)->rGetLocation()[1];
            }
            
            // for all cells - move up if below the bottom surface
            if (new_point.rGetLocation()[1] < BottomSurfaceProfile(new_point.rGetLocation()[0]) )
            {
                new_point.rGetLocation()[1] = BottomSurfaceProfile(new_point.rGetLocation()[0]);
            }
            
//          // move the cell
            mrCrypt.MoveCell(cell_iter, new_point);                    
        }
        else
        {   // 1D or 3D
            mrCrypt.MoveCell(cell_iter, new_point);    
        }
    }
}


/**
 * Change the state of cells
 *
 * At the moment this turns cells to be differentiated
 * dependent on a protein concentration when using the Wnt model.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::UpdateCellTypes()
{
    // Designate cells as proliferating (transit) or
    // quiescent (differentiated) according to protein concentrations
    for (typename Crypt<DIM>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        cell_iter->UpdateCellType();
    }
    
}


/**
 * Set the timestep of the simulation
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetDt(double dt)
{
    assert(dt>0);
    mDt=dt;
}

/**
 * Get the timestep of the simulation
 */
template<unsigned DIM> 
double TissueSimulation<DIM>::GetDt()
{
    return mDt;
}

/**
 * Sets the end time and resets the timestep to be endtime/100
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetEndTime(double endTime)
{
    assert(endTime>0);
    mEndTime=endTime;
}

/**
 * Set the output directory of the simulation.
 * 
 * Note that tabulated results (for test comparison) go into a /tab_results subfolder
 * And visualizer results go into a /vis_results subfolder.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
}

/**
 * Sets the maximum number of cells that the simulation will contain (for use by the datawriter)
 * default value is set to 10x the initial mesh value by the constructor.
 */
template<unsigned DIM>  
void TissueSimulation<DIM>::SetMaxCells(unsigned maxCells)
{
    mMaxCells = maxCells;
    if (maxCells<mrCrypt.rGetMesh().GetNumAllNodes())
    {
        EXCEPTION("mMaxCells is less than the number of cells in the mesh.");
    }
    mrCrypt.SetMaxCells(maxCells);
}

/**
 * Sets the maximum number of elements that the simulation will contain (for use by the datawriter)
 * default value is set to 10x the initial mesh value by the constructor.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetMaxElements(unsigned maxElements)
{
    mMaxElements = maxElements;
    if (maxElements<mrCrypt.rGetMesh().GetNumAllElements())
    {
        EXCEPTION("mMaxElements is less than the number of elements in the mesh.");
    }
    mrCrypt.SetMaxElements(maxElements);
}


template<unsigned DIM> 
Crypt<DIM>& TissueSimulation<DIM>::rGetCrypt()
{
    return mrCrypt;
}

template<unsigned DIM> 
const Crypt<DIM>& TissueSimulation<DIM>::rGetCrypt() const
{
    return mrCrypt;
}


/**
 * Set whether the mesh should be remeshed at every time step.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetReMeshRule(bool remesh)
{
    mReMesh = remesh;
}


/**
 * Set the simulation to run with no birth.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetNoBirth(bool nobirth)
{
    mNoBirth = nobirth;
}

/**
 * Set the simulation to Count and store the number of each cell type.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetOutputCellTypes(bool output_cell_types)
{
    mOutputCellTypes = output_cell_types;
}


/**
 * Use a cutoff point, ie specify zero force if two cells are greater 
 * than the cutoff distance apart
 */
template<unsigned DIM>
void TissueSimulation<DIM>::UseCutoffPoint(double cutoffPoint)
{
    assert(cutoffPoint > 0.0);
    mUseCutoffPoint = true;
    mCutoffPoint = cutoffPoint;
}


/**
 * This automatically sets this to be a wnt dependent simulation.
 * You should supply cells with a wnt cell cycle...
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetWntGradient(WntGradientType wntGradientType)
{
    mWntIncluded = true;
    mWntGradient = WntGradient(wntGradientType);
}


/**
 * Add a cell killer to be used in this simulation
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::AddCellKiller(AbstractCellKiller<DIM>* pCellKiller)
{
    mCellKillers.push_back(pCellKiller);
}

/**
 * Get a node's location (ONLY FOR TESTING)
 *
 * @param the node index
 * @return the co-ordinates of this node.
 */
template<unsigned DIM> 
std::vector<double> TissueSimulation<DIM>::GetNodeLocation(const unsigned& rNodeIndex)
{
    std::vector<double> location;
    for(unsigned i=0; i<DIM; i++)
    {
        location.push_back( mrCrypt.rGetMesh().GetNode(rNodeIndex)->rGetLocation()[i] );
    }
    return location;
}

/**
 * Main Solve method.
 *
 * Once CryptSimulation object has been set up, call this to run simulation
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::Solve()
{ 
    // Set up the simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetDimensionalisedTime();
    
    unsigned num_time_steps = (unsigned) ((mEndTime-current_time)/mDt+0.5);

    if (current_time>0)//use the reset function if necessary
    {
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
    }
    else
    {
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
    }
    

    if (mOutputDirectory=="")
    {
        EXCEPTION("OutputDirectory not set");
    }
    

    double time_now = p_simulation_time->GetDimensionalisedTime();
    std::ostringstream time_string;
    time_string << time_now;
    
    std::string results_directory = mOutputDirectory +"/results_from_time_" + time_string.str();
    
    
    ///////////////////////////////////////////////////////////
    // Set up Simulation
    ///////////////////////////////////////////////////////////
    
    // Data writers for tabulated results data, used in tests
    // first construction clears out the folder
    ColumnDataWriter tabulated_node_writer(results_directory+"/tab_results", "tabulated_node_results",true);
    ColumnDataWriter tabulated_element_writer(results_directory+"/tab_results", "tabulated_element_results",false);
    
    mrCrypt.SetupTabulatedWriters(tabulated_node_writer, tabulated_element_writer);//, element_writer_ids);
    
    // This keeps track of when tabulated results were last output
    unsigned tabulated_output_counter = 0;
    
    // Create output files for the visualizer
    OutputFileHandler output_file_handler(results_directory+"/vis_results/",false);
    out_stream p_node_file = output_file_handler.OpenOutputFile("results.viznodes");
    out_stream p_element_file = output_file_handler.OpenOutputFile("results.vizelements");
    out_stream p_setup_file = output_file_handler.OpenOutputFile("results.vizsetup");
    // Creates output file to store number of different cells
    
    out_stream p_cell_types_file = output_file_handler.OpenOutputFile("celltypes.dat");
    if (mOutputCellTypes)
    {
        *p_cell_types_file <<   "Time\t Healthy\t Labelled\t APC_1\t APC_2\t BETA_CAT \n";
    }
    
    
    
    /* Age the cells to the correct time (cells set up with negative birth dates
     * to give some that are almost ready to divide).
     * 
     * TODO:For some strange reason this seems to take about 3 minutes for a realistic Wnt-Crypt.
     * Not sure why - when the same code was evaluated in a test it seemed almost instant.
     */
    
    for (typename Crypt<DIM>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        double y = cell_iter.rGetLocation()[1];
        std::vector<double> cell_cycle_influences;
        if (mWntIncluded)
        {
            double wnt_stimulus = mWntGradient.GetWntLevel(y);
            cell_cycle_influences.push_back(wnt_stimulus);
        }
        // We don't use the result; this call is just to force the cells to age to time 0,
        // running their cell cycle models to get there.
        cell_iter->ReadyToDivide(cell_cycle_influences);
    }
    
    UpdateCellTypes();
    
    // Write initial conditions to file for the visualizer.
    if(DIM==2)
    {
        WriteVisualizerSetupFile(*p_setup_file);
    }
    
    mrCrypt.WriteResultsToFiles(tabulated_node_writer, 
                               tabulated_element_writer,
                               *p_node_file, *p_element_file, *p_cell_types_file,
                               false,
                               true,
                               mOutputCellTypes);
                               
                               
    /////////////////////////////////////////////////////////////////////
    // Main time loop
    /////////////////////////////////////////////////////////////////////
    while (p_simulation_time->GetTimeStepsElapsed() < num_time_steps)
    {        
        LOG(1, "--TIME = " << p_simulation_time->GetDimensionalisedTime() << "\n");
        
        // remove dead cells before doing birth
        // neither of these functions use any element information so they 
        // just delete and create nodes
        mNumDeaths += DoCellRemoval();
        LOG(1, "\tNum deaths = " << mNumDeaths << "\n");

        mNumBirths += DoCellBirth();
        LOG(1, "\tNum births = " << mNumBirths << "\n");
        

        if( (mNumBirths>0) || (mNumDeaths>0) )
        {   
            // If any nodes have been deleted or added we MUST call a ReMesh
            assert(mReMesh);
        }

        if(mReMesh)
        {
            LOG(1, "\tRemeshing..\n");
            mrCrypt.ReMesh();
        }

        //  calculate node velocities
        std::vector<c_vector<double, DIM> > drdt = CalculateVelocitiesOfEachNode();

        // update node positions
        UpdateNodePositions(drdt);
        
        // Change the state of some cells
        // Only active for WntCellCycleModel at the moment
        // but mutations etc. could occur in this function
        UpdateCellTypes();
                
        
        // Increment simulation time here, so results files look sensible
        p_simulation_time->IncrementTimeOneStep();
        
        // Write results to file
        mrCrypt.WriteResultsToFiles(tabulated_node_writer, 
                                    tabulated_element_writer, 
                                    *p_node_file, *p_element_file, *p_cell_types_file,
                                    tabulated_output_counter%80==0,
                                    true,
                                    mOutputCellTypes);
        
        
        
        
                            
        tabulated_output_counter++;
    }
    
    // Write end state to tabulated files (not visualizer - this
    // is taken care of in the main loop).
    // Doesn't need to count cell types again as it is done in the last loop
    mrCrypt.WriteResultsToFiles(tabulated_node_writer, 
                                tabulated_element_writer, 
                                *p_node_file, *p_element_file, *p_cell_types_file,
                                true,
                                false,
                                false);
                        
    tabulated_node_writer.Close();
    tabulated_element_writer.Close();
}



/**
 * Saves the whole crypt simulation for restarting later.
 *
 * Puts it in the folder mOutputDirectory/archive/
 * and the file "2dCrypt_at_time_<SIMULATION TIME>.arch"
 *
 * First archives simulation time then the simulation itself.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::Save()
{
    SimulationTime* p_sim_time = SimulationTime::Instance();
    assert(p_sim_time->IsStartTimeSetUp());
    
    std::string archive_directory = mOutputDirectory + "/archive/";
    
    std::ostringstream time_stamp;
    time_stamp << p_sim_time->GetDimensionalisedTime();
    
    // create an output file handler in order to get the full path of the
    // archive directory. Note the false is so the handler doesn't clean
    // the directory
    OutputFileHandler handler(archive_directory, false);
    std::string archive_filename = handler.GetTestOutputDirectory() + "2dCrypt_at_time_"+time_stamp.str()+".arch";
    std::string mesh_filename = std::string("mesh_") + time_stamp.str();
    
    if(mReMesh)
    {
        mrCrypt.ReMesh();
    }

    
    // the false is so the directory isn't cleaned
    TrianglesMeshWriter<DIM,DIM> mesh_writer(archive_directory, mesh_filename, false);
    mesh_writer.WriteFilesUsingMesh(mrCrypt.rGetMesh());
    
    std::ofstream ofs(archive_filename.c_str());
    boost::archive::text_oarchive output_arch(ofs);
    
    // cast to const.
    const SimulationTime* p_simulation_time = SimulationTime::Instance();
    output_arch << *p_simulation_time;
    TissueSimulation<DIM> * p_sim = this;
    output_arch & p_sim;
}

/**
 * Loads a saved crypt simulation to run further.
 *
 * @param rArchiveDirectory the name of the simulation to load
 * (specified originally by simulator.SetOutputDirectory("wherever"); )
 * @param rTimeStamp the time at which to load the simulation (this must
 * be one of the times at which the simulation.Save() was called)
 */
template<unsigned DIM> 
TissueSimulation<DIM>* TissueSimulation<DIM>::Load(const std::string& rArchiveDirectory, const double& rTimeStamp)
{
    // Find the right archive and mesh to load
    std::ostringstream time_stamp;
    time_stamp << rTimeStamp;
    
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    
    OutputFileHandler any_old_handler("",false);
    std::string test_output_directory = any_old_handler.GetTestOutputDirectory();
    
    std::string archive_filename = test_output_directory + rArchiveDirectory + "/archive/2dCrypt_at_time_"+time_stamp.str() +".arch";
    std::string mesh_filename = test_output_directory + rArchiveDirectory + "/archive/mesh_" + time_stamp.str();
    Crypt<DIM>::meshPathname = mesh_filename;
    
    // Create an input archive
    std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
    boost::archive::text_iarchive input_arch(ifs);

    // Read the archive
    assert(p_simulation_time->IsStartTimeSetUp());
    input_arch >> *p_simulation_time;
    TissueSimulation<DIM>* p_sim;
    input_arch >> p_sim;
    
    if (p_sim->rGetCrypt().rGetMesh().GetNumNodes()!=p_sim->rGetCrypt().rGetCells().size())
    {
        #define COVERAGE_IGNORE
        std::stringstream string_stream;
        string_stream << "Error in Load(), number of nodes (" << p_sim->rGetCrypt().rGetMesh().GetNumNodes()
                      << ") is not equal to the number of cells (" << p_sim->rGetCrypt().rGetCells().size() 
                      << ")";
        EXCEPTION(string_stream.str());
        #undef COVERAGE_IGNORE
    }
    
    return p_sim;
}


template<unsigned DIM> 
double TissueSimulation<DIM>::BottomSurfaceProfile(double x)
{
    if(mUseNonFlatBottomSurface)
    {
        double crypt_width = mrCrypt.rGetMesh().GetWidth(0u);
        double frequency = floor(crypt_width/4.0)+1.0;
        
        return 0.05*(sin(2*frequency*M_PI*x/crypt_width)+1);
    }
    else
    {
        return 0;
    }
}

template<unsigned DIM>
void TissueSimulation<DIM>::UseNonFlatBottomSurface()
{
    assert(DIM==2);
    mUseNonFlatBottomSurface = true;
}

