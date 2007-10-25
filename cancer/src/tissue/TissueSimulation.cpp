#ifndef _TISSUESIMULATION_CPP_
#define _TISSUESIMULATION_CPP_

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
#include "LogFile.hpp"
#include "VoronoiTessellation.cpp"

template<unsigned DIM> 
TissueSimulation<DIM>::TissueSimulation(Tissue<DIM>& rTissue, bool deleteTissue)
  :  mrTissue(rTissue)
{
    srandom(0);
    mDeleteTissue = deleteTissue;
    mpParams = CancerParameters::Instance();
    mpRandomGenerator = RandomNumberGenerator::Instance();
    
    mDt = 1.0/120.0; // Timestep of 30 seconds (as per Meineke)
    mEndTime = 0.0; // hours - this is set later on.
    
    // defaults
    mOutputDirectory = "";
    mReMesh = true;
    mOutputCellTypes = false ;
    mNoBirth = false;
    mMaxCells = 10*mrTissue.rGetMesh().GetNumNodes();
    mMaxElements = 10*mrTissue.rGetMesh().GetNumElements();
    mNumBirths = 0;
    mNumDeaths = 0;
    mUseEdgeBasedSpringConstant = false;
    mUseAreaBasedViscosity = false;
	mUseMutantSprings = false;
    mWriteVoronoiData = false;
    mCreateVoronoiTessellation = false;
    mFollowLoggedCell = false;
    
    // whether to use a cutoff point, ie specify zero force if two 
    // cells are greater than a certain distance apart. By default 
    // we don't have a cutoff point
    mUseCutoffPoint = false;
    mCutoffPoint = 1e10;
    
    mrTissue.SetMaxCells(mMaxCells);
    mrTissue.SetMaxElements(mMaxElements);
}

/**
 * Free any memory allocated by the constructor.
 * This frees the tissue and cell killers, if they were created by de-serialization.
 */
template<unsigned DIM> 
TissueSimulation<DIM>::~TissueSimulation()
{
    if (mDeleteTissue)
    {
        delete &mrTissue;
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
    rSetupFile << "MeshWidth\t" << mrTissue.rGetMesh().GetWidth(0u);// get furthest distance between nodes in the x-direction
    rSetupFile.close();
}

template<unsigned DIM>  
unsigned TissueSimulation<DIM>::DoCellBirth()
{
    if (mNoBirth)
    {
        return 0;
    }

    unsigned num_births_this_step = 0;

    // Iterate over all cells, seeing if each one can be divided
    for (typename Tissue<DIM>::Iterator cell_iter = mrTissue.Begin();
         cell_iter != mrTissue.End();
         ++cell_iter)
    {
        TissueCell& cell = *cell_iter;

        // check if this cell is ready to divide - if so create a new cell etc.
        if (cell.ReadyToDivide())
        {
            // Create new cell
            TissueCell new_cell = cell.Divide();
        
            // Add a new node to the mesh
            c_vector<double, DIM> new_location = CalculateDividingCellCentreLocations(cell_iter);
            
            TissueCell *p_new_cell=mrTissue.AddCell(new_cell, new_location);
            mrTissue.MarkSpring(cell, *p_new_cell);
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
    // tissue.RemoveDeadCells() needs to be called for this.
    for(unsigned killer_index = 0; killer_index<mCellKillers.size(); killer_index++)
    {
        mCellKillers[killer_index]->TestAndLabelCellsForApoptosisOrDeath();
    }
    
    num_deaths_this_step += mrTissue.RemoveDeadCells(); 
    
    return num_deaths_this_step;
}



template<unsigned DIM> 
c_vector<double, DIM> TissueSimulation<DIM>::CalculateDividingCellCentreLocations(typename Tissue<DIM>::Iterator parentCell)
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
        
        parent_coords = parent_coords-random_vector;
        daughter_coords = parent_coords+random_vector;        
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
    mrTissue.MoveCell(parentCell, parent_coords_point);
    return daughter_coords;
}



template<unsigned DIM>  
std::vector<c_vector<double, DIM> > TissueSimulation<DIM>::CalculateVelocitiesOfEachNode()
{
    std::vector<c_vector<double, DIM> > drdt(mrTissue.rGetMesh().GetNumAllNodes());
    for (unsigned i=0; i<drdt.size(); i++)
    {
        drdt[i]=zero_vector<double>(DIM);
    }

    for(typename Tissue<DIM>::SpringIterator spring_iterator=mrTissue.SpringsBegin();
        spring_iterator!=mrTissue.SpringsEnd();
        ++spring_iterator)
    {
        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
         
        double damping_constantA;
        double damping_constantB;
        
        if (mUseAreaBasedViscosity)
        {
            assert(DIM==2);
            double rest_length = 1.0;
            double d0 = 0.1;
            // this number is such that d0+A*d1=1, where A is the area of a equilibrium
            // cell (=sqrt(3)/4 = a third of the area of a hexagon with edges of size 1)
            double d1 = 1.8/(sqrt(3)*rest_length*rest_length); 

            VoronoiTessellation<DIM>& tess = mrTissue.rGetVoronoiTessellation();
        
            double area_cell_A = tess.GetFaceArea(nodeA_global_index);
            double area_cell_B = tess.GetFaceArea(nodeB_global_index);
            
            // the areas should be order 1, this is just to avoid getting infinite areas
            // if an area based viscosity option is chosen without ghost nodes.
            assert(area_cell_A < 1000);
            assert(area_cell_B < 1000);
            
            damping_constantA = (d0 + area_cell_A*d1)*mpParams->GetDampingConstantNormal();
            damping_constantB = (d0 + area_cell_B*d1)*mpParams->GetDampingConstantNormal();
        }
        else
        {
            damping_constantA = mpParams->GetDampingConstantNormal();
            damping_constantB = mpParams->GetDampingConstantNormal();
        }

        if(   (spring_iterator.rGetCellA().GetMutationState()!=HEALTHY)
           && (spring_iterator.rGetCellA().GetMutationState()!=APC_ONE_HIT))
        {
            // needs refactoring - this if isn't actually needed
            if (mUseAreaBasedViscosity)
            {
                damping_constantA *= (mpParams->GetDampingConstantMutant()/mpParams->GetDampingConstantNormal());
            }
            else
            {
                damping_constantA = mpParams->GetDampingConstantMutant();
            }
        }


        if(   (spring_iterator.rGetCellB().GetMutationState()!=HEALTHY)
           && (spring_iterator.rGetCellB().GetMutationState()!=APC_ONE_HIT))
        {
            // needs refactoring - this if isn't actually needed
            if (mUseAreaBasedViscosity)
            {
                damping_constantB *= (mpParams->GetDampingConstantMutant()/mpParams->GetDampingConstantNormal());
            }
            else
            {
                damping_constantB = mpParams->GetDampingConstantMutant();
            }
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
    c_vector<double, DIM> node_a_location = mrTissue.rGetMesh().GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = mrTissue.rGetMesh().GetNode(nodeBGlobalIndex)->rGetLocation();
    
    // there is reason not to substract one position from the other (cyclidrical meshes). clever gary
    unit_difference = mrTissue.rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);   
    
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
        
    double ageA = mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex).GetAge();
    double ageB = mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex).GetAge();
    
    TissueCell& r_cell_A = mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex);
    TissueCell& r_cell_B = mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex);
    
    if (ageA<CancerParameters::Instance()->GetMDuration() && ageB<CancerParameters::Instance()->GetMDuration() )
    {
        // Spring Rest Length Increases to normal rest length from ???? to normal rest length, 1.0, over 1 hour
        if (mrTissue.IsMarkedSpring(r_cell_A, r_cell_B))
        {   
            double lambda=CancerParameters::Instance()->GetDivisionRestingSpringLength();
            rest_length=(lambda+(1.0-lambda)*(ageA/(CancerParameters::Instance()->GetMDuration())));           
        }
        
        if (ageA+ SimulationTime::Instance()->GetTimeStep() >= CancerParameters::Instance()->GetMDuration())
        {
            // This spring is about to go out of scope
            mrTissue.UnmarkSpring(r_cell_A, r_cell_B);
        }
    }
    
    double a_rest_length=rest_length*0.5;
    double b_rest_length=a_rest_length;    
    
    if (mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex).HasApoptosisBegun())
    {
        double time_until_death_a = mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex).TimeUntilDeath();
        a_rest_length = a_rest_length*(time_until_death_a)/(CancerParameters::Instance()->GetApoptosisTime());
    }
    if (mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex).HasApoptosisBegun())
    {
        double time_until_death_b = mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex).TimeUntilDeath();
        b_rest_length = b_rest_length*(time_until_death_b)/(CancerParameters::Instance()->GetApoptosisTime());
    }
    
    rest_length = a_rest_length + b_rest_length;
    
    assert(rest_length<=1.0+1e-12);
    
    double multiplication_factor = 1.0;
    
    if (mUseEdgeBasedSpringConstant)
    {
        VoronoiTessellation<DIM>& tess = mrTissue.rGetVoronoiTessellation();
        
        multiplication_factor = tess.GetEdgeLength(nodeAGlobalIndex,nodeBGlobalIndex)*sqrt(3);
    }
    
    if (mUseMutantSprings)
    {
        unsigned number_of_mutants=0;
        
        if (r_cell_A.GetMutationState() == APC_TWO_HIT || r_cell_A.GetMutationState() == BETA_CATENIN_ONE_HIT)
        {   // if cell A is mutant
            number_of_mutants++;
        }
        
        if (r_cell_B.GetMutationState() == APC_TWO_HIT || r_cell_B.GetMutationState() == BETA_CATENIN_ONE_HIT)
        {   // if cell B is mutant
            number_of_mutants++;
        }
        
        switch  (number_of_mutants)
        {
            case 1u:
            {
                multiplication_factor*= mNormalMutantMultiplier;
                break;
            }
            case 2u:
            {
                multiplication_factor*= mMutantMutantMultiplier;
                break;
            }
        }
                    
    }
    return multiplication_factor * mpParams->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
}

//
//
//
//
//
//template<unsigned DIM> 
//c_vector<double, DIM> TissueSimulation<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex)
//{
//    assert(nodeAGlobalIndex!=nodeBGlobalIndex);
//    c_vector<double, DIM> unit_difference;
//    c_vector<double, DIM> node_a_location = mrTissue.rGetMesh().GetNode(nodeAGlobalIndex)->rGetLocation();
//    c_vector<double, DIM> node_b_location = mrTissue.rGetMesh().GetNode(nodeBGlobalIndex)->rGetLocation();
//    
//    // there is reason not to substract one position from the other (cyclidrical meshes). clever gary
//    unit_difference = mrTissue.rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);   
//    
//    double distance_between_nodes = norm_2(unit_difference);
//    
//    unit_difference /= distance_between_nodes;
//    
//    if(mUseCutoffPoint)
//    {
//        if( distance_between_nodes >= mCutoffPoint )
//        {
//            return zero_vector<double>(DIM);
//            //c_vector<double,DIM>() is not guaranteed to be fresh memory
//            // ie return zero force;
//        }
//    }
//    
//    double rest_length = 1.0;
//        
//    double ageA = mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex).GetAge();
//    double ageB = mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex).GetAge();
//    
//    if (ageA<CancerParameters::Instance()->GetMDuration() && ageB<CancerParameters::Instance()->GetMDuration() )
//    {
//        // Spring Rest Length Increases to normal rest length from ???? to normal rest length, 1.0, over 1 hour
//        TissueCell& r_cell_A = mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex);
//        TissueCell& r_cell_B = mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex);
//        if (mrTissue.IsMarkedSpring(r_cell_A, r_cell_B))
//        {   
//            double lambda=CancerParameters::Instance()->GetDivisionRestingSpringLength();
//            rest_length=(lambda+(1.0-lambda)*(ageA/(CancerParameters::Instance()->GetMDuration())));           
//        }
//        
//        if (ageA+ SimulationTime::Instance()->GetTimeStep() >= CancerParameters::Instance()->GetMDuration())
//        {
//            // This spring is about to go out of scope
//            mrTissue.UnmarkSpring(r_cell_A, r_cell_B);
//        }
//    }
//    
//    double a_rest_length=rest_length*0.5;
//    double b_rest_length=a_rest_length;    
//    
//    if (mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex).HasApoptosisBegun())
//    {
//        double time_until_death_a = mrTissue.rGetCellAtNodeIndex(nodeAGlobalIndex).TimeUntilDeath();
//        a_rest_length = a_rest_length*(time_until_death_a)/(CancerParameters::Instance()->GetApoptosisTime());
//    }
//    if (mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex).HasApoptosisBegun())
//    {
//        double time_until_death_b = mrTissue.rGetCellAtNodeIndex(nodeBGlobalIndex).TimeUntilDeath();
//        b_rest_length = b_rest_length*(time_until_death_b)/(CancerParameters::Instance()->GetApoptosisTime());
//    }
//    
//    rest_length = a_rest_length + b_rest_length;
//    
//    assert(rest_length<=1.0+1e-12);
//    
//    double multiplication_factor = 1.0;
//    
//    if (mUseEdgeBasedSpringConstant)
//    {
//        VoronoiTessellation<DIM>& tess = mrTissue.rGetVoronoiTessellation();
//        
//        multiplication_factor = tess.GetEdgeLength(nodeAGlobalIndex,nodeBGlobalIndex)*sqrt(3);
//    }
//    
//    return multiplication_factor * mpParams->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
//}
//


template<unsigned DIM> 
void TissueSimulation<DIM>::UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rDrDt)
{
    // update ghost positions first because they do not affect the real cells
    mrTissue.UpdateGhostPositions(mDt);

    // Iterate over all cells to update their positions.
    for (typename Tissue<DIM>::Iterator cell_iter = mrTissue.Begin();
         cell_iter != mrTissue.End();
         ++cell_iter)
    {
        TissueCell& cell = *cell_iter;
        unsigned index = cell.GetNodeIndex();
        
        ChastePoint<DIM> new_point(mrTissue.rGetMesh().GetNode(index)->rGetLocation() + mDt*rDrDt[index]);
        mrTissue.MoveCell(cell_iter, new_point);    
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
    if (maxCells<mrTissue.rGetMesh().GetNumAllNodes())
    {
        EXCEPTION("mMaxCells is less than the number of cells in the mesh.");
    }
    mrTissue.SetMaxCells(maxCells);
}

/**
 * Sets the maximum number of elements that the simulation will contain (for use by the datawriter)
 * default value is set to 10x the initial mesh value by the constructor.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetMaxElements(unsigned maxElements)
{
    mMaxElements = maxElements;
    if (maxElements<mrTissue.rGetMesh().GetNumAllElements())
    {
        EXCEPTION("mMaxElements is less than the number of elements in the mesh.");
    }
    mrTissue.SetMaxElements(maxElements);
}


template<unsigned DIM> 
Tissue<DIM>& TissueSimulation<DIM>::rGetTissue()
{
    return mrTissue;
}

template<unsigned DIM> 
const Tissue<DIM>& TissueSimulation<DIM>::rGetTissue() const
{
    return mrTissue;
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
void TissueSimulation<DIM>::SetOutputCellTypes(bool outputCellTypes)
{
    mOutputCellTypes = outputCellTypes;
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
 * Use an edge-based spring constant
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant)
{
    assert(DIM == 2);
    mUseEdgeBasedSpringConstant = useEdgeBasedSpringConstant;
    mCreateVoronoiTessellation = useEdgeBasedSpringConstant;
}

/**
 * Use an area based viscosity
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetAreaBasedViscosity(bool useAreaBasedViscosity)
{
    assert(DIM == 2);
    mUseAreaBasedViscosity = useAreaBasedViscosity;
    mCreateVoronoiTessellation = useAreaBasedViscosity;
}


/**
 * Use Different spring strengths depending on two cells:
 * Normal-normal, Normal-mutant, mutant-mutant
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetMutantSprings(bool useMutantSprings, double mutantMutantMultiplier=2, double normalMutantMultiplier=1.5)
{
    mUseMutantSprings = useMutantSprings;
    mMutantMutantMultiplier = mutantMutantMultiplier;
    mNormalMutantMultiplier = normalMutantMultiplier;
}


template<unsigned DIM> 
void TissueSimulation<DIM>::SetWriteVoronoiData(bool writeVoronoiData, bool followLoggedCell)
{
    assert(DIM == 2);
    mWriteVoronoiData = writeVoronoiData;
    mCreateVoronoiTessellation = writeVoronoiData;
    mFollowLoggedCell = followLoggedCell;
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
        location.push_back( mrTissue.rGetMesh().GetNode(rNodeIndex)->rGetLocation()[i] );
    }
    return location;
}

/**
 * Main Solve method.
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
    
    mrTissue.SetupTabulatedWriters(tabulated_node_writer, tabulated_element_writer);//, element_writer_ids);
    
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
    for (typename Tissue<DIM>::Iterator cell_iter = mrTissue.Begin();
         cell_iter != mrTissue.End();
         ++cell_iter)
    {
        /* We don't use the result; this call is just to force the cells 
         * to age to current time running their cell cycle models to get there.
         */
        cell_iter->ReadyToDivide();
    }
        
    // Write initial conditions to file for the visualizer.
    if(DIM==2)
    {
        WriteVisualizerSetupFile(*p_setup_file);
    }
    
    mrTissue.WriteResultsToFiles(tabulated_node_writer, 
                                tabulated_element_writer,
                                *p_node_file, *p_element_file, *p_cell_types_file,
                                false,
                                true,
                                mOutputCellTypes);

    CryptVoronoiDataWriter<DIM>* p_voronoi_data_writer = NULL;
    if(mWriteVoronoiData)
    {
        p_voronoi_data_writer = new CryptVoronoiDataWriter<DIM>(mrTissue,
                                                                mOutputDirectory,
                                                                "VoronoiAreaAndPerimeter.dat",
                                                                mFollowLoggedCell);
    }

                               
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
            mrTissue.ReMesh();
        }

        if (mCreateVoronoiTessellation)
        {
            mrTissue.CreateVoronoiTessellation();
        }

        //  calculate node velocities
        std::vector<c_vector<double, DIM> > drdt = CalculateVelocitiesOfEachNode();

        // update node positions
        UpdateNodePositions(drdt);
     
        // Increment simulation time here, so results files look sensible
        p_simulation_time->IncrementTimeOneStep();
        
        // Write results to file
        mrTissue.WriteResultsToFiles(tabulated_node_writer, 
                                    tabulated_element_writer, 
                                    *p_node_file, *p_element_file, *p_cell_types_file,
                                    tabulated_output_counter%80==0,
                                    true,
                                    mOutputCellTypes);
                                    
        if(mWriteVoronoiData)
        {
            p_voronoi_data_writer->WriteData();
        }

        tabulated_output_counter++;
        
        PostSolve();
    }
    
    // Write end state to tabulated files (not visualizer - this
    // is taken care of in the main loop).
    // Doesn't need to count cell types again as it is done in the last loop
    mrTissue.WriteResultsToFiles(tabulated_node_writer, 
                                tabulated_element_writer, 
                                *p_node_file, *p_element_file, *p_cell_types_file,
                                true,
                                false,
                                false);
                        
    tabulated_node_writer.Close();
    tabulated_element_writer.Close();
    
    if(p_voronoi_data_writer!=NULL)
    {
        delete p_voronoi_data_writer;
    }
}



/**
 * Saves the whole tissue simulation for restarting later.
 *
 * Puts it in the folder mOutputDirectory/archive/
 * and the file "tissue_sim_at_time_<SIMULATION TIME>.arch"
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
    std::string archive_filename = handler.GetOutputDirectoryFullPath() + "tissue_sim_at_time_"+time_stamp.str()+".arch";
    std::string mesh_filename = std::string("mesh_") + time_stamp.str();
    
    if(mReMesh)
    {
        mrTissue.ReMesh();
    }

    
    // the false is so the directory isn't cleaned
    TrianglesMeshWriter<DIM,DIM> mesh_writer(archive_directory, mesh_filename, false);
    mesh_writer.WriteFilesUsingMesh(mrTissue.rGetMesh());
    
    std::ofstream ofs(archive_filename.c_str());
    boost::archive::text_oarchive output_arch(ofs);
    
    // cast to const.
    const SimulationTime* p_simulation_time = SimulationTime::Instance();
    output_arch << *p_simulation_time;
    TissueSimulation<DIM> * p_sim = this;
    output_arch & p_sim;
}

/**
 * Loads a saved tissue simulation to run further.
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
    
    std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
    
    std::string archive_filename = test_output_directory + rArchiveDirectory + "/archive/tissue_sim_at_time_"+time_stamp.str() +".arch";
    std::string mesh_filename = test_output_directory + rArchiveDirectory + "/archive/mesh_" + time_stamp.str();
    Tissue<DIM>::meshPathname = mesh_filename;
    
    // Create an input archive
    std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
    boost::archive::text_iarchive input_arch(ifs);
        
    // Read the archive
    assert(p_simulation_time->IsStartTimeSetUp());
    input_arch >> *p_simulation_time;

    TissueSimulation<DIM>* p_sim;

    input_arch >> p_sim;

    if (p_sim->rGetTissue().rGetMesh().GetNumNodes()!=p_sim->rGetTissue().rGetCells().size())
    {
        #define COVERAGE_IGNORE
        std::stringstream string_stream;
        string_stream << "Error in Load(), number of nodes (" << p_sim->rGetTissue().rGetMesh().GetNumNodes()
                      << ") is not equal to the number of cells (" << p_sim->rGetTissue().rGetCells().size() 
                      << ")";
        EXCEPTION(string_stream.str());
        #undef COVERAGE_IGNORE
    }
    
    return p_sim;
}

/**
 * Find out how many cells of each mutation state there are
 * 
 * @return The number of cells of each type (evaluated at each visualizer output)
 * [0] = healthy count
 * [1] = labelled cells
 * [2] = APC one hit
 * [3] = APC two hit
 * [4] = beta catenin one hit
 */
template<unsigned DIM>
c_vector<unsigned,5> TissueSimulation<DIM>::GetCellTypeCount()
{
    return mrTissue.GetCellTypeCount();
}

#endif //_TISSUESIMULATION_CPP_
