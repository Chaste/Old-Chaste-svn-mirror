#ifndef _TISSUESIMULATION_CPP_
#define _TISSUESIMULATION_CPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <set>
#include "TissueSimulation.hpp"
#include "SimpleTissueMechanicsSystem.hpp"
#include "Exception.hpp"
#include "RandomNumberGenerator.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.cpp"
#include "WntConcentration.hpp"
#include "CellwiseData.cpp"
#include "LogFile.hpp"


template<unsigned DIM> 
TissueSimulation<DIM>::TissueSimulation(AbstractTissue<DIM>& rTissue, 
                                        AbstractDiscreteTissueMechanicsSystem<DIM>* pMechanicsSystem, 
                                        bool deleteTissueAndMechanicsSystem,
                                        bool initialiseCells)
  :  mrTissue(rTissue)
{
    #define COVERAGE_IGNORE
    assert(DIM==2 || DIM==3); // there are no instances of TissueSimulation<1>
    #undef COVERAGE_IGNORE        
    
    CancerEventHandler::BeginEvent(CANCER_EVERYTHING);

    mDeleteTissue = deleteTissueAndMechanicsSystem;
    mInitialiseCells = initialiseCells;
    
    mpParams = CancerParameters::Instance();
    // This line sets a random seed of 0 if it wasn't specified earlier.
    mpRandomGenerator = RandomNumberGenerator::Instance();
    
    mDt = 1.0/120.0; // Timestep of 30 seconds (as per Meineke)
    mEndTime = 0.0; // hours - this is set later on.
    
    // Defaults
    mOutputDirectory = "";
    mSimulationOutputDirectory = mOutputDirectory;
    
    if (mrTissue.HasMesh())
    {
        mReMesh = true;
    }
    else
    {
        mReMesh = false;
    }
        
    mOutputCellMutationStates = false;
    mOutputCellAncestors = false;
    mOutputCellTypes = false;
    mOutputCellVariables = false;
    mOutputCellCyclePhases = false;
    mNoBirth = false;
    mNumBirths = 0;
    mNumDeaths = 0;
    mSamplingTimestepMultiple = 1;
    mAllocatedMemoryForMechanicsSystem = deleteTissueAndMechanicsSystem;
    
    if (pMechanicsSystem == NULL)
    {
        mAllocatedMemoryForMechanicsSystem = true;
        if (mrTissue.HasMesh())
        {
            pMechanicsSystem = new Meineke2001SpringSystem<DIM>(*(static_cast<MeshBasedTissue<DIM>*>(&mrTissue)));
        }
        else
        {
            pMechanicsSystem = new SimpleTissueMechanicsSystem<DIM>(*(static_cast<SimpleTissue<DIM>*>(&mrTissue)));
        }
    }
    mpMechanicsSystem = pMechanicsSystem;    
    
    if (mInitialiseCells)
    {
        mrTissue.InitialiseCells();
    }
}

/**
 * Free any memory allocated by the constructor.
 * This frees the tissue and cell killers, if they were created by de-serialization.
 */
template<unsigned DIM> 
TissueSimulation<DIM>::~TissueSimulation()
{
    if (mAllocatedMemoryForMechanicsSystem)
    {
        delete mpMechanicsSystem;
    }
    if (mDeleteTissue)
    {
        for (typename std::vector<AbstractCellKiller<DIM>*>::iterator it=mCellKillers.begin();
             it != mCellKillers.end();
             ++it)
        {
            delete *it;
        }
        delete &mrTissue;
    }
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
    for (typename AbstractTissue<DIM>::Iterator cell_iter = mrTissue.Begin();
         cell_iter != mrTissue.End();
         ++cell_iter)
    {
        TissueCell& cell = *cell_iter;

        // Check if this cell is ready to divide - if so create a new cell etc.
        if (cell.GetAge() > 0.0)
        {
            if (cell.ReadyToDivide())
            {
                // Create new cell
                TissueCell new_cell = cell.Divide();
            
                // Add a new node to the mesh
                c_vector<double, DIM> new_location = CalculateDividingCellCentreLocations(cell_iter);
                
                TissueCell *p_new_cell = mrTissue.AddCell(new_cell, new_location);
                
                if (mrTissue.HasMesh())
                {
                    (static_cast<MeshBasedTissue<DIM>*>(&mrTissue))->MarkSpring(cell, *p_new_cell);
                }
                num_births_this_step++;
            } 
        }
    } 
   
    return num_births_this_step;
}

template<unsigned DIM> 
unsigned TissueSimulation<DIM>::DoCellRemoval()
{
    unsigned num_deaths_this_step=0;
        
    // This labels cells as dead or apoptosing. It does not actually remove the cells, 
    // tissue.RemoveDeadCells() needs to be called for this.
    for(unsigned killer_index = 0; killer_index<mCellKillers.size(); killer_index++)
    {
        mCellKillers[killer_index]->TestAndLabelCellsForApoptosisOrDeath();
    }
    
    num_deaths_this_step += mrTissue.RemoveDeadCells(); 
    
    return num_deaths_this_step;
}

template<unsigned DIM> 
c_vector<double, DIM> TissueSimulation<DIM>::CalculateDividingCellCentreLocations(typename AbstractTissue<DIM>::Iterator parentCell)
{
    double separation = CancerParameters::Instance()->GetDivisionSeparation();
    c_vector<double, DIM> parent_coords = parentCell.rGetLocation();
    c_vector<double, DIM> daughter_coords;
        
    // Pick a random direction and move the parent cell backwards by 0.5*sep in that
    // direction and return the position of the daughter cell (0.5*sep forwards in the
    // random vector direction

    // Make a random direction vector of the required length
    c_vector<double, DIM> random_vector;
    
    if (DIM==2)
    {
        double random_angle = RandomNumberGenerator::Instance()->ranf();
        random_angle *= 2.0*M_PI;
        
        random_vector(0) = 0.5*separation*cos(random_angle);
        random_vector(1) = 0.5*separation*sin(random_angle);
        
        parent_coords = parent_coords-random_vector;
        daughter_coords = parent_coords+random_vector;        
    }
    else if (DIM==3)
    {
        double random_zenith_angle = RandomNumberGenerator::Instance()->ranf(); // phi 
        random_zenith_angle *= M_PI;
        double random_azimuth_angle = RandomNumberGenerator::Instance()->ranf(); // theta
        random_azimuth_angle *= 2*M_PI;
        
        random_vector(0) = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
        random_vector(1) = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
        random_vector(2) = 0.5*separation*cos(random_zenith_angle);

        daughter_coords = parent_coords+random_vector;
        parent_coords = parent_coords-random_vector;
    }
        
    // Set the parent to use this location
    ChastePoint<DIM> parent_coords_point(parent_coords);
    mrTissue.MoveCell(parentCell, parent_coords_point);
    return daughter_coords;
}

template<unsigned DIM> 
void TissueSimulation<DIM>::UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rDrDt)
{
    if (mrTissue.HasGhostNodes())
    {
        // Update ghost positions first because they do not affect the real cells
        (static_cast<MeshBasedTissueWithGhostNodes<DIM>*>(&mrTissue))->UpdateGhostPositions(mDt);
    }

    // Iterate over all cells to update their positions.
    for (typename AbstractTissue<DIM>::Iterator cell_iter = mrTissue.Begin();
         cell_iter != mrTissue.End();
         ++cell_iter)
    {
        TissueCell& cell = *cell_iter;
        unsigned index = cell.GetNodeIndex();
        
        ChastePoint<DIM> new_point(mrTissue.GetNode(index)->rGetLocation() + mDt*rDrDt[index]);
        mrTissue.MoveCell(cell_iter, new_point);    
    }
}

/**
 * Set the timestep of the simulation
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetDt(double dt)
{
    assert(dt > 0);
    mDt = dt;
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
    assert(endTime > 0);
    mEndTime = endTime;
}

/**
 * Set the output directory of the simulation.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
    mSimulationOutputDirectory = mOutputDirectory;
}
    
/**
 * Sets the ratio of the number of actual timesteps to the number of timesteps 
 * at which results are written to file. Default value is set to 1 by the constructor.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetSamplingTimestepMultiple(unsigned samplingTimestepMultiple)
{
    assert(samplingTimestepMultiple > 0);
    mSamplingTimestepMultiple = samplingTimestepMultiple;
}

template<unsigned DIM> 
AbstractTissue<DIM>& TissueSimulation<DIM>::rGetTissue()
{
    return mrTissue;
}

template<unsigned DIM> 
const AbstractTissue<DIM>& TissueSimulation<DIM>::rGetTissue() const
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
 * Set the simulation to count and store the number of each cell mutation state.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetOutputCellMutationStates(bool outputCellMutationStates)
{
    mOutputCellMutationStates = outputCellMutationStates;
}

/**
 * Set the simulation to count and store the number of each cell type.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetOutputCellTypes(bool outputCellTypes)
{
    mOutputCellTypes = outputCellTypes;
}

/**
 * Set the simulation to output the cell ancestors if they are set.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetOutputCellAncestors(bool outputCellAncestors)
{
    mOutputCellAncestors = outputCellAncestors;
}

/**
 * Set the simulation to output the cell variables.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetOutputCellVariables(bool outputCellVariables)
{
    mOutputCellVariables = outputCellVariables;
}

/**
 * Set the simulation to output the cell cycle phases.
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::SetOutputCellCyclePhases(bool outputCellCyclePhases)
{
    mOutputCellCyclePhases = outputCellCyclePhases;
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
        location.push_back( mrTissue.GetNode(rNodeIndex)->rGetLocation()[i] );
    }
    return location;
}

/**
 * Main Solve method
 */
template<unsigned DIM> 
void TissueSimulation<DIM>::Solve()
{        
    CancerEventHandler::BeginEvent(SETUP);

    // Set up the simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetDimensionalisedTime();
    
    unsigned num_time_steps = (unsigned) ((mEndTime-current_time)/mDt+0.5);

    if (current_time > 0) // use the reset function if necessary
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
    mSimulationOutputDirectory = results_directory;
        
    ///////////////////////////////////////////////////////////
    // Set up Simulation
    ///////////////////////////////////////////////////////////
   
    // Create output files for the visualizer
    OutputFileHandler output_file_handler(results_directory+"/", false);
    
    mrTissue.CreateOutputFiles(results_directory+"/",
                               false,
                               mOutputCellMutationStates,
                               mOutputCellTypes,
                               mOutputCellVariables,
                               mOutputCellCyclePhases,
                               mOutputCellAncestors);
    
    mpSetupFile = output_file_handler.OpenOutputFile("results.vizsetup");

    SetupSolve();
    
    // Age the cells to the correct time. Note that cells are created with
    // negative birth times so that some are initially almost ready to divide.
    LOG(1, "Setting up cells...");
    for (typename AbstractTissue<DIM>::Iterator cell_iter = mrTissue.Begin();
         cell_iter != mrTissue.End();
         ++cell_iter)
    {
        // We don't use the result; this call is just to force the cells to age 
        // to the current time running their cell cycle models to get there
        cell_iter->ReadyToDivide();
    }
    LOG(1, "\tdone\n");
     
    // Write initial conditions to file for the visualizer
    if (DIM==2)
    {
        WriteVisualizerSetupFile();
    }
    mpSetupFile->close();
    
    mrTissue.WriteResultsToFiles(mOutputCellMutationStates,
                                 mOutputCellTypes,
                                 mOutputCellVariables,
                                 mOutputCellCyclePhases,
                                 mOutputCellAncestors);

    CancerEventHandler::EndEvent(SETUP);
                               
    /////////////////////////////////////////////////////////////////////
    // Main time loop
    /////////////////////////////////////////////////////////////////////
    while (p_simulation_time->GetTimeStepsElapsed() < num_time_steps)
    {          
        LOG(1, "--TIME = " << p_simulation_time->GetDimensionalisedTime() << "\n");
        
        // Remove dead cells then implement cell birth. Note that neither
        // of these methods use any element information, they just delete
        // and create nodes
        CancerEventHandler::BeginEvent(DEATH);
        mNumDeaths += DoCellRemoval();
        LOG(1, "\tNum deaths = " << mNumDeaths << "\n");
        CancerEventHandler::EndEvent(DEATH);

        CancerEventHandler::BeginEvent(BIRTH);
        mNumBirths += DoCellBirth();
        LOG(1, "\tNum births = " << mNumBirths << "\n");
        CancerEventHandler::EndEvent(BIRTH);
        
        // If the tissue has a mesh, then we currently must call a ReMesh at
        // each timestep. Otherwise, we only need to call a ReMesh after there 
        // has been any cell birth or cell death.        
        if (mrTissue.HasMesh())
        {
            assert(mReMesh);
        }
        else
        {
            mReMesh = false;
            if ( (mNumBirths>0) || (mNumDeaths>0) )
            {   
                mReMesh = true;
            }
        }
        
        // Remesh
        CancerEventHandler::BeginEvent(REMESH);
        if (mReMesh)
        {
            LOG(1, "\tRemeshing...");
            mrTissue.ReMesh();
            LOG(1, "\tdone.\n");
        }
        CancerEventHandler::EndEvent(REMESH);
        
        CancerEventHandler::BeginEvent(TESSELLATION);
        if (mrTissue.HasMesh())
        {
            if ( (static_cast<MeshBasedTissue<DIM>*>(&mrTissue))->GetWriteVoronoiData() 
                 || mpMechanicsSystem->NeedsVoronoiTessellation() 
                 || (static_cast<MeshBasedTissue<DIM>*>(&mrTissue))->GetWriteTissueAreas() )
            {
                (static_cast<MeshBasedTissue<DIM>*>(&mrTissue))->CreateVoronoiTessellation();
            }
        }
        CancerEventHandler::EndEvent(TESSELLATION);

        // Calculate node velocities
        CancerEventHandler::BeginEvent(VELOCITY);
        std::vector<c_vector<double, DIM> >& drdt = mpMechanicsSystem->rCalculateVelocitiesOfEachNode();
        CancerEventHandler::EndEvent(VELOCITY);

        // Update node positions
        CancerEventHandler::BeginEvent(POSITION);
        UpdateNodePositions(drdt);
        CancerEventHandler::EndEvent(POSITION);
     
        PostSolve();
                          
        // Increment simulation time here, so results files look sensible
        p_simulation_time->IncrementTimeOneStep();
        
        CancerEventHandler::BeginEvent(OUTPUT);
                
        // Write results to file
        if (p_simulation_time->GetTimeStepsElapsed()%mSamplingTimestepMultiple==0)
        {
            mrTissue.WriteResultsToFiles(mOutputCellMutationStates,
                                         mOutputCellTypes,
                                         mOutputCellVariables,
                                         mOutputCellCyclePhases,
                                         mOutputCellAncestors);
        }
        
        CancerEventHandler::EndEvent(OUTPUT);
    }

    AfterSolve();
    
    CancerEventHandler::BeginEvent(OUTPUT);
    mrTissue.CloseOutputFiles(mOutputCellMutationStates,
                              mOutputCellTypes,
                              mOutputCellVariables,
                              mOutputCellCyclePhases,
                              mOutputCellAncestors);
                              
    CancerEventHandler::EndEvent(OUTPUT);
    
    CancerEventHandler::EndEvent(CANCER_EVERYTHING);
}

template<unsigned DIM>
void TissueSimulation<DIM>::AfterSolve()
{
    LOG(1, "--TIME = " << SimulationTime::Instance()->GetDimensionalisedTime() << "\n");
    
    // Remove dead cells then implement cell birth
    CancerEventHandler::BeginEvent(DEATH);
    mNumDeaths += DoCellRemoval();
    LOG(1, "\tNum deaths = " << mNumDeaths << "\n");
    CancerEventHandler::EndEvent(DEATH);

    CancerEventHandler::BeginEvent(BIRTH);
    mNumBirths += DoCellBirth();
    LOG(1, "\tNum births = " << mNumBirths << "\n");
    CancerEventHandler::EndEvent(BIRTH);
    
    // Carry out a final remesh if necessary
    if (mrTissue.HasMesh())
    {
        assert(mReMesh);
    }
    else
    {
        mReMesh = false;
        if ( (mNumBirths>0) || (mNumDeaths>0) )
        {   
            mReMesh = true;
        }
    }
    
    if (mReMesh)
    {
        LOG(1, "\tRemeshing...");
        mrTissue.ReMesh();
        LOG(1, "\tdone.\n");
    }
}

/**
 * Saves the whole tissue simulation for restarting later.
 *
 * Puts it in the folder mOutputDirectory/archive/
 * and the file "tissue_sim_at_time_<SIMULATION TIME>.arch".
 * The mesh is written to files in the same folder.
 *
 * First archives simulation time then the simulation itself.
 */
template<unsigned DIM>
void TissueSimulation<DIM>::Save()
{
    CommonSave(this);
}

/**
 * The function that does the actual work.  Templated over the type
 * of simulation that's actually being saved, since normal inheritance
 * doesn't seem to work for this - we get segfaults on loading.
 * Not sure why this is so, but...
 *
 * @param pSim = this, but with the type explicit
 */
template<unsigned DIM>
template<class SIM>
void TissueSimulation<DIM>::CommonSave(SIM* pSim)
{
    // Get the simulation time as a string
    const SimulationTime* p_sim_time = SimulationTime::Instance();
    assert(p_sim_time->IsStartTimeSetUp());
    std::ostringstream time_stamp;
    time_stamp << p_sim_time->GetDimensionalisedTime();
    
    // Create an output file handler in order to get the full path of the
    // archive directory.  Note the false is so the handler doesn't clean
    // the directory.
    std::string archive_directory = mOutputDirectory + "/archive/";
    OutputFileHandler handler(archive_directory, false);
    std::string archive_filename = handler.GetOutputDirectoryFullPath() + "tissue_sim_at_time_"+time_stamp.str()+".arch";
    std::string mesh_filename = std::string("mesh_") + time_stamp.str();

    // Write the mesh to file.  Remesh first to ensure it's in a good state.
    mrTissue.ReMesh();   
    
    if (mrTissue.HasMesh())
    {
        // The false is so the directory isn't cleaned
        TrianglesMeshWriter<DIM,DIM> mesh_writer(archive_directory, mesh_filename, false);
        mesh_writer.WriteFilesUsingMesh((static_cast<MeshBasedTissue<DIM>*>(&mrTissue))->rGetMesh());
    }
    
    // Create a new archive
    std::ofstream ofs(archive_filename.c_str());
    boost::archive::text_oarchive output_arch(ofs);
    
    // Save the simulation.  We save the time directly first to maintain its
    // singleton-ness on load.
    output_arch << *p_sim_time;
    
    // Archive the Wnt concentration if it's used
    bool archive_wnt = WntConcentration::Instance()->IsWntSetUp();
    output_arch & archive_wnt;
    if (archive_wnt)
    {
        WntConcentration* p_wnt = WntConcentration::Instance();
        output_arch & *p_wnt;
    }
    
    // Archive the CellwiseData if it's used
    bool archive_cellwise_data = CellwiseData<DIM>::Instance()->IsSetUp();
    output_arch & archive_cellwise_data;
    if (archive_cellwise_data)
    {
        CellwiseData<DIM>* p_cellwise_data = CellwiseData<DIM>::Instance();
        output_arch & *p_cellwise_data;
    }    
    output_arch & pSim; // const-ness would be a pain here
}

/**
 * Loads a saved tissue simulation to run further.
 *
 * @param rArchiveDirectory the name of the simulation to load
 * (specified originally by simulation.SetOutputDirectory("wherever"); )
 * @param rTimeStamp the time at which to load the simulation (this must
 * be one of the times at which simulation.Save() was called)
 */
template<unsigned DIM> 
TissueSimulation<DIM>* TissueSimulation<DIM>::Load(const std::string& rArchiveDirectory, const double& rTimeStamp)
{    
    std::string archive_filename = TissueSimulation<DIM>::GetArchivePathname(rArchiveDirectory, rTimeStamp);

    // Create an input archive
    std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
    boost::archive::text_iarchive input_arch(ifs);
    
    TissueSimulation<DIM>::CommonLoad(input_arch);
    
    TissueSimulation<DIM>* p_sim;
    input_arch >> p_sim;

    if (p_sim->rGetTissue().GetNumNodes()!=p_sim->rGetTissue().rGetCells().size())
    {
        #define COVERAGE_IGNORE
        std::stringstream string_stream;
        string_stream << "Error in Load(), number of nodes (" << p_sim->rGetTissue().GetNumNodes()
                      << ") is not equal to the number of cells (" << p_sim->rGetTissue().rGetCells().size() 
                      << ")";
        EXCEPTION(string_stream.str());
        #undef COVERAGE_IGNORE
    }

    return p_sim;
}

/**
 * Find the right archive (and mesh) to load.  The files are contained within
 * the 'archive' folder in rArchiveDirectory, with the archive itself called
 * 'tissue_sim_at_time_`rTimeStamp`.arch'.  The path to this file is returned.
 *
 * The path to the mesh is stored as Tissue<DIM>::meshPathname for use by the
 * Tissue de-serialization routines.
 */
template<unsigned DIM>
std::string TissueSimulation<DIM>::GetArchivePathname(const std::string& rArchiveDirectory, const double& rTimeStamp)
{
    // Find the right archive and mesh to load
    std::ostringstream time_stamp;
    time_stamp << rTimeStamp;
    
    std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
    
    std::string archive_filename = test_output_directory + rArchiveDirectory + "/archive/tissue_sim_at_time_"+time_stamp.str() +".arch";
    std::string mesh_filename = test_output_directory + rArchiveDirectory + "/archive/mesh_" + time_stamp.str();
    MeshBasedTissue<DIM>::meshPathname = mesh_filename;
    return archive_filename;
}


/**
 * Load any data from the archive that isn't the simulation class itself.
 * At present this is just the simulation time.
 */
template<unsigned DIM>
template<class Archive>
void TissueSimulation<DIM>::CommonLoad(Archive& rInputArch)
{
    // Load simulation time
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    assert(p_simulation_time->IsStartTimeSetUp());
    rInputArch >> *p_simulation_time;
    
    // Load Wnt concentration if it's used
    bool archive_wnt;
    rInputArch & archive_wnt;
    if (archive_wnt)
    {
        WntConcentration* p_wnt = WntConcentration::Instance();
        rInputArch & *p_wnt;
    }
    
    // Load CellwiseData if it's used
    bool archive_cellwise_data;
    rInputArch & archive_cellwise_data;
    if (archive_cellwise_data)
    {
        CellwiseData<DIM>* p_cellwise_data = CellwiseData<DIM>::Instance();
        rInputArch & *p_cellwise_data;
    }
}

/**
 * Find out how many cells of each mutation state there are
 * 
 * @return The number of cells of each mutation state (evaluated at each visualizer output)
 * [0] = healthy count
 * [1] = labelled cells
 * [2] = APC one hit
 * [3] = APC two hit
 * [4] = beta catenin one hit
 */
template<unsigned DIM>
c_vector<unsigned, NUM_CELL_MUTATION_STATES> TissueSimulation<DIM>::GetCellMutationStateCount()
{
    return mrTissue.GetCellMutationStateCount();
}

/**
 * Find out how many cells of each type there are
 * 
 * @return The number of cells of each type (evaluated at each visualizer output)
 * [0] = STEM
 * [1] = TRANSIT
 * [2] = DIFFERENTIATED
 * [3] = NECROTIC
 */
template<unsigned DIM>
c_vector<unsigned, NUM_CELL_TYPES> TissueSimulation<DIM>::GetCellTypeCount()
{
    return mrTissue.GetCellTypeCount();
}

/**
 * Find out how many cells in each cell cycle phase there are
 * 
 * @return The number of cells of each phase (evaluated at each visualizer output)
 * [0] = G_ZERO_PHASE
 * [1] = G_ONE_PHASE
 * [2] = S_PHASE
 * [3] = G_TWO_PHASE
 * [4] = M_PHASE
 */
template<unsigned DIM>
c_vector<unsigned, 5> TissueSimulation<DIM>::GetCellCyclePhaseCount()
{
    return mrTissue.GetCellCyclePhaseCount();
}


#endif //_TISSUESIMULATION_CPP_
