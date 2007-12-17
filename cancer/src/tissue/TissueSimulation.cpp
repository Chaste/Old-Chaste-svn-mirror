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
#include "CellwiseData.cpp"

template<unsigned DIM> 
TissueSimulation<DIM>::TissueSimulation(Tissue<DIM>& rTissue, 
                                        AbstractDiscreteTissueMechanicsSystem<DIM>* pMechanicsSystem, 
                                        bool deleteTissue,
                                        bool initialiseCells)
  :  mrTissue(rTissue)
{
    #define COVERAGE_IGNORE
    assert(DIM==2 || DIM==3); // there are no instances of TissueSimulation<1>
    #undef COVERAGE_IGNORE
        
    
    CancerEventHandler::BeginEvent(CANCER_EVERYTHING);

    mDeleteTissue = deleteTissue;
    mInitialiseCells = initialiseCells;
    
    mpParams = CancerParameters::Instance();
    // this line sets a random seed of 0 if it wasn't specified earlier.
    mpRandomGenerator = RandomNumberGenerator::Instance();
    
    mDt = 1.0/120.0; // Timestep of 30 seconds (as per Meineke)
    mEndTime = 0.0; // hours - this is set later on.
    
    // defaults
    mOutputDirectory = "";
    mSimulationOutputDirectory = mOutputDirectory;
    mReMesh = true;
    mOutputCellTypes = false ;
    mNoBirth = false;
    mMaxCells = 10*mrTissue.rGetMesh().GetNumNodes();
    mMaxElements = 10*mrTissue.rGetMesh().GetNumElements();
    mNumBirths = 0;
    mNumDeaths = 0;

    mWriteVoronoiData = false;
    mFollowLoggedCell = false;
        
    mrTissue.SetMaxCells(mMaxCells);
    mrTissue.SetMaxElements(mMaxElements);
    
    if (pMechanicsSystem == NULL)
    {
        pMechanicsSystem = new Meineke2001SpringSystem<DIM>(mrTissue);
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
    delete mpMechanicsSystem;
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
        if (cell.GetAge()>0.0)
        {
            if (cell.ReadyToDivide())
            {
                // Create new cell
                TissueCell new_cell = cell.Divide();
            
                // Add a new node to the mesh
                c_vector<double, DIM> new_location = CalculateDividingCellCentreLocations(cell_iter);
                
                TissueCell *p_new_cell=mrTissue.AddCell(new_cell, new_location);
                mrTissue.MarkSpring(cell, *p_new_cell);
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
        
    // Pick a random direction and move the parent cell backwards by 0.5*sep in that
    // direction and return the position of the daughter cell (0.5*sep forwards in the
    // random vector direction

    // Make a random direction vector of the required length
    c_vector<double, DIM> random_vector;
    
//    if(DIM==1)
//    {
//        random_vector(0) = 0.5*separation;
//
//        daughter_coords = parent_coords+random_vector;
//        parent_coords = parent_coords-random_vector;
//    }   
//    else if(DIM==2)
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
    mSimulationOutputDirectory = mOutputDirectory;
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


template<unsigned DIM> 
void TissueSimulation<DIM>::SetWriteVoronoiData(bool writeVoronoiData, bool followLoggedCell)
{
    assert(DIM == 2);
    mWriteVoronoiData = writeVoronoiData;
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
    CancerEventHandler::BeginEvent(SETUP);

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
    mSimulationOutputDirectory = results_directory;
    
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
    mpSetupFile = output_file_handler.OpenOutputFile("results.vizsetup");
    
    // Creates output file to store number of different cells
    out_stream p_cell_types_file = output_file_handler.OpenOutputFile("celltypes.dat");
    if (mOutputCellTypes)
    { 
        *p_cell_types_file <<   "Time\t Healthy\t Labelled\t APC_1\t APC_2\t BETA_CAT \n";
    }
    
    SetupSolve();
    /* 
     * Age the cells to the correct time (cells set up with negative birth dates
     * to give some that are almost ready to divide).
     * 
     * TODO:For some strange reason this seems to take about 3 minutes for a realistic Wnt-Crypt.
     * Not sure why - when the same code was evaluated in a test it seemed almost instant.
     */
    LOG(1, "Setting up cells...");
    for (typename Tissue<DIM>::Iterator cell_iter = mrTissue.Begin();
         cell_iter != mrTissue.End();
         ++cell_iter)
    {
        /* We don't use the result; this call is just to force the cells 
         * to age to current time running their cell cycle models to get there.
         */
        cell_iter->ReadyToDivide();
    }
    LOG(1, "\tdone\n");
     
    // Write initial conditions to file for the visualizer.
    if(DIM==2)
    {
        WriteVisualizerSetupFile();
    }
    mpSetupFile->close();
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
                                                                mSimulationOutputDirectory+"/vis_results/",
                                                                "results.visvoronoi",
                                                                mFollowLoggedCell);
    }

    CancerEventHandler::EndEvent(SETUP);
                               
    /////////////////////////////////////////////////////////////////////
    // Main time loop
    /////////////////////////////////////////////////////////////////////
    while (p_simulation_time->GetTimeStepsElapsed() < num_time_steps)
    {        
        LOG(1, "--TIME = " << p_simulation_time->GetDimensionalisedTime() << "\n");
        
        // remove dead cells before doing birth
        // neither of these functions use any element information so they 
        // just delete and create nodes
        CancerEventHandler::BeginEvent(DEATH);
        mNumDeaths += DoCellRemoval();
        LOG(1, "\tNum deaths = " << mNumDeaths << "\n");
        CancerEventHandler::EndEvent(DEATH);


        CancerEventHandler::BeginEvent(BIRTH);
        mNumBirths += DoCellBirth();
        LOG(1, "\tNum births = " << mNumBirths << "\n");
        CancerEventHandler::EndEvent(BIRTH);
        

        if( (mNumBirths>0) || (mNumDeaths>0) )
        {   
            // If any nodes have been deleted or added we MUST call a ReMesh
            assert(mReMesh);
        }

        CancerEventHandler::BeginEvent(REMESH);
        if(mReMesh)
        {
            LOG(1, "\tRemeshing...");
            mrTissue.ReMesh();
            LOG(1, "\tdone.\n");
        }
        CancerEventHandler::EndEvent(REMESH);


        CancerEventHandler::BeginEvent(TESSELLATION);
        if(mWriteVoronoiData || mpMechanicsSystem->NeedsVoronoiTessellation())
        {
            mrTissue.CreateVoronoiTessellation();
        }
        CancerEventHandler::EndEvent(TESSELLATION);

        //  calculate node velocities
        CancerEventHandler::BeginEvent(VELOCITY);
        std::vector<c_vector<double, DIM> >& drdt = mpMechanicsSystem->rCalculateVelocitiesOfEachNode();
        CancerEventHandler::EndEvent(VELOCITY);


        // update node positions
        CancerEventHandler::BeginEvent(POSITION);
        UpdateNodePositions(drdt);
        CancerEventHandler::EndEvent(POSITION);
     
        // Increment simulation time here, so results files look sensible
        p_simulation_time->IncrementTimeOneStep();
        
        CancerEventHandler::BeginEvent(OUTPUT);
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
        CancerEventHandler::EndEvent(OUTPUT);
    }

    AfterSolve();
    
    // Write end state to tabulated files (not visualizer - this
    // is taken care of in the main loop).
    // Doesn't need to count cell types again as it is done in the last loop
    CancerEventHandler::BeginEvent(OUTPUT);
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
    CancerEventHandler::EndEvent(OUTPUT);
    CancerEventHandler::EndEvent(CANCER_EVERYTHING);
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
    if(mReMesh)
    {
        mrTissue.ReMesh();
    }
    
    // the false is so the directory isn't cleaned
    TrianglesMeshWriter<DIM,DIM> mesh_writer(archive_directory, mesh_filename, false);
    mesh_writer.WriteFilesUsingMesh(mrTissue.rGetMesh());
    
    // Create a new archive
    std::ofstream ofs(archive_filename.c_str());
    boost::archive::text_oarchive output_arch(ofs);
    
    // Save the simulation.  We save the time directly first to maintain its
    // singleton-ness on load.
    output_arch << *p_sim_time;
    
    // Archive the Wnt gradient if it's used
    bool archive_wnt = WntGradient::Instance()->IsGradientSetUp();
    output_arch & archive_wnt;
    if (archive_wnt)
    {
        WntGradient* p_wnt_gradient = WntGradient::Instance();
        output_arch & *p_wnt_gradient;
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
    Tissue<DIM>::meshPathname = mesh_filename;
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
    
    // Load Wnt gradient if it's used
    bool archive_wnt;
    rInputArch & archive_wnt;
    if (archive_wnt)
    {
        WntGradient* p_wnt_gradient = WntGradient::Instance();
        rInputArch & *p_wnt_gradient;
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
