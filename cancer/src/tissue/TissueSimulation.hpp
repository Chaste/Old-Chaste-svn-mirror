/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef TISSUESIMULATION_HPP_
#define TISSUESIMULATION_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include <vector>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <set>

#include "MeshBasedTissueWithGhostNodes.hpp"
#include "AbstractForce.hpp"
#include "CellwiseData.hpp"
#include "WntConcentration.hpp"
#include "AbstractCellKiller.hpp"
#include "RandomNumberGenerator.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "CancerEventHandler.hpp"
#include "LogFile.hpp"


/**
 * Run a 2D or 3D tissue simulation, currently based on the Meineke Paper
 * (doi:10.1046/j.0960-7722.2001.00216.x)
 *
 * Cells are represented by their centres in space, they are connected by
 * springs defined by the cells' Delaunay/Voronoi tesselation.
 *
 * The spring lengths are governed by the equations
 * dr/dt = stem_cycle_time*(mu/eta) sum_j r_hat_i,j*(|r_i,j|-s0)
 *       = alpha sum_j r_hat_i,j*(|r_i,j|-s0)
 *
 * where alpha = stem_cycle_time*(mu/eta) = stem_cycle_time*meineke_lambda.
 *       s0    = natural length of the spring.

 * Length is scaled by natural length.
 * Time is in hours.
 *
 * meineke_lambda = mu (spring constant) / eta (damping) = 0.01 (from Meineke - note
 * that the value we use for Meineke lambda is completely different because we have
 * nondimensionalised)
 *
 * The TissueSimulation accepts a tissue (facade class), containing either a mesh or
 * just a vector of nodes, where nodes are associated with TissueCells or are ghost
 * nodes. The TissueSimulation then accesses only the TissueCells via an iterator in
 * the tissue facade class.
 *
 * If simulating a crypt with a mesh-based tissue, the mesh should be surrounded by at
 * least one layer of ghost nodes. These are nodes which do not correspond to a cell,
 * but are necessary for remeshing (because the remesher tries to create a convex hull
 * of the set of nodes) and visualization purposes. The MeshBasedTissueWithGhostNodes
 * class deals with ghost nodes. SetGhostNodes() should have been called on it.
 *
 * Cells can divide (at a time governed by their cell cycle models).
 *
 * Cells can die - at a time/position specified by cell killers which can be
 * added to the simulation.
 *
 */
/// \todo Some methods in this class need documenting (see #736)
template<unsigned DIM>
class TissueSimulation
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestCryptSimulation2d;
    friend class TestTissueSimulation3d;

protected:

    /** TimeStep */
    double mDt;

    /** Time to run the Solve() method up to */
    double mEndTime;

    /** Facade encapsulating cells in the tissue being simulated */
    AbstractTissue<DIM>& mrTissue;
    
    /** Whether to delete the facade in our destructor */ 
    bool mDeleteTissue;

    bool mAllocatedMemoryForForceCollection;
    
    /** Whether to initialise the cells */
    bool mInitialiseCells;

    /** Whether to run the simulation with no birth (defaults to false). */
    bool mNoBirth;

    /** Whether to remesh at each timestep or not (defaults to true).*/
    bool mReMesh;

    /** Whether to count the number of each cell mutation state and output to file*/
    bool mOutputCellMutationStates;

    /** Whether to output the ancestor of each cell to a visualizer file*/
    bool mOutputCellAncestors;

    /** Whether to count the number of each cell type and output to file*/
    bool mOutputCellTypes;

    /** Whether to write the cell variables to a file */
    bool mOutputCellVariables;

    /** Whether to write the cell cycle phases to a file */
    bool mOutputCellCyclePhases;

    /** Output directory (a subfolder of tmp/<USERNAME>/testoutput) */
    std::string mOutputDirectory;

    /** Simulation Output directory either the same as mOutputDirectory or includes mOutputDirectory/results_from_time_<TIME> */
    std::string mSimulationOutputDirectory;

    /** Visualiser setup file */
    out_stream mpSetupFile;

    /** The Meineke and cancer parameters */
    CancerParameters *mpParams;

    /** The singleton RandomNumberGenerator */
    RandomNumberGenerator *mpRandomGenerator;

    /** Counts the number of births during the simulation */
    unsigned mNumBirths;

    /** Counts the number of deaths during the simulation */
    unsigned mNumDeaths;

    /**
     * The ratio of the number of actual timesteps to the number
     * of timesteps at which results are written to file
     * */
    unsigned mSamplingTimestepMultiple;

    /** List of cell killers */
    std::vector<AbstractCellKiller<DIM>*> mCellKillers;

    /** The mechanics used to determine the new location of the cells */
    std::vector<AbstractForce<DIM>*> mForceCollection;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        mpParams = CancerParameters::Instance();
        archive & *mpParams;
        archive & mpParams;

        mpRandomGenerator = RandomNumberGenerator::Instance();
        archive & *mpRandomGenerator;
        archive & mpRandomGenerator;

        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & mDt;
        archive & mEndTime;
        archive & mNoBirth;
        archive & mReMesh;
        archive & mOutputDirectory;
        archive & mNumBirths;
        archive & mNumDeaths;
        archive & mCellKillers;
        archive & mOutputCellMutationStates;
        archive & mOutputCellAncestors;
        archive & mOutputCellTypes;
        archive & mOutputCellVariables;
        archive & mOutputCellCyclePhases;
        archive & mSamplingTimestepMultiple;
        archive & mForceCollection;
    }

    /**
     * Writes out special information about the mesh to the visualizer.
     */
    virtual void WriteVisualizerSetupFile()
    {
    }

    /**
     * During a simulation time step, process any cell divisions that need to occur.
     * If the simulation includes cell birth, causes (almost) all cells that are ready to divide
     * to produce daughter cells.
     *
     * @return the number of births that occurred.
     */
    unsigned DoCellBirth();

    /**
     * Calculates the new locations of a dividing cell's cell centres.
     * Moves the dividing node a bit and returns co-ordinates for the new node.
     * It does this by picking a random direction (0->2PI) and placing the parent
     * and daughter in opposing directions on this axis.
     *
     * @param node_index The parent node index
     *
     * @return daughter_coords The coordinates for the daughter cell.
     *
     */
    virtual c_vector<double, DIM> CalculateDividingCellCentreLocations(typename AbstractTissue<DIM>::Iterator parentCell);

    /**
     * During a simulation time step, process any cell sloughing or death
     *
     * This uses the cell killers to remove cells and associated nodes from the
     * facade class.
     *
     * @return the number of deaths that occurred.
     */
    unsigned DoCellRemoval();

    /**
     * Moves each node to a new position for this timestep
     *
     * @param nodeForces the x and y force components on each node.
     */
    virtual void UpdateNodePositions(const std::vector< c_vector<double, DIM> >& nodeForces);

    /**
     *  A method for subclasses to do something at the end of each timestep
     */
    virtual void PostSolve()
    {
    }

    /**
     *  A method for subclasses to do something at before the start of the time loop
     */
    virtual void SetupSolve()
    {
    }

    /**
     *  Implements out cell birth, cell death and a remesh if necessary for a
     *  final time. This method may be overridden in subclasses to do something
     *  at the end of each time loop. Note that each subclass should also call
     *  the base class method.
     */
    virtual void AfterSolve();
    
    /**
     *  A child class can overload this if they want the simulation to stop 
     *  based on certain conditions before the specified end time (for example,
     *  run until a crypt becomes monoclonal).
     */
    virtual bool StoppingEventHasOccurred();

public:

    /**
     *  Constructor
     *
     *  @param rTissue A tissue facade class (contains a mesh and cells)
     *  @param forceCollection The mechanics to use in the simulation
     *  @param deleteTissueAndForceCollection Whether to delete the tissue and force collection on destruction to free up memory
     *  @param initialiseCells Whether to initialise cells (set to false when loading from an archive)
     */
    TissueSimulation(AbstractTissue<DIM>& rTissue,
                     std::vector<AbstractForce<DIM>*> forceCollection,
                     bool deleteTissueAndForceCollection=false,
                     bool initialiseCells=true);

    /**
     * Destructor.
     */
    virtual ~TissueSimulation();

    /**
     * Get methods.
     */
    std::vector<double> GetNodeLocation(const unsigned& rNodeIndex);
    c_vector<unsigned, NUM_CELL_MUTATION_STATES> GetCellMutationStateCount();
    c_vector<unsigned, NUM_CELL_TYPES> GetCellTypeCount();
    c_vector<unsigned, 5> GetCellCyclePhaseCount();

    double GetDt();
    unsigned GetNumBirths();
    unsigned GetNumDeaths();
    
    std::string GetOutputDirectory();
    
    /**
     * Set methods.
     */
    void SetDt(double dt);
    void SetEndTime(double endTime);
    void SetOutputDirectory(std::string outputDirectory);
    void SetSamplingTimestepMultiple(unsigned samplingTimestepMultiple);
    void SetNoBirth(bool nobirth);
    void SetOutputCellMutationStates(bool outputCellMutationStates);
    void SetOutputCellAncestors(bool outputCellAncestors);
    void SetOutputCellTypes(bool outputCellTypes);
    void SetOutputCellVariables(bool outputCellVariables);
    void SetOutputCellCyclePhases(bool outputCellCyclePhases);
    void SetReMeshRule(bool remesh);

    void AddCellKiller(AbstractCellKiller<DIM>* pCellKiller);

    void Solve();

    AbstractTissue<DIM>& rGetTissue();
    const AbstractTissue<DIM>& rGetTissue() const;
        
    const std::vector<AbstractForce<DIM>*> rGetForceCollection() const;
    
    /**
     * Apply any tissue boundary conditions. Can be overridden in subclasses.
     */
    virtual void ApplyTissueBoundaryConditions(TissueCell& rCell, ChastePoint<DIM>& rPoint)
    {}

};


template<unsigned DIM>
TissueSimulation<DIM>::TissueSimulation(AbstractTissue<DIM>& rTissue,
                                        std::vector<AbstractForce<DIM>*> forceCollection,
                                        bool deleteTissueAndForceCollection,
                                        bool initialiseCells)
  :  mrTissue(rTissue)
{
    #define COVERAGE_IGNORE
    assert(DIM==2 || DIM==3); // there are no instances of TissueSimulation<1>
    #undef COVERAGE_IGNORE
    
    mDeleteTissue = deleteTissueAndForceCollection;
    
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
    
    mAllocatedMemoryForForceCollection = deleteTissueAndForceCollection;
    
    mForceCollection = forceCollection;
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
    if (mAllocatedMemoryForForceCollection) 
    {    
        for (typename std::vector<AbstractForce<DIM>*>::iterator force_iter = mForceCollection.begin();
             force_iter != mForceCollection.end();
             ++force_iter)
        {     
            delete *force_iter;
        }
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
    for (typename std::vector<AbstractCellKiller<DIM>*>::iterator killer_iter = mCellKillers.begin();
                 killer_iter != mCellKillers.end();
                 ++killer_iter)
    {
        (*killer_iter)->TestAndLabelCellsForApoptosisOrDeath();
    }

    num_deaths_this_step += mrTissue.RemoveDeadCells();

    return num_deaths_this_step;
}

template<unsigned DIM>
const std::vector<AbstractForce<DIM>*> TissueSimulation<DIM>::rGetForceCollection() const
{
    return mForceCollection;
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
void TissueSimulation<DIM>::UpdateNodePositions(const std::vector< c_vector<double, DIM> >& nodeForces)
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
        unsigned index = cell.GetLocationIndex();
        double damping_const = mrTissue.GetDampingConstant(cell);
        
        ChastePoint<DIM> new_point(mrTissue.GetNode(index)->rGetLocation() + mDt*nodeForces[index]/damping_const);
        
        ApplyTissueBoundaryConditions(cell, new_point);
        
        // Move the cell
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
 * Get the number of births that have occurred in the entire simulation (since t=0)
 */
template<unsigned DIM>
unsigned TissueSimulation<DIM>::GetNumBirths()
{
    return mNumBirths;
}

/**
 * Get the number of deaths that have occurred in the entire simulation (since t=0).
 */
template<unsigned DIM>
unsigned TissueSimulation<DIM>::GetNumDeaths()
{
    return mNumDeaths;
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
 * Get the output directory of the simulation.
 */
template<unsigned DIM>
std::string TissueSimulation<DIM>::GetOutputDirectory()
{
    return mOutputDirectory;
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
 *
 * test for this is in TestCryptSimulation2d::TestStandardResultForArchivingTestsBelow().
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
    for (unsigned i=0; i<DIM; i++)
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
    CancerEventHandler::BeginEvent(CANCER_EVERYTHING);
    CancerEventHandler::BeginEvent(SETUP);

    // Set up the simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();

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

    double time_now = p_simulation_time->GetTime();
    std::ostringstream time_string;
    time_string << time_now;

    std::string results_directory = mOutputDirectory +"/results_from_time_" + time_string.str();
    mSimulationOutputDirectory = results_directory;

    ///////////////////////////////////////////////////////////
    // Set up Simulation
    ///////////////////////////////////////////////////////////

    // Create output files for the visualizer
    OutputFileHandler output_file_handler(results_directory+"/", true);

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
    *mpSetupFile << std::flush;

    mrTissue.WriteResultsToFiles(mOutputCellMutationStates,
                                 mOutputCellTypes,
                                 mOutputCellVariables,
                                 mOutputCellCyclePhases,
                                 mOutputCellAncestors);

    CancerEventHandler::EndEvent(SETUP);

    /////////////////////////////////////////////////////////////////////
    // Main time loop
    /////////////////////////////////////////////////////////////////////
    
    while ((p_simulation_time->GetTimeStepsElapsed() < num_time_steps) && !(StoppingEventHasOccurred()) )
    {
        LOG(1, "--TIME = " << p_simulation_time->GetTime() << "\n");

        // Remove dead cells then implement cell birth. Note that neither
        // of these methods use any element information, they just delete
        // and create nodes
        CancerEventHandler::BeginEvent(DEATH);
        unsigned deaths_this_step = DoCellRemoval();
        mNumDeaths += deaths_this_step;
        LOG(1, "\tNum deaths = " << mNumDeaths << "\n");
        CancerEventHandler::EndEvent(DEATH);

        CancerEventHandler::BeginEvent(BIRTH);
        unsigned births_this_step = DoCellBirth();
        mNumBirths += births_this_step;
        LOG(1, "\tNum births = " << mNumBirths << "\n");
        CancerEventHandler::EndEvent(BIRTH);

        // If the tissue has a mesh, then we currently must call a ReMesh at
        // each timestep. Otherwise, we only need to call a ReMesh after there
        // has been any cell birth or cell death.
        if (mrTissue.HasMesh())
        {
            //This assertion is not necessarily true
            //assert(mReMesh);
            //See TestCryptSimulation2dNightly::Test2DSpringSystem where the
            //default value is over-written
        }
        else
        {
            mReMesh = false;
            if ( (births_this_step>0) || (deaths_this_step>0) )
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
                 || (static_cast<MeshBasedTissue<DIM>*>(&mrTissue))->UseAreaBasedDampingConstant()
                 || (static_cast<MeshBasedTissue<DIM>*>(&mrTissue))->GetWriteTissueAreas() )
            {
                (static_cast<MeshBasedTissue<DIM>*>(&mrTissue))->CreateVoronoiTessellation();
            }
        }
        CancerEventHandler::EndEvent(TESSELLATION);

        // Calculate Forces
        CancerEventHandler::BeginEvent(FORCE);
        
        // Initialise a vector of forces on node
        std::vector<c_vector<double, DIM> > forces;
        forces.reserve(mrTissue.GetNumNodes());
        
        for (unsigned i=0; i<mrTissue.GetNumNodes(); i++)
        {
             forces.push_back(zero_vector<double>(DIM));
        }
        
        for (typename std::vector<AbstractForce<DIM>*>::iterator iter = mForceCollection.begin();
             iter !=mForceCollection.end();
             iter++)
        {
            (*iter)->AddForceContribution(forces, mrTissue);
        }
        CancerEventHandler::EndEvent(FORCE);

        // Update node positions
        CancerEventHandler::BeginEvent(POSITION);
        UpdateNodePositions(forces);
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

    *mpSetupFile << "\nComplete\n";
    mpSetupFile->close();

    CancerEventHandler::EndEvent(OUTPUT);

    CancerEventHandler::EndEvent(CANCER_EVERYTHING);
}

template<unsigned DIM>
void TissueSimulation<DIM>::AfterSolve()
{
    LOG(1, "--TIME = " << SimulationTime::Instance()->GetTime() << "\n");

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
        //This assertion is not necessarily true
        //assert(mReMesh);
        //See TestCryptSimulation2dNightly::Test2DSpringSystem where the
        //default value is over-written
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


template<unsigned DIM>
bool TissueSimulation<DIM>::StoppingEventHasOccurred()
{
    return false;
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
    if (!mOutputCellMutationStates)
    {
        EXCEPTION("Call simulator.SetOutputCellMutationStates before using this function");
    }
    return mrTissue.GetCellMutationStateCount();
}

/**
 * Find out how many cells of each type there are
 *
 * @return The number of cells of each type (evaluated at each visualizer output)
 * [0] = STEM
 * [1] = TRANSIT
 * [2] = DIFFERENTIATED
 * [3] = APOPTOTIC
 */
template<unsigned DIM>
c_vector<unsigned, NUM_CELL_TYPES> TissueSimulation<DIM>::GetCellTypeCount()
{
    if (!mOutputCellTypes)
    {
        EXCEPTION("Call simulator.SetOutputCellTypes() before using this function");
    }
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
    if (!mOutputCellCyclePhases)
    {
        EXCEPTION("Call simulator.SetOutputCellCyclePhases() before using this function");
    }
    return mrTissue.GetCellCyclePhaseCount();
}

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSimulation.
 *
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const TissueSimulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<DIM> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
    const std::vector<AbstractForce<DIM>*> force_collection = t->rGetForceCollection();
    ar & force_collection;
}

/**
 * De-serialize constructor parameters and initialise Tissue.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, TissueSimulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<DIM>* p_tissue;
    ar >> p_tissue;
    std::vector<AbstractForce<DIM>*> force_collection;
    ar >> force_collection;

    // Invoke inplace constructor to initialize instance
    ::new(t)TissueSimulation<DIM>(*p_tissue, force_collection, true);
}
}
} // namespace ...


#endif /*TISSUESIMULATION_HPP_*/
