/*

Copyright (C) University of Oxford, 2005-2009

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

#include <climits> // work around boost bug
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include <vector>

#include "AbstractForce.hpp"
#include "AbstractCellKiller.hpp"
#include "AbstractTissue.hpp"
#include "RandomNumberGenerator.hpp"
#include "ChastePoint.hpp"


/**
 * Run a 2D or 3D tissue simulation, currently based on the Meineke Paper
 * (doi:10.1046/j.0960-7722.2001.00216.x)
 *
 * Cells are represented by their centres in space, they are connected by
 * springs defined by the cells' Delaunay/Voronoi tessellation.
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
 * class deals with ghost nodes.
 *
 * Cells can divide (at a time governed by their cell cycle models).
 *
 * Cells can die - at a time/position specified by cell killers which can be
 * added to the simulation.
 *
 */
template<unsigned DIM>
class TissueSimulation
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestCryptSimulation2d;
    friend class TestTissueSimulation3d;

protected:

    /** Time step. */
    double mDt;

    /** Time to run the Solve() method up to. */
    double mEndTime;

    /** Facade encapsulating cells in the tissue being simulated. */
    AbstractTissue<DIM>& mrTissue;

    /** Whether to delete the facade in the destructor. */
    bool mDeleteTissue;

    /** Whether delete the collection of force laws in the destructor. */
    bool mAllocatedMemoryForForceCollection;

    /** Whether to initialise the cells. */
    bool mInitialiseCells;

    /** Whether to run the simulation with no birth (defaults to false). */
    bool mNoBirth;

    /** Whether to update the topology of the tissue at each time step (defaults to true).*/
    bool mUpdateTissue;

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

    /** Whether to write the cell ages to a file */
    bool mOutputCellAges;

    /** Output directory (a subfolder of tmp/[USERNAME]/testoutput) */
    std::string mOutputDirectory;

    /** Simulation Output directory either the same as mOutputDirectory or includes mOutputDirectory/results_from_time_[TIME] */
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

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * Serialization of singleton objects must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     *
     * @param archive
     * @param version
     */
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
        archive & mUpdateTissue;
        archive & mOutputDirectory;
        archive & mNumBirths;
        archive & mNumDeaths;
        archive & mCellKillers;
        archive & mOutputCellMutationStates;
        archive & mOutputCellAncestors;
        archive & mOutputCellTypes;
        archive & mOutputCellVariables;
        archive & mOutputCellCyclePhases;
        archive & mOutputCellAges;
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
     * Calculate the new locations of the cell centres of a dividing cell, move 
     * the parent cell and return the location of the daughter cell.
     * 
     * The new locations are found by picking a random direction and placing the 
     * parent and daughter in opposing directions along this axis.
     *
     * @param pParentCell pointer to the parent cell
     *
     * @return daughter_coords The coordinates for the daughter cell.
     *
     */
    virtual c_vector<double, DIM> CalculateDividingCellCentreLocations(TissueCell* pParentCell);

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
     * Moves each node to a new position for this timestep by
     * calling the tissue UpdateNodeLocations() method then
     * applying any boundary conditions.
     *
     * @param rNodeForces the forces on nodes
     */
    virtual void UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rNodeForces);

    /**
     * Apply any tissue boundary conditions. Can be overridden in subclasses.
     *
     * @param rOldLocations
     */
    virtual void ApplyTissueBoundaryConditions(const std::vector< c_vector<double, DIM> >& rOldLocations)
    {
    }

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
     *  This method may be overridden in subclasses to do something
     *  at the end of each time loop.
     */
    virtual void AfterSolve(){};

    /**
     *  A child class can overload this if they want the simulation to stop
     *  based on certain conditions before the specified end time (for example,
     *  run until a crypt becomes monoclonal).
     */
    virtual bool StoppingEventHasOccurred();

    /**
     * Calls the deaths, births and (if mUpdateTissue is true) Tissue::Update() methods.
     */
    void UpdateTissue();

public:

    /**
     *  Constructor.
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
     *
     * This frees the tissue and cell killers, if they were created by de-serialization.
     */
    virtual ~TissueSimulation();

    /**
     * Get a node's location (ONLY FOR TESTING)
     *
     * @param rNodeIndex  the node index
     * @return the co-ordinates of this node.
     */
    std::vector<double> GetNodeLocation(const unsigned& rNodeIndex);

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
    c_vector<unsigned, NUM_CELL_MUTATION_STATES> GetCellMutationStateCount();

    /**
     * Find out how many cells of each type there are
     *
     * @return The number of cells of each type (evaluated at each visualizer output)
     * [0] = STEM
     * [1] = TRANSIT
     * [2] = DIFFERENTIATED
     * [3] = APOPTOTIC
     */
    c_vector<unsigned, NUM_CELL_TYPES> GetCellTypeCount();

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
    c_vector<unsigned, 5> GetCellCyclePhaseCount();

    /**
     * @return the timestep of the simulation
     */
    double GetDt();

    /**
     * @return the number of births that have occurred in the entire simulation (since t=0)
     */
    unsigned GetNumBirths();

    /**
     * @return the number of deaths that have occurred in the entire simulation (since t=0).
     */
    unsigned GetNumDeaths();

    /**
     * Get the output directory of the simulation.
     */
    std::string GetOutputDirectory();

    /**
     * Set the timestep of the simulation.
     *
     * @param dt the timestep to use
     */
    void SetDt(double dt);

    /**
     * Set the end time and resets the timestep to be endtime/100.
     *
     * @param endTime the end time to use
     */
    void SetEndTime(double endTime);

    /**
     * Set the output directory of the simulation.
     *
     * @param outputDirectory the output directory to use
     */
    void SetOutputDirectory(std::string outputDirectory);

    /**
     * Set the ratio of the number of actual timesteps to the number of timesteps
     * at which results are written to file. Default value is set to 1 by the constructor.
     *
     * @param samplingTimestepMultiple the ratio to use
     */
    void SetSamplingTimestepMultiple(unsigned samplingTimestepMultiple);

    /**
     * Set the simulation to run with no birth.
     *
     * @param noBirth whether to run with no birth
     */
    void SetNoBirth(bool noBirth);

    /**
     * Set the simulation to count and store the number of each cell mutation state.
     *
     * @param outputCellMutationStates whether to output cell mutation states
     */
    void SetOutputCellMutationStates(bool outputCellMutationStates);

    /**
     * Set the simulation to output the cell ancestors if they are set.
     *
     * @param outputCellAncestors whether to output cell ancestors
     */
    void SetOutputCellAncestors(bool outputCellAncestors);

    /**
     * Set the simulation to count and store the number of each cell type.
     *
     * @param outputCellTypes whether to output cell types
     */
    void SetOutputCellTypes(bool outputCellTypes);

    /**
     * Set the simulation to output the cell-cycle variables.
     *
     * @param outputCellVariables whether to output cell-cycle variables
     */
    void SetOutputCellVariables(bool outputCellVariables);

    /**
     * Set the simulation to output the cell cycle phases.
     *
     * The test for this method is in TestCryptSimulation2d::TestStandardResultForArchivingTestsBelow().
     *
     * @param outputCellCyclePhases whether to output cell-cycle phases
     */
    void SetOutputCellCyclePhases(bool outputCellCyclePhases);

    /**
     * Set the simulation to output the cell ages.
     * 
     * @param outputCellAges whether to output cell ages
     */
    void SetOutputCellAges(bool outputCellAges);

    /**
     * Set whether to update the topology of the tissue at each time step.
     *
     * @param updateTissue  whether to update the tissue each time step
     */
    void SetUpdateTissueRule(bool updateTissue);

    /**
     * Add a cell killer to be used in this simulation.
     *
     * @param pCellKiller pointer to a cell killer
     */
    void AddCellKiller(AbstractCellKiller<DIM>* pCellKiller);

    /**
     * Main solve method.
     *
     * This method sets up the simulation time, creates output files, and initialises the
     * tissue. It then iterates through a time loop. At each time step, first any cell death
     * or birth is implemented, then the tissue topology is updated, then the forces are
     * recalculated and the tissue evolved according to whatever force laws are present in
     * the simulation, and finally the results for that time step are output to file. At the
     * end of the time loop, the method closes any output files.
     */
    void Solve();

    /**
     * @return reference to the tissue.
     */
    AbstractTissue<DIM>& rGetTissue();

    /**
     * @return const reference to the tissue (used in archiving).
     */
    const AbstractTissue<DIM>& rGetTissue() const;

    /**
     * @return const reference to mForceCollection (used in archiving).
     */
    const std::vector<AbstractForce<DIM>*> rGetForceCollection() const;
};


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
 * De-serialize constructor parameters and initialise a TissueSimulation.
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

    // Invoke inplace constructor to initialise instance
    ::new(t)TissueSimulation<DIM>(*p_tissue, force_collection, true);
}
}
} // namespace


#endif /*TISSUESIMULATION_HPP_*/
