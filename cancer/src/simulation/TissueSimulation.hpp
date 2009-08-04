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
 * Run an off-lattice 2D or 3D cell-based simulation using a cell-centre- or vertex-based
 * tissue.
 *
 * In cell-centre-based tissues, each cell is represented by a single node (corresponding
 * to its centre), and connectivity is defined either by a Delaunay triangulation or
 * a radius of influence. In vertex-based tissues, each cell is represented by a polytope
 * (corresponding to its membrane) with a variable number of vertices.
 *
 * The TissueSimulation is constructed with a Tissue, which updates the correspondence
 * between each TissueCell and its spatial representation and handles cell division (governed
 * by the CellCycleModel associated with each cell); and one or more Force laws, which define
 * the mechanical properties of the Tissue. It is also possible to add one or more CellKiller
 * objects to the TissueSimulation, which specify the conditions under which a TissueCell dies.
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

    /** Whether to delete the tissue facade in the destructor. */
    bool mDeleteTissue;

    /** Whether delete the collection of force laws in the destructor. */
    bool mAllocatedMemoryForForceCollection;

    /** Whether to initialise the cells. */
    bool mInitialiseCells;

    /** Whether to run the simulation with no birth (defaults to false). */
    bool mNoBirth;

    /** Whether to update the topology of the tissue at each time step (defaults to true).*/
    bool mUpdateTissue;

    /** Output directory (a subfolder of tmp/[USERNAME]/testoutput) */
    std::string mOutputDirectory;

    /** Simulation Output directory either the same as mOutputDirectory or includes mOutputDirectory/results_from_time_[TIME] */
    std::string mSimulationOutputDirectory;

    /** Visualiser setup file */
    out_stream mpSetupFile;

    /** The cancer tissue configuration */
    TissueConfig *mpConfig;

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
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        mpConfig = TissueConfig::Instance();
        archive & *mpConfig;
        archive & mpConfig;

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
     * Method for determining how cell division occurs. This method returns a vector
     * which is then passed into the Tissue method AddCell(). This method may be
     * overridden by subclasses.
     * 
     * For a cell-centre tissue, this method calculates the new locations of the cell
     * centres of a dividing cell, moves the parent cell and returns the location of
     * the daughter cell. The new locations are found by picking a random direction
     * and placing the parent and daughter in opposing directions along this axis.
     * 
     * For a vertex tissue, the method returns the zero vector.
     *
     * @param pParentCell pointer to the parent cell
     *
     * @return a vector containing information on cell division.
     *
     */
    virtual c_vector<double, DIM> CalculateCellDivisionVector(TissueCell* pParentCell);

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
     * @param rOldLocations the node locations before any boundary conditions are applied
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
     * Get a node's location (ONLY FOR TESTING).
     *
     * @param rNodeIndex the node index
     * @return the co-ordinates of this node.
     */
    std::vector<double> GetNodeLocation(const unsigned& rNodeIndex);

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
    const AbstractTissue<DIM> *p_tissue = &(t->rGetTissue());
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
    AbstractTissue<DIM> *p_tissue;
    ar >> p_tissue;
    std::vector<AbstractForce<DIM>*> force_collection;
    ar >> force_collection;

    // Invoke inplace constructor to initialise instance
    ::new(t)TissueSimulation<DIM>(*p_tissue, force_collection, true);
}
}
} // namespace


#endif /*TISSUESIMULATION_HPP_*/
