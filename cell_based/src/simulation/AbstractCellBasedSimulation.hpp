/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef ABSTRACTCELLBASEDSIMULATION_HPP_
#define ABSTRACTCELLBASEDSIMULATION_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include <vector>

#include "AbstractForce.hpp"
#include "AbstractCellKiller.hpp"
#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "AbstractCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "Identifiable.hpp"

/**
 * An abstract cell based simulation class. This class contains common functionality
 * from off lattice and on lattice simulations.
 *
 * The AbstractCellBasedSimulation is constructed with a CellPopulation, which
 * updates the correspondence between each Cell and its spatial representation
 * and handles cell division (governed by the CellCycleModel associated
 * with each cell). Once constructed,  one or more CellKillers may be passed
 * to the AbstractCellBasedSimulation object to specify conditions in which Cells
 * may die,
 *
 * Subclasses use one or more Force laws or update rules (Which are passed
 * to the child class object) to define the mechanical properties of the CellPopulation.
 */
template<unsigned DIM>
class AbstractCellBasedSimulation : public Identifiable
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestCryptSimulation2d;
    friend class TestOffLatticeSimulation3d;
    friend class TestOffLatticeSimulation;

protected:

    /** Time step. */
    double mDt;

    /** Time to run the Solve() method up to. */
    double mEndTime;

    /** Facade encapsulating cells in the cell population being simulated. */
    AbstractCellPopulation<DIM>& mrCellPopulation;

    /** Whether to delete the cell population, force laws and boundary conditions in the destructor. */
    bool mDeleteCellPopulationAndForcesAndBCsInDestructor;

    /** Whether to initialise the cells. */
    bool mInitialiseCells;

    /** Whether to run the simulation with no birth (defaults to false). */
    bool mNoBirth;

    /** Whether to update the topology of the cell population at each time step (defaults to true).*/
    bool mUpdateCellPopulation;

    /** Output directory (a subfolder of tmp/[USERNAME]/testoutput). */
    std::string mOutputDirectory;

    /** Simulation Output directory either the same as mOutputDirectory or includes mOutputDirectory/results_from_time_[TIME]. */
    std::string mSimulationOutputDirectory;

    /** Visualiser setup file. */
    out_stream mpVizSetupFile;

    /** The singleton RandomNumberGenerator. */
    RandomNumberGenerator* mpRandomGenerator;

    /** Counts the number of births during the simulation. */
    unsigned mNumBirths;

    /** Counts the number of deaths during the simulation. */
    unsigned mNumDeaths;

    /**
     * The ratio of the number of actual timesteps to the number
     * of timesteps at which results are written to file.
     */
    unsigned mSamplingTimestepMultiple;

    /** List of cell killers. */
    std::vector<AbstractCellKiller<DIM>*> mCellKillers;

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
        mpRandomGenerator = RandomNumberGenerator::Instance();
        archive & *mpRandomGenerator;
        archive & mpRandomGenerator;

        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & mDt;
        archive & mEndTime;
        archive & mNoBirth;
        archive & mUpdateCellPopulation;
        archive & mOutputDirectory;
        archive & mNumBirths;
        archive & mNumDeaths;
        archive & mCellKillers;
        archive & mSamplingTimestepMultiple;
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
    virtual unsigned DoCellBirth();

    /**
     * Method for determining how cell division occurs. This method returns a vector
     * which is then passed into the CellPopulation method AddCell(). This method may be
     * overridden by subclasses.
     *
     * For a centre-based cell population, this method calculates the new locations of the cell
     * centres of a dividing cell, moves the parent cell and returns the location of
     * the daughter cell. The new locations are found by picking a random direction
     * and placing the parent and daughter in opposing directions along this axis.
     *
     * For a vertex-based cell population, the method returns the zero vector.
     *
     * @param pParentCell the parent cell
     *
     * @return a vector containing information on cell division.
     *
     */
    virtual c_vector<double, DIM> CalculateCellDivisionVector(CellPtr pParentCell);

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
     * A method for subclasses to do something at the end of each timestep
     */
    virtual void PostSolve()
    {
    }

    /**
     * A method for subclasses to do something at before the start of the time loop.
     */
    virtual void SetupSolve()
    {
    }

    /**
     * This method may be overridden in subclasses to do something
     * at the end of each time loop.
     */
    virtual void AfterSolve()
    {
    }

    /**
     * A child class can overload this if they want the simulation to stop
     * based on certain conditions before the specified end time (for example,
     * run until a crypt becomes monoclonal).
     */
    virtual bool StoppingEventHasOccurred();

    /**
     * Calls the deaths, births and (if mUpdateCellPopulation is true) CellPopulation::Update() methods.
     */
    void UpdateCellPopulation();

    /**
     * Updates the locations and topology of cells, must be defined in subclasses.
     *
     * For example this calculates forces and updates node positions in off lattice simulations.
     */
    virtual void UpdateCellLocationsAndTopology()=0;

    /**
     * Helper method to output all the simulations parameters and information to file.
     */
    void OutputSimulationSetup();

    /**
     * Helper method to output additional simulations parameters and information defined in
     * subclasses to file.
     */
    virtual void OutputAdditionalSimulationSetup(out_stream& rParamsFile)
    {
    }


public:

    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationAndForceCollection Whether to delete the cell population and cell killer collection on destruction to free up memory
     * @param initialiseCells Whether to initialise cells (set to false when loading from an archive)
     */
    AbstractCellBasedSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                        bool deleteCellPopulationAndForceCollection=false,
                        bool initialiseCells=true);

    /**
     * Destructor.
     *
     * This frees the cell population and cell killers, if they were created by de-serialization.
     */
    virtual ~AbstractCellBasedSimulation();

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
     * Set whether to update the topology of the cell population at each time step.
     *
     * @param updateCellPopulation  whether to update the cell population each time step
     */
    void SetUpdateCellPopulationRule(bool updateCellPopulation);

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
     * cell population. It then iterates through a time loop. At each time step, first any cell death
     * or birth is implemented, then the cell population topology is updated, then the forces are
     * recalculated and the cell population evolved according to whatever force laws are present in
     * the simulation, and finally the results for that time step are output to file. At the
     * end of the time loop, the method closes any output files.
     */
    void Solve();

    /**
     * @return reference to the cell population.
     */
    AbstractCellPopulation<DIM>& rGetCellPopulation();

    /**
     * @return const reference to the cell population (used in archiving).
     */
    const AbstractCellPopulation<DIM>& rGetCellPopulation() const;

    /**
     * Outputs simulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSimulationParameters(out_stream& rParamsFile)=0;
};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AbstractCellBasedSimulation)


namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a AbstractCellBasedSimulation.
 *
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const AbstractCellBasedSimulation<DIM>* t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a AbstractCellBasedSimulation.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, AbstractCellBasedSimulation<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)AbstractCellBasedSimulation<DIM>(*p_cell_population, true, false);
}
}
} // namespace


#endif /*ABSTRACTCELLBASEDSIMULATION_HPP_*/
