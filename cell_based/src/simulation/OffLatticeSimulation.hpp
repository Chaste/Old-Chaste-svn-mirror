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

#ifndef OFFLATTICESIMULATION_HPP_
#define OFFLATTICESIMULATION_HPP_

#include "AbstractCellBasedSimulation.hpp"
#include "AbstractForce.hpp"
#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

/**
 * Run an off-lattice 2D or 3D cell-based simulation using a cell-centre-
 * or vertex-based cell population.
 *
 * In cell-centre-based cell populations, each cell is represented by a
 * single node (corresponding to its centre), and connectivity is defined
 * either by a Delaunay triangulation or a radius of influence. In vertex-
 * based cell populations, each cell is represented by a polytope
 * (corresponding to its membrane) with a variable number of vertices.
 *
 * The OffLatticeSimulation is constructed with a CellPopulation, which
 * updates the correspondence between each Cell and its spatial representation
 * and handles cell division (governed by the CellCycleModel associated
 * with each cell). Once constructed, one or more Force laws may be passed
 * to the OffLatticeSimulation object, to define the mechanical properties
 * of the CellPopulation. Similarly, one or more CellKillers may be passed
 * to the OffLatticeSimulation object to specify conditions in which Cells
 * may die, and one or more CellPopulationBoundaryConditions to specify
 * regions in space beyond which Cells may not move.
 */
template<unsigned DIM>
class OffLatticeSimulation : public AbstractCellBasedSimulation<DIM>
{

protected:
    /** The mechanics used to determine the new location of the cells, a list of the forces. */
    std::vector<AbstractForce<DIM>*> mForceCollection;

    /** List of boundary conditions. */
    std::vector<AbstractCellPopulationBoundaryCondition<DIM>*> mBoundaryConditions;

    /** Whether to write the node velocities to a file. */
    bool mOutputNodeVelocities;

    /** Results file node velocities. */
    out_stream mpNodeVelocitiesFile;

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
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<AbstractCellBasedSimulation<DIM> >(*this);
        archive & mForceCollection;
        archive & mBoundaryConditions;
        archive & mOutputNodeVelocities;
    }

    /**
     * Overridden this UpdateCellLocationsAndTopology method.
     *
     * In Off lattice simulations this method calculates forces and updates node positions.
     */
    virtual void UpdateCellLocationsAndTopology();

    /**
     * Moves each node to a new position for this timestep by
     * calling the CellPopulation::UpdateNodeLocations() method then
     * applying any boundary conditions.
     *
     * @param rNodeForces the forces on nodes
     */
    virtual void UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rNodeForces);

    /**
     * Overridden SetupSolve() method to setup the node velocities file.
     */
    virtual void SetupSolve();

    /**
     * Overridden AfterSolve() method to close the node velocities file.
     */
    virtual void AfterSolve();


public:
    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationAndCellKillersInDestructor Whether to delete the cell population cell killer, force, and boundary conditions collections on destruction to free up memory
     * @param initialiseCells Whether to initialise cells (set to false when loading from an archive)
     */
    OffLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                         bool deleteCellPopulationAndCellKillersInDestructor=false,
                         bool initialiseCells=true);

    /**
     * Destructor.
     *
     * This frees the Forces and Boundary Conditions, if they were created by de-serialization.
     */
    virtual ~OffLatticeSimulation();

    /**
     * Add a force to be used in this simulation (use this to set the mechanics system).
     *
     * @param pForce pointer to a force law
     */
    void AddForce(AbstractForce<DIM>* pForce);

    /**
     * Add a cell population boundary condition to be used in this simulation.
     *
     * @param pBoundaryCondition pointer to a boundary condition
     */
    void AddCellPopulationBoundaryCondition(AbstractCellPopulationBoundaryCondition<DIM>* pBoundaryCondition);

    /**
     * @return mOutputNodeVelocities
     */
    bool GetOutputNodeVelocities();

    /**
     * Set mOutputNodeVelocities.
     *
     * @param outputNodeVelocities the new value of mOutputNodeVelocities
     */
    void SetOutputNodeVelocities(bool outputNodeVelocities);

    /**
     * Overridden OutputAdditionalSimulationSetup method to output the force and cell
     * population boundary condition information
     */
    void OutputAdditionalSimulationSetup(out_stream& rParamsFile);

    /**
     * Outputs simulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationParameters(out_stream& rParamsFile);
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OffLatticeSimulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct an OffLatticeSimulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const OffLatticeSimulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise an OffLatticeSimulation.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, OffLatticeSimulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance, last two variables set extra
    // member variables to be deleted as they are loaded from archive and to not initialise sells.
    ::new(t)OffLatticeSimulation<DIM>(*p_cell_population, true, false);
}
}
} // namespace

#endif /*OFFLATTICESIMULATION_HPP_*/
