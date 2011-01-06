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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTCREATINGACELLBASEDSIMULATIONWITHBOUNDARYCONDITIONSTUTORIAL_HPP_
#define TESTCREATINGACELLBASEDSIMULATIONWITHBOUNDARYCONDITIONSTUTORIAL_HPP_

/*
 * = An example showing how to create a cell-based simulation with boundary conditions =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how to create a new cell-based simulation class in which cells
 * are constrained to lie within a fixed domain.
 *
 * EMPTYLINE
 *
 * == 1. Including header files ==
 *
 * EMPTYLINE
 *
 * The first thing to do is include the following header, which allows us
 * to use certain methods in our test (this header file should be included
 * in any Chaste test):
 */
#include <cxxtest/TestSuite.h>

/* Any test in which the {{{GetIdentifier()}}} method is used, 
 * even via the main cell_based code ({{{AbstraceCellPopulation}}} output methods), must 
 * include {{{CheckpointArchiveTypes.hpp}}} 
 * or {{{CellBasedSimulationArchiver.hpp}}} as the first Chaste header included. 
 */
#include "CheckpointArchiveTypes.hpp" 

/* The next header defines a base class for cell-based simulations, from which
 * the new class will inherit. */
#include "CellBasedSimulation.hpp"
/* The remaining header files define classes that will be used in the cell population
 * simulation test: {{{HoneycombMeshGenerator}}} defines a helper class for
 * generating a suitable mesh;  {{{CellsGenerator}}}
 * defines a helper class for generating a vector of cells and
 * {{{FixedDurationGenerationBasedCellCycleModel}}} makes them have fixed cell
 * cycle models; {{{GeneralisedLinearSpringForce}}} defines a force law for
 * describing the mechanical interactions between neighbouring cells in the
 * cell population; and {{{CellBasedSimulation}}} defines the class that simulates the
 * evolution of a cell population. */
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"

/*
 * EMPTYLINE
 *
 * == Defining the cell-based simulation class ==
 *
 * As an example, let us consider a two-dimensional cell-based simulation in which
 * all cells are constrained to lie within the domain given in Cartesian coordinates
 * by 0 <= y <= 5. To implement this we define a cell-based simulation class,
 * {{{MyCellBasedSimulation}}}, which inherits from {{{CellBasedSimulation}}} and overrides
 * the {{{ApplyCellPopulationBoundaryConditions()}}} method.
 */
class MyCellBasedSimulation : public CellBasedSimulation<2>
{
/* The first public method is a default constructor, which just calls the base
 * constructor. There are four input argument: a reference to a cell population object,
 * {{{rCellPopulation}}}; an optional flag, {{{deleteCellPopulationAndForceCollection}}},
 * telling the simulation whether to delete the cell population and force collection
 * on destruction to free up memory; and another optional flag, {{{initialiseCells}}},
 * telling the simulation whether to initialise cells. */
public:

    MyCellBasedSimulation(AbstractCellPopulation<2>& rCellPopulation,
                       bool deleteCellPopulationAndForceCollection=false,
                       bool initialiseCells=true)
        : CellBasedSimulation<2>(rCellPopulation, deleteCellPopulationAndForceCollection, initialiseCells)
    {
    }

    /* The second public method overrides {{{ApplyCellPopulationBoundaryConditions()}}}.
     * This method is called during the {{{Solve()}}} method at the end of each
     * timestep, just after the position of each node in the cell population has been
     * updated according to its equation of motion. We iterate over all nodes
     * associated with real cells and update their positions according to any
     * boundary conditions. In our case, this means that if any node moves
     * This method iterates over all cells in the cell_population, and moves
     * any cell whose centre has y coordinate less than 0 or greater than 5 back
     * into the domain. */
    void ApplyCellPopulationBoundaryConditions(const std::vector<c_vector<double,2> >& rOldLocations)
    {
        for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
             cell_iter != mrCellPopulation.End();
             ++cell_iter)
        {
            c_vector<double, 2> cell_location = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

            unsigned node_index = mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<2>* p_node = mrCellPopulation.GetNode(node_index);

            if (cell_location[1] > 5.0)
            {
                p_node->rGetModifiableLocation()[1] = 5.0;
            }
            else if (cell_location[1] < 0.0)
            {
                p_node->rGetModifiableLocation()[1] = 0.0;
            }

            assert(p_node->rGetLocation()[1] <= 5.0);
            assert(p_node->rGetLocation()[1] >= 0.0);
        }
    }
};

/* You only need to include the next block of code if you want to be able to
 * archive (save or load) the cell-based simulation. We start by including a
 * serialization header, then define {{{save_construct_data}}} and
 * {{{load_construct_data}}} methods, which archive the cell-based simulation
 * constructor input argument(s) (in this case, a cell population). */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MyCellBasedSimulation)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const MyCellBasedSimulation * t, const BOOST_PFTO unsigned int file_version)
        {
            // Save data required to construct instance
            const AbstractCellPopulation<2> * p_cell_population = &(t->rGetCellPopulation());
            ar & p_cell_population;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, MyCellBasedSimulation * t, const unsigned int file_version)
        {
            // Retrieve data from archive required to construct new instance
            AbstractCellPopulation<2>* p_cell_population;
            ar >> p_cell_population;

            // Invoke inplace constructor to initialise instance
            ::new(t)MyCellBasedSimulation(*p_cell_population, true, false);
        }
    }
}


/*
 * This completes the code for {{{MyCellBasedSimulation}}}. Note that usually this code
 * would be separated out into a separate declaration in a .hpp file and definition
 * in a .cpp file.
 *
 * EMPTYLINE
 *
 * === The Tests ===
 *
 * EMPTYLINE
 *
 * We now define the test class, which inherits from {{{CxxTest::TestSuite}}}.
 */
class TestCreatingACellBasedSimulationWithBoundaryConditionsTutorial : public CxxTest::TestSuite
{
public:

    /*
     * EMPTYLINE
     *
     * == Testing the cell-based simulation ==
     *
     * EMPTYLINE
     *
     * We now test that our new cell-based simulation is implemented correctly.
     */
    void TestMyCellBasedSimulation() throw(Exception)
    {
        /* The first thing to do, as before, is to set up the start time. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* We use the honeycomb mesh generator to create a honeycomb mesh. */
        HoneycombMeshGenerator generator(5, 5, 0, false);
        /* Get the mesh using the {{{GetMesh()}}} method. */
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * To do this, we can use a static method on the {{{CellsGenerator}}} helper class.
         * The {{{<FixedDurationGenerationBasedCellCycleModel, 2>}}} defines the
         * cell-cycle model and that it is in 2d. We create an empty vector of cells
         * and pass this into the method along with the mesh. The {{{cells}}} vector is
         * populated once the method {{{GenerateBasic}}} is called. */
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        /* Now that we have defined the mesh and cells, we can define the cell population. The
         * constructor takes in the mesh and the cells vector. */
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);


        /* We pass in the cell population into a {{{CellBasedSimulation}}}. */
        MyCellBasedSimulation simulator(cell_population);

        /* We set the output directory and end time. */
        simulator.SetOutputDirectory("TestMyCellBasedSimulation");
        simulator.SetEndTime(30.0);

        /* We must now create one or more force laws, which determine the mechanics of
         * the cell population. For this test, we assume that a cell experiences a force from each
         * neighbour that can be represented as a linear overdamped spring, and so use
         * a {{{GeneralisedLinearSpringForce}}} object. We pass a pointer to this force
         * into a vector. Note that we have called the method {{{SetCutOffLength}}} on the
         * {{{GeneralisedLinearSpringForce}}} before passing it into the collection of force
         * laws - this modifies the force law so that two neighbouring cells do not impose
         * a force on each other if they are located more than 3 units (=3 cell widths)
         * away from each other. This modification is necessary when no ghost nodes are used,
         * for example to avoid artificially large forces between cells that lie close together
         * on the cell population boundary.
         */

        /* We create a force law and pass it to the {{{CellBasedSimulation}}}. */
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(3);
        simulator.AddForce(&linear_force);

        /* Test that the Solve() method does not throw any exceptions: */
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        /* Test that the boundary conditions have been implemented properly: */
        for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
             cell_iter != simulator.rGetCellPopulation().End();
             ++cell_iter)
        {
            unsigned node_index = simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter);
            Node<2>* p_node = simulator.rGetCellPopulation().GetNode(node_index);

            TS_ASSERT_LESS_THAN_EQUALS(p_node->rGetModifiableLocation()[1], 5.0);
            TS_ASSERT_LESS_THAN_EQUALS(0.0, p_node->rGetModifiableLocation()[1]);
        }

        /* We conclude the test by calling {{{Destroy()}}} on any singleton classes. */
        SimulationTime::Destroy();
    }
};

#endif /*TESTCREATINGACELLBASEDSIMULATIONWITHBOUNDARYCONDITIONSTUTORIAL_HPP_*/
