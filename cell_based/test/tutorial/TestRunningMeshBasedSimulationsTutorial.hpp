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

#ifndef TESTRUNNINGMESHBASEDSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGMESHBASEDSIMULATIONSTUTORIAL_HPP_

/*
 * = Examples showing how to create, run and visualize mesh-based simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to create, run and visualize mesh-based simulations.
 * Full details of the mathematical model can be found in van Leeuwen ''et al.'' (2009) [doi:10.1111/j.1365-2184.2009.00627.x].
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

/* The remaining header files define classes that will be used in the cell population
 * simulation test. We have encountered some of these header files in previous cell-based
 * Chaste tutorials. */
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the cell cycle model. */
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
/* The next header file defines a helper class for generating a suitable mesh. */
#include "HoneycombMeshGenerator.hpp"
/* The next header file defines the class that simulates the evolution of a {{{CellPopulation}}}
 * for a vertex mesh. */
#include "OffLatticeSimulation.hpp"
/* The next header files define a mesh-based {{{CellPopulation}}} with and without ghost nodes class.*/
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
/* The next header file defines a force law for describing the mechanical interactions
 * between neighbouring cells in the cell population.
 */
#include "GeneralisedLinearSpringForce.hpp"
/* Next, we define the test class, which inherits from {{{CxxTest::TestSuite}}}
 * and defines some test methods.
 */
class TestRunningMeshBasedSimulationsTutorial : public CxxTest::TestSuite
{
public:
    /* EMPTYLINE
    *
    * == Test 1 - a basic mesh-based simulation ==
    *
    * EMPTYLINE
    *
    * In the first test, we run a simple mesh-based simulation, in which we create a monolayer
    * of cells, using a mutable mesh. Each cell is assigned a stochastic cell-cycle model.
    */
    void TestMonolayer() throw(Exception)
    {
        /* As in previous cell-based Chaste tutorials, we begin by setting up the start time. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* Next, we generate a mutable mesh. To create a {{{MutableMesh}}}, we can use
        * the {{{HoneycombMeshGenerator}}}. This generates a honeycomb-shaped mesh,
        * in which all nodes are equidistant. Here the first and second arguments
        * define the size of the mesh - we have chosen a mesh that is 2 nodes (i.e.
        * cells) wide, and 2 nodes high.
        */
        HoneycombMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
        * To do this, we the `CellsGenerator` helper class, which is templated over the type
        * of cell model required (here {{{StiochasticDurationGenerationBasedCellCycleModel}}})
        * and the dimension. We create an empty vector of cells and pass this into the
        * method along with the mesh. The second argument represents the size of that the vector
        * {{{cells}}} should become - one cell for each node, the third argument specifies
        * the proliferative type of the cell STEM TRANSIT or DIFFERENTIATED. */
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(),TRANSIT);

        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
        * In general, this class associates a collection of cells with a set of elements or a mesh.
        * For this test, because we have a {{{MutableMesh}}}, we use a particular type of
        * cell population called a {{{MeshBasedCellPopulation}}}.
        */
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* We then pass in the cell population into a {{{OffLatticeSimulation}}},
         * and set the output directory and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("MeshBasedMonolayer");
        simulator.SetEndTime(10.0);

        /*
         * For longer simulations, we may not want to output the results
         * every time step. In this case we can use the following method,
         * to print results every 12 time steps instead. As the time step
         * used by the simulator, is 30 seconds, this method will cause the
         * simulator to print results every 6 minutes.
         */
        simulator.SetSamplingTimestepMultiple(12);

        /* We must now create one or more force laws, which determine the mechanics of the centres
        * of each cell in a cell population. For this test, we use one force law, based on the
        * spring based model, and pass it to the {{{OffLatticeSimulation}}}.
        * For a list of possible update rules see subclasses of {{{AbstractForce}}}.
        * These can be found in the inheritance diagram, here, [class:AbstractForce AbstractForce].
        * Note that some of these forces are not compatible with node based simulations see the specific class documentation for details,
        * if you try to use an incompatible class then you will receive a warning.
        */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* {{{SimulationTime::Destroy()}}} '''must''' be called at the end of the test.
        * If not, when {{{SimulationTime::Instance()->SetStartTime(0.0);}}} is called
        * at the beginning of the next test in this file, an assertion will be triggered.
        */
        SimulationTime::Destroy();
    }

    /*
    * EMPTYLINE
    *
    * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
    * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/MonolayerFixedCellCycle/results_from_time_0}}}.
    * We may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
    * java executable.
    *
    * EMPTYLINE
    *
    * == Test 2 - a basic mesh-based simulation with ghost nodes ==
    *
    * EMPTYLINE
    *
    * In the first test, we run a simple mesh-based simulation, in which we create a monolayer
    * of cells, using a mutable mesh. Each cell is assigned a stochastic cell-cycle model.
    */
    void TestMonolayerWithGhostNodes() throw(Exception)
    {

    }

};

#endif /* TESTRUNNINGMESHBASEDSIMULATIONSTUTORIAL_HPP_ */
