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

#ifndef TESTRUNNINGNODEBASEDSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGNODEBASEDSIMULATIONSTUTORIAL_HPP_

/*
 * = Examples showing how to create, run and visualize node-based simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to create, run and visualize node-based simulations.
 * Full details of the mechanical model can be found in Pathamathan et al "A computational study of
 * discrete mechanical tissue models", Physical Biology. Vol. 6. No. 3. 2009.. DOI (10.1088/1478-3975/6/3/036001).
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
/* The next header file defines the class for storing the spatial information of cells. */
#include "NodesOnlyMesh.hpp"
/* The next header file defines the class that simulates the evolution of a {{{CellPopulation}}}
 * for a vertex mesh. */
#include "OffLatticeSimulation.hpp"
/* The next header file defines a vertex-based {{{CellPopulation}}} class.*/
#include "NodeBasedCellPopulation.hpp"
/* The next header file defines a force law for describing the mechanical interactions
 * between neighbouring cells in the cell population.
 */
#include "GeneralisedLinearSpringForce.hpp"
/* Next, we define the test class, which inherits from {{{CxxTest::TestSuite}}}
 * and defines some test methods.
 */
class TestRunningNodeBasedSimulationsTutorial : public CxxTest::TestSuite
{
public:
    /* EMPTYLINE
    *
    * == Test 1 - a basic node-based simulation ==
    *
    * EMPTYLINE
    *
    * In the first test, we run a simple node-based simulation, in which we create a monolayer
    * of cells, using a nodes only mesh. Each cell is assigned a stochastic cell-cycle model.
    */
    void TestMonolayer() throw(Exception)
    {
        /* As in previous cell-based Chaste tutorials, we begin by setting up the start time. */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* Next, we generate a nodes only mesh. To do this we first create a {{{MutableMesh}}}
        * to use as a generating mesh.
        * To do this we can use the {{{HoneycombMeshGenerator}}}. This generates a honeycomb-shaped mesh,
        * in which all nodes are equidistant. Here the first and second arguments
        * define the size of the mesh - we have chosen a mesh that is 2 nodes (i.e.
        * cells) wide, and 2 nodes high.
        */
        HoneycombMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        /* Once we have a {{{MutableMesh}}} we can generate a {{{NodesOnlyMesh}}} from it using the
         * following commands. Note you can also generate the {{{NodesOnlyMesh}}} from a collection of
         * nodes, see  [class:NodesOnlyMesh NodesOnlyMesh] for details.
        */
        NodesOnlyMesh<2>* p_mesh = new(NodesOnlyMesh<2>);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh);


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
        * For this test, because we have a {{{NodesOnlyMesh}}}, we use a particular type of
        * cell population called a {{{NodeBasedCellPopulation}}}.
        */
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        /* To run node based simulations you need to define a cut of length, which
         * defines the connectivity of the nodes by defining a radius of interaction.
         */
        cell_population.SetMechanicsCutOffLength(1.5);

        /* We then pass in the cell population into a {{{OffLatticeSimulation}}},
         * and set the output directory and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("NodeBasedMonolayer");
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
    * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/NodeBasedMonolayer/results_from_time_0}}}.
    * You will need to select the 'Cells as circles` option to be able to visualise the cells, as opposed
    * to just the centres.
    * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
    * java executable.
    *
    * EMPTYLINE
    *
    * == Test 2 - a basic node-based simulation in 3d ==
    *
    * EMPTYLINE
    *
    * In the second test we run a simple node based simulation in 3d
    */
    void TestSpheroid() throw(Exception)
    {
    }

    /*
    * EMPTYLINE
    *
    * To visualize the results, use paraview
    *
    * EMPTYLINE
    *
    * == Test 3 - a node-based simulation on a restricted geometry ==
    *
    * EMPTYLINE
    *
    * In the third test we run a node based simulation restricted to the surface of a sphere.
    */
    void TestOnSurfaceOfSphere() throw(Exception)
    {

    }
    /*
    * EMPTYLINE
    *
    * To visualize the results, use paraview
    *
    * EMPTYLINE
    */
};

#endif /* TESTRUNNINGNODEBASEDSIMULATIONSTUTORIAL_HPP_ */
