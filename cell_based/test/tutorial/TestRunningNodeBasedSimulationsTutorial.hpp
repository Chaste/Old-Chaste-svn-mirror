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
 * As in previous cell-based Chaste tutorials (UserTutorials/RunningMeshBasedSimulations), we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

/* The following header is usually included in all cell-based test suites. It enables us to write tests where the {{{SimulationTime}}} is handled automatically and simplifies the tests.*/
#include "AbstractCellBasedTestSuite.hpp"

/* The remaining header files define classes that will be used in the cell population
 * simulation test. We encountered some of these header files in 
 * UserTutorials/RunningMeshBasedSimulations. */
#include "CellsGenerator.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the class for storing the spatial information of cells. */
#include "NodesOnlyMesh.hpp"
/* The next header file defines a node-based {{{CellPopulation}}} class.*/
#include "NodeBasedCellPopulation.hpp"
/* The next header file defines a boundary condition to be used in the third test.*/
#include "SphereGeometryBoundaryCondition.hpp"
/* Next, we define the test class. This time we inherit from {{{AbstractCellBasedTestSuite}}} rather than {{{CxxTest::TestSuite}}}.
 * When using this class the singleton objects are set up and destroyed for us:
 * {{{SimulationTime}}} is initialised to zero at the beginning of the test and destroyed at the end of the test;
 * {{{RandomNumberGenerator}}} is re-seeded with zero at the begining and destroyed at the end of the test; and
 * {{{CellPropertyRegistry}}} (which stores {{{CellProperties}}}, you learn about these in a later tutorial UserTutorials/CreatingAndUsingANewCellProperty)  is cleared at the beginning of the test.
 * This makes for cleaner code.
 */
class TestRunningNodeBasedSimulationsTutorial : public AbstractCellBasedTestSuite
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
        /* We no longer need to set up the start time as this is fone in the {{{AbstractCellBasedTestSuite}}}.
         * The first thing we do is generate a nodes only mesh. To do this we first create a {{{MutableMesh}}}
         * to use as a generating mesh.
         * To do this we can use the {{{HoneycombMeshGenerator}}}. This generates a honeycomb-shaped mesh,
         * in which all nodes are equidistant. Here the first and second arguments
         * define the size of the mesh - we have chosen a mesh that is 2 nodes (i.e.
         * cells) wide, and 2 nodes high.
         */
        HoneycombMeshGenerator generator(2, 2);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        /* Once we have a {{{MutableMesh}}} we can generate a {{{NodesOnlyMesh}}} from it using the
         * following commands. Note you can also generate the {{{NodesOnlyMesh}}} from a collection of
         * nodes, see  [class:NodesOnlyMesh NodesOnlyMesh] for details.
         */
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh);

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * To do this, we the `CellsGenerator` helper class, which is templated over the type
         * of cell model required (here {{{StochasticDurationCellCycleModel}}})
         * and the dimension. We create an empty vector of cells and pass this into the
         * method along with the mesh. The second argument represents the size of that the vector
         * {{{cells}}} should become - one cell for each node, the third argument specifies
         * the proliferative type of the cell STEM, TRANSIT or DIFFERENTIATED. */
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(),TRANSIT);

        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
        * In general, this class associates a collection of cells with a mesh.
        * For this test, because we have a {{{NodesOnlyMesh}}}, we use a particular type of
        * cell population called a {{{NodeBasedCellPopulation}}}.
        */
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        /* To run node-based simulations you need to define a cut off length, which
         * defines the connectivity of the nodes by defining a radius of interaction. */
        cell_population.SetMechanicsCutOffLength(1.5);

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory, output multiple and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("NodeBasedMonolayer");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        /* We now pass a force law to the simulation. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* We conclude by calling {{{SimulationTime::Destroy()}}}. To avoid memory leaks, 
         * we also delete any pointers that we created in the test. */
        delete p_mesh;
    }

    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/NodeBasedMonolayer/results_from_time_0}}}.
     * we need to select the 'Cells as circles` option to be able to visualise the cells, as opposed
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
     * In the second test we run a simple node-based simulation in 3D. This is very similar
     * to the 2D test with the dimension template (<2,2> and <2>) changed from 2 to 3 and instead of using a mesh
     * generator we generate the nodes directly.
     */
    void TestSpheroid() throw(Exception)
    {
        /*
         * First, we generate a nodes only mesh. This time we specify the nodes manually by first
         * creating a vector of nodes. */
        std::vector<Node<3>*> nodes;
        /* We then create some nodes to add to this vector. */
        nodes.push_back(new Node<3>(0u,  false,  0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
        nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));
        /* Finally a {{{NodesOnlyMesh}}} is created and the vector of nodes is passed to
         * the {{{ConstructNodesWithoutMesh}}} method. */
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        /*
         * Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * As before, we do this with the `CellsGenerator` helper class (this time with dimension 3).
         */
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(),TRANSIT);

        /* We make a {{{NodeBasedCellPopulation}}} (this time with dimension 3) as before and define the cut off length.
         */
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * (this time with dimension 3) and set the output directory, output multiple and end time. */
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("NodeBasedSpheroid");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        /* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        simulator.AddForce(p_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    /*
     * EMPTYLINE
     *
     * To visualize the results, use Paraview. See the UserTutorials/VisualizingWithParaview tutorial for more information.
     *
     * Load the file {{{/tmp/$USER/testoutput/NodeBasedSpheroid/results_from_time_0/results.pvd}}},
     * and add spherical glyphs to represent cells.
     *
     * EMPTYLINE
     *
     * == Test 3 - a node-based simulation on a restricted geometry ==
     *
     * EMPTYLINE
     *
     * In the third test we run a node-based simulation restricted to the surface of a sphere.
     */
    void TestOnSurfaceOfSphere() throw(Exception)
    {
        /*
         * We begin with exactly the same code as the previous test: we create a cell population
         * from a mesh and vector of cells, and use this in turn to create
         * a simulation object.
         */

        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u,  false,  0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
        nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(),TRANSIT);

        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("NodeBasedOnSphere");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        /* As before, we create a linear spring force and pass it to the simulation object. */
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        simulator.AddForce(p_force);

        /*
         * This time we create a {{{CellPopulationBoundaryCondition}}} and pass this to
         * the {{{OffLatticeSimulation}}}. Here we use a {{{SphereGeometryBoundaryCondition}}}
         * which restricts cells to lie on a sphere (in 3D) or circle (in 2D).
         *
         * For a list of possible boundary conditions see subclasses of {{{AbstractCellPopulationBoundaryCondition}}}.
         * These can be found in the inheritance diagram, here, [class:AbstractCellPopulationBoundaryCondition AbstractCellPopulationBoundaryCondition].
         * Note that some of these boundary conditions are not compatible with node-based simulations see the specific class documentation for details,
         * if you try to use an incompatible class then you will receive a warning.
         *
         * First we set the centre (0,0,1) and radius of the sphere (1).
         */
        c_vector<double,3> centre = zero_vector<double>(3);
        centre(2) = 1.0;
        double radius = 1.0;
        /* We then make a pointer to the boundary condition using the MAKE_PTR_ARGS macro, and pass
         * it to the {{{OffLatticeSimulation}}}. */
        MAKE_PTR_ARGS(SphereGeometryBoundaryCondition<3>, p_boundary_condition, (&cell_population, centre, radius));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
    /*
     * EMPTYLINE
     *
     * To visualize the results, use Paraview. See the UserTutorials/VisualizingWithParaview tutorial for more information.
     *
     * Load the file {{{/tmp/$USER/testoutput/NodeBasedOnSphere/results_from_time_0/results.pvd}}},
     * and add spherical glyphs to represent cells.
     *
     * EMPTYLINE
     */
};

#endif /* TESTRUNNINGNODEBASEDSIMULATIONSTUTORIAL_HPP_ */
