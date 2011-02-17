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
#ifndef TESTRUNNINGCRYPTSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGCRYPTSIMULATIONSTUTORIAL_HPP_
/*
 * = Examples showing how to run crypt simulations on periodic meshes with different cell cycle models =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to simulate a cylindrical model of an
 * intestinal crypt. Full details of the computational model can be found in the paper by
 * van Leeuwen ''et al'' (2009) [doi:10.1111/j.1365-2184.2009.00627.x].
 *
 * The first thing to do is include the following header, which allows us
 * to use certain methods in our test. This header file should be included
 * in any Chaste test.
 */
#include <cxxtest/TestSuite.h>
/*
 * Any test in which the {{{GetIdentifier()}}} method is used, even via the main
 * `cell_based` code (through calls to {{{AbstractCellPopulation}}} output methods),
 * must include {{{CheckpointArchiveTypes.hpp}}} or {{{CellBasedSimulationArchiver.hpp}}}
 * as the first Chaste header included.
 */
#include "CheckpointArchiveTypes.hpp" 

/* The next header file defines a helper class for generating cells for crypt simulations. */
#include "CryptCellsGenerator.hpp"
/*
 * The next two header files define two different types of cell-cycle model.
 * In a {{{FixedDurationGenerationBasedCellCycleModel}}}, the duration of each phase
 * of the cell cycle is fixed. In a {{{WntCellCycleModel}}}, the duration of G1 phase
 * is determined by a system of nonlinear ODEs describing a cell's response to Wnt,
 * a secreted cell–cell signalling molecule that is known to play a key role in cell
 * proliferation in the crypt.
 */
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
/* The next header file defines a helper class for generating a suitable mesh
 * for the crypt simulation, such that the cell corresponding to each node is initially
 * in mechanical equilibrium with its neighours. */
#include "CylindricalHoneycombMeshGenerator.hpp"
/* The next header file defines a {{{CellPopulation}}} class that uses a tetrahedral mesh, and allows
 * for the inclusion of 'ghost nodes'. These are nodes in the mesh that do not correspond
 * to cells; instead they help ensure that a sensible Delaunay triangulation is generated
 * at each timestep (since the triangulation algorithm requires a convex hull). */
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
/*
 * The next header file defines a force law, based on a linear spring, for describing
 * the mechanical interactions between neighbouring cells in the crypt.
 */
#include "GeneralisedLinearSpringForce.hpp"
/*
 * The next header file defines the class that simulates the evolution of a {{{CellPopulation}}},
 * specialized to deal with the cylindrical crypt model.
 */
#include "CryptSimulation2d.hpp"
/*
 * The next header file defines a Wnt singleton class, which (if used) deals with the
 * imposed Wnt gradient in our crypt model. This affects cell proliferation in the case
 * where we construct each cell with a {{{WntCellCycleModel}}.
 */
#include "WntConcentration.hpp"
/*
 * The final header file defines a cell killer class, which implements sloughing of cells
 * into the lumen once they reach the top of the crypt.
 */
#include "SloughingCellKiller.hpp"
/*
 * Next, we define the test class, which inherits from {{{CxxTest::TestSuite}}}
 * and defines some test methods.
 */
class TestRunningCryptSimulationsTutorial : public CxxTest::TestSuite
{
public:
    /* EMPTYLINE
     *
     * == Test 1 - a basic crypt simulation ==
     *
     * EMPTYLINE
     *
     * In the first test, we demonstrate how to create a crypt simulation using a
     * cylindrical mesh, with each cell progressing through a fixed cell-cycle model,
     * and sloughing enforced at the top of the crypt.
     */
    void TestCryptFixedCellCycle() throw(Exception)
    {
        /* As in '''all''' cell-based simulations, we must first set the start time.
         */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* Next, we generate a mesh. The basic Chaste mesh is {{{TetrahedralMesh}}}.
         * To enforce periodicity at the left- and right-hand sides of the mesh, we
         * use a subclass called {{{Cylindrical2dMesh}}}, which has extra methods for
         * maintaining periodicity. To create a {{{Cylindrical2dMesh}}}, we can use
         * the {{{CylindricalHoneycombMeshGenerator}}}. This generates a periodic honeycomb-shaped mesh,
         * in which all nodes are equidistant to their neighbours. Here the first and second arguments
         * define the size of the mesh - we have chosen a mesh that is 6 nodes (i.e.
         * cells) wide, and 9 nodes high. The third argument indicates that we require
         * a double layer of ghost nodes around the mesh (technically, just above
         * and below the mesh, since it is periodic). We call {{{GetCylindricalMesh()}}} on the {{{CylindricalHoneycombMeshGenerator}}} to
         * return our {{{Cylindrical2dMesh}}}, and call {{{ GetCellLocationIndices()}}}
         * to return a {{{std::vector}}} of indices of nodes in the mesh that correspond to real cells (as opposed
         * to ghost nodes).
         */
        CylindricalHoneycombMeshGenerator generator(6, 9, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        /*
         * Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * To do this, we the `CryptCellsGenerator` helper class, which is templated over the type
         * of cell-cycle model required (here {{{FixedDurationGenerationBasedCellCycleModel}}})
         * and the dimension. We create an empty vector of cells and pass this into the
         * method {{{Generate()}}} along with the mesh. The third argument 'true' indicates that the cells
         * should be assigned random birth times, to avoid synchronous division. The
         * {{{cells}}} vector is populated once the method {{{Generate()}}} is
         * called. Note that we only ever deal with shared pointers to cells, named {{{CellPtr}}}s.
         */
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        /* 
         * Now we have a mesh, a set of cells to go with it, and a vector of node indices
         * corresponding to real cells, we can create a {{{CellPopulation}}} object. In general,
         * this class associates a collection of cells with a set of nodes or a mesh.
         * For this test, because we have a mesh and ghost nodes, we use a particular type of
         * cell population called a {{{MeshBasedCellPopulationWithGhostNodes}}}.
         */
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        /*
         * Next we use the ''CellPopulation'' object to construct a {{{CryptSimulation2d}}} object,
         * which will be used to simulate the crypt model. */
        CryptSimulation2d simulator(cell_population);

        /*
         * We must set the output directory on the simulator (relative to
         * "/tmp/<USER_NAME>/testoutput") and the end time (in hours).
         */
        simulator.SetOutputDirectory("CryptTutorialFixedCellCycle");
        simulator.SetEndTime(1);
        /*
         * For longer simulations, you may not want to output the results
         * every time step. In this case you can use the following method,
         * to print results every 10 time steps instead. As the time step
         * used by the simulator, is 30 seconds, this method will cause the
         * simulator to print results every 5 minutes.
         */
        simulator.SetSamplingTimestepMultiple(10);

        /*
         * Before running the simulation, we must add one or more force laws, which determine the mechanical
         * behaviour of the cell population. For this test, we use a {{{GeneralisedLinearSpringForce}}}, which assumes
         * that every cell experiences a force from each of its neighbours that can be represented as a linear overdamped
         * spring.
         */
        GeneralisedLinearSpringForce<2> linear_force;
        simulator.AddForce(&linear_force);

        /*
         * We also add a cell killer to the simulator. This object
         * dictates under what conditions cells die. For this test, we use
         * a {{{SloughingCellKiller}}}, which kills cells above a certain
         * height (passed as an argument to the constructor).
         */
        double crypt_height = 8.0;
        SloughingCellKiller<2> killer(&cell_population, crypt_height);
        simulator.AddCellKiller(&killer);

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
     * Finally, to visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/CryptTutorialFixedCellCycle/results_from_time_0}}}.
     * You may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
     * java executable.
     *
     * EMPTYLINE
     *
     * == Test 2 - a Wnt-dependent crypt simulation ==
     *
     * EMPTYLINE
     *
     * The next test is very similar to Test 1, except that instead of
     * using a fixed cell-cycle model, we use a Wnt-dependent cell cycle model,
     * with the Wnt concentration varying within the crypt in a predefined manner.
     */
    void TestCryptWntCellCycle() throw(Exception)
    {
        /* First re-initialize time to zero and reseed the random number generator. */
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(0);

        /* Create a cylindrical mesh, and get the cell location indices, exactly as before. */
        CylindricalHoneycombMeshGenerator generator(6, 9, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        /* Create the cells, using the same method as before. Here, though, we use a {{{WntCellCycleModel}}}.*/
        std::vector<CellPtr> cells;
        CryptCellsGenerator<WntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        /* Create the cell_population, as before. */
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        /*
         * Set the height of the crypt. As well as passing this variable into the {{{sloughingCellKiller}}},
         * we will pass it to the {{{WntConcentration}}} object (see below).
         */
        double crypt_height = 8.0;

        /*
         * When using a {{{WntCellCycleModel}}}, we need a way of telling each cell what the Wnt concentration
         * is at its location. To do this, we set up a {{{WntConcentration}}} object. Like {{{SimulationTime}}},
         * {{{WntConcentration}}} is a singleton class, so when instantiated it is accessible from anywhere in
         * the code (and in particular, all cells and cell cycle models can access it). We need to say what 
         * the profile of the Wnt concentation should be up the crypt: here, we say it is {{{LINEAR}}} (linear
         * decreasing from 1 to 0 from the bottom of the crypt to the top). We also need to inform the
         * {{{WntConcentration}}} of the cell population and the height of the crypt.
         */
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_height);

        /* Create a simulator as before (except setting a different output directory). */
        CryptSimulation2d simulator(cell_population);
        simulator.SetOutputDirectory("CryptTutorialWntCellCycle");
        simulator.SetEndTime(1);

        /* As before, we create a force law and cell killer and pass these objects to the simulator, then call
         * Solve(). */
        GeneralisedLinearSpringForce<2> linear_force;
        simulator.AddForce(&linear_force);
        SloughingCellKiller<2> killer(&cell_population, crypt_height);
        simulator.AddCellKiller(&killer);

        simulator.Solve();

        /*
         * Finally, we must tidy up by destroying the {{{SimulationTime}}} and the {{{WntConcentration}}}
         * singleton objects. This avoids memory leaks occurring. */
        WntConcentration<2>::Destroy();
        SimulationTime::Destroy();
    }
};
    /*
     * EMPTYLINE
     *
     * The results of this test can be visualized as in Test 1, just with the different output directory.
     */
#endif /*TESTRUNNINGCRYPTSIMULATIONSTUTORIAL_HPP_*/
