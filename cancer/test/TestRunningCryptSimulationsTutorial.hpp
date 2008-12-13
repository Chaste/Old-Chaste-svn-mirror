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
 * In this tutorial we show how Chaste is used to run crypt simulations. 
 * Full details of the computational model can be found in the paper by 
 * van Leeuwen ''et al'' (to appear in Cell Prolif.)  
 *
 * The first thing to do is include the following header, which allows us 
 * to use certain methods in our test (this header file should be included 
 * in any Chaste test):
 */
#include <cxxtest/TestSuite.h>
/* The next two header files define helper classes for generating a vector of
 * cells with fixed, and Wnt-dependent, cell cycle models: */
#include "FixedCellCycleModelCellsGenerator.hpp"
#include "WntCellCycleModelCellsGenerator.hpp"
/* This header file defines a helper class for generating a suitable mesh: */
#include "HoneycombMeshGenerator.hpp"
/* These are the classes that will be used in these tests */
#include "MeshBasedTissueWithGhostNodes.hpp"
#include "MeinekeInteractionForce.hpp"
#include "CryptSimulation2d.hpp"
#include "WntConcentration.hpp"
#include "SloughingCellKiller.hpp"

/* Next, we define the test class, which inherits from {{{CxxTest::TestSuite}}}
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
     * In the first test, we run a simple crypt simulation, in which we use
     * a cylindrical mesh, give each cell a fixed cell cycle model, and enforce 
     * sloughing at the top of the crypt.
     */
    void TestCryptFixedCellCycle() throw(Exception)
    {
        /* As in '''all''' tissue simulations, we must first set the start time. 
         * In addition, it is advisable to reset the values of all model parameters. 
         * {{{SimulationTime}}} and {{{CancerParameters}}} are ''singleton'' classes; this 
         * means that one and only one of each of these objects is instantiated at 
         * any time, and that that single object is accessible from anywhere in the 
         * code. As a result, we do not need to keep passing round the current time or 
         * model parameter values.
         */
        SimulationTime::Instance()->SetStartTime(0.0);
        CancerParameters::Instance()->Reset();

        /* Next, we generate a mesh. The basic Chaste mesh is {{{TetrahedralMesh}}}. 
         * To enforce periodicity at the left and right hand sides of the mesh, we 
         * use a sublcass called {{{Cylindrical2dMesh}}}, which has extra methods for
         * maintaining periodicity. To create a {{{Cylindrical2dMesh}}}, we can use
         * the {{{HoneycombMeshGenerator}}}. This generates a honeycomb-shaped mesh, 
         * in which all nodes are equidistant. Here the first and second arguments 
         * define the size of the mesh - we have chosen a mesh that is 6 nodes (i.e. 
         * cells) wide, and 9 nodes high. The third argument indicates that we require 
         * a double layer of ghost nodes around the mesh (technically, just above 
         * and below the mesh, since it is periodic).
         */
        HoneycombMeshGenerator generator(6, 9, 2, true); // params are: cells across, cells up, thickness of ghost layer, whether to be cylindrical
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{TissueCell}}}s.
         * To do this, we can use a static method on the {{{FixedCellCycleModelCellsGenerator}}}
         * helper class. The {{{<2>}}} below denotes the dimension. We create an empty vector 
         * of cells and pass this into the method along with the mesh. The third argument 
         * 'true' indicates that the cells should be assigned random birth times, to avoid 
         * synchronous division. The {{{cells}}} vector is populated once the method 
         * {{{GenerateForCrypt}}} is called. */
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);

        /* Now we have a mesh, a set of cells to go with it, and ghost nodes indices, 
         * we can create a ''Tissue''. In general, this class associates a collection 
         * of cells with a set of nodes or a mesh. For this test, because we have a 
         * mesh and ghost nodes, we use aparticular type of tissue called a 
         * {{{MeshBasedTissueWithGhostNodes}}}.
         */
        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, ghost_node_indices);

        /* We must now create one or more force laws, which determine the mechanics of
         * the tissue. For this test, we assume that a cell experiences a force from each
         * neighbour that can be represented as a linear overdamped spring. Since this 
         * model was first proposed in the context of crypt modelling by Meineke ''et al'' 
         * (Cell Prolif. 34:253-266, 2001), we call this object a 
         * {{{MeinekeInteractionForce}}}. We pass a pointer to this force into a vector.
         */
        MeinekeInteractionForce<2> meineke_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&meineke_force);

        /* Now we define the tissue simulation object, passing in the tissue and collection
         * of force laws: */
        CryptSimulation2d simulator(tissue, force_collection);

        /* Set the output directory on the simulator (relative to
         * "/tmp/<USER_NAME>/testoutput") and the end time (in hours).
         */
        simulator.SetOutputDirectory("CryptTutorialFixedCellCycle");
        simulator.SetEndTime(1);
        /* For longer simulations, you may not want to output the results
         * every time step. In this case you can use the following method, 
         * to print results every 10 time steps instead. As the time step 
         * used by the simulator, is 30 s, this method will cause the 
         * simulator to print results every 5 min.
         */
        //simulator.SetSamplingTimestepMultiple(10);

        /* Before running the simulation, we add a cell killer. This object 
         * dictates conditions under which cells die. For this test, we use
         * a {{{SloughingCellKiller}}}, which kills cells above a certain height.
         */
        SloughingCellKiller killer(&tissue);
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
     * To visualize the results, open a new terminal, cd to the Chaste directory,
     * then cd to 'anim'. Then do: {{{java Visualize2dCells /tmp/<USER_NAME>/testoutput/CryptTutorialFixedCellCycle/results_from_time_0}}}.
     * You may have to do: {{{javac Visualize2dCells.java}}} beforehand to create the
     * java executable.
     *
     * EMPTYLINE
     *
     * == Test 2 - using Wnt based cell-cycle models ==
     *
     * EMPTYLINE
     *
     * The next test is very similar (almost identical in fact), except instead of
     * using a fixed cell cycle model, we use a Wnt (a protein) based cell cycle model,
     * with the Wnt concentration depending on the position of the cell within the crypt.
     */
    void TestCryptWntCellCycle() throw(Exception)
    {
        /* First reinitialise time to 0, and reset the cancer parameters, again. */
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(0);
        CancerParameters::Instance()->Reset();

        /* Create a cylindrical mesh, and get the ghost node indices, exactly as before. */
        HoneycombMeshGenerator generator(6, 9, 2, true);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        /* Create the cells, using the same method as before. Here, though, we pass
         * in 'WNT' as the third parameters, saying the cells should have a
         * Wnt based cell-cycle. This is an ODE based cell cycle. */
        std::vector<TissueCell> cells;
        WntCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);

        /* Create the tissue, as before. */
        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, ghost_node_indices);

        /* The other change needed: Cells with a Wnt-based cell cycle need to know
         * the concentration of Wnt wherever they are. To do this, we set up a {{{WntConcentration}}}
         * class. This is another singleton class (ie accessible from anywhere), so all
         * cells and cell cycle models can access it. We need to say what the profile of the
         * Wnt concentation should be - here, we say it is linear (linear decreasing from 1 to 0
         * from the bottom of the crypt to the top). We also need to inform the {{{WntConcentration}}}
         * of the tissue.*/
        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(tissue);

        /* '''TODO''' add comment */
        MeinekeInteractionForce<2> meineke_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&meineke_force);

        /* Create a simulator as before (except setting a different output directory). */
        CryptSimulation2d simulator(tissue,force_collection);
        simulator.SetOutputDirectory("CryptTutorialWntCellCycle");
        simulator.SetEndTime(1);

        /* Create a killer, as before. */
        SloughingCellKiller killer(&tissue);
        simulator.AddCellKiller(&killer);

        /* Solve. */
        simulator.Solve();

        /* Destroy the time, and the {{{WntConcentration}}} object. The solution can be visualised using the
         * Visualizer as before, just with the different output directory. */
        WntConcentration::Destroy();
        SimulationTime::Destroy();
    }
};
#endif /*TESTRUNNINGCRYPTSIMULATIONSTUTORIAL_HPP_*/
