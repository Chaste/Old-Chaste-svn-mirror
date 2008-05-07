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
 * = Introduction =
 * 
 * In this tutorial we show how Chaste is used to run discrete crypt 
 * simulations.
 * 
 * The first thing that needs to be done, when writing any Chaste test,
 * is to include the following header
 */
#include <cxxtest/TestSuite.h>
/* The following have to be included, for technical reasons....... */
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
/* This header defines a helper class that is useful for generating a 
 * vector of cells */
#include "CellsGenerator.hpp"
/* These are the classes that will be used in these tests */ 
#include "MeshBasedTissueWithGhostNodes.hpp"
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
     * = Test 1 - a basic crypt simulation =
     * 
     * EMPTYLINE
     * 
     * In the first test, we solve run a simple crypt simulation, where we use
     * a cylindrical mesh, fixed cell cycles, and sloughing at the top of the 
     * crypt.
     */ 
    void TestCryptFixedCellCycle() throw(Exception)
    {
        /* The first thing that has to be done, in '''all tissue simulations''', in
         * the following. First, the start time is set (has to be done), and the 
         * cancer parameters are reset (ought to be done). (These are singleton classes, 
         * which means there is one and only one of these objects instantiated at any time, 
         * and that that single object is accessible from anywhere in the code (this means 
         * that time does not have to be passed about)).
         */
        SimulationTime::Instance()->SetStartTime(0.0);
        CancerParameters::Instance()->Reset();
        
        /* Next, we generate a mesh. The basic Chaste mesh is {{{ConformingTetrahedralMesh}}}. To enforce 
         * periodicity, {{{Cylindrical2dMesh}}} was created, which is basically a normal mesh
         * that knows how to keep itself periodic. To create a {{{Cylindrical2dMesh}}}, we can use
         * the {{{HoneycombMeshGenerator}}} which generates a honeycomb-shaped mesh, ie all nodes
         * equidistant. Here we create a honeymesh which is 6 nodes (ie cells) wide, and 9 cells long 
         * (representing a mouse crypt). The {{{true}}} indicates that we want a cylindrical mesh. 
         * ''Ghost Nodes:'' Ghost nodes are 'fake' nodes around the edge of a mesh, which are needed
         * when cell-connectivity is defined using a triangulation method, and when we don't want
         * cells on the edge of a mesh that are far away getting connected. The '2' below says we want a 
         * layer of ghost nodes of thickness around the mesh (or actually, above and below the mesh, as
         * it is periodic).
         */
        HoneycombMeshGenerator generator(6, 9, 2, true); // params are: cells across, cells up, thickness of ghost layer, whether to be cylindrical
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        /* Now, we create the ''Tissue'', which, in general, is a collection of cells
         * together with nodes or a mesh. First, we need to define a {{{std::vector}}} of 
         * {{{TissueCell}}}s. To do this, we can use a static method on the {{{CellsGenerator}}} 
         * class. Note that the {{{<2>}}} below denotes the dimension. Here, we create an empty 
         * vector of cells, pass that into the method along with the 
         * mesh, 'FIXED' saying we want cells with a fixed cell cycle, and 'true' indicating
         * we want random birth times for the cells. The {{{cells}}} vector will be populated
         * once the method is called. */
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);

        /* Now we have a mesh, a set of cells to go with it, and ghost nodes indices, we can 
         * create the tissue, which for this test is of type {{{MeshBasedTissueWithGhostNodes}}}. 
         */
        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, ghost_node_indices);
        
        /* That's most of the setup. Now we just define the Simulation object, passing in the tissue. */
        CryptSimulation2d simulator(tissue);
        
        /* Set the output directory on the simulator (NOTE: this is relative to
         * "/tmp/<USER_NAME>/testoutput"), and the end time (NOTE: in hours). 
         */
        simulator.SetOutputDirectory("CryptTutorialFixedCellCycle");
        simulator.SetEndTime(1);
        /* Note: for longer simulations, you may not want to output the results
         * every timestep, in which case you can do the following, to print
         * every 10 timesteps instead. (The timestep which is used, internally in 
         * the simulator, is about 0.01hrs, so this will print every ~0.1hrs).
         */
        //simulator.SetSamplingTimestepMultiple(10);

        /* Before calling solve, we want to add a cell killer. These are objects
         * which provide rules on when cells should be killed. We will use
         * a {{{SloughingCellKiller}}}, which kills cells above a certain height.
         */
        SloughingCellKiller killer(&tissue);
        simulator.AddCellKiller(&killer);

        /* To solve, just call Solve */
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
     * To visualise the results, open a new terminal, cd to the Chaste directory,
     * then cd to 'anim'. Then do: {{{java Visualize2dCells /tmp/<USER_NAME>/testoutput/CryptTutorialFixedCellCycle/results_from_time_0}}}. 
     * You may have to do: {{{javac Visualize2dCells.java}}} beforehand to create the 
     * java executable.
     * 
     * EMPTYLINE
     * 
     * = Test 2 - using Wnt based cell-cycle models =
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
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, WNT, true);
        
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
        
        /* Create a simulator as before (except setting a different output directory). */
        CryptSimulation2d simulator(tissue);
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
