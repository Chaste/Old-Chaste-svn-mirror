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

#ifndef TESTRUNNINGCRYPTSIMULATIONSWITHMUTATIONSTUTORIAL_HPP_
#define TESTRUNNINGCRYPTSIMULATIONSWITHMUTATIONSTUTORIAL_HPP_

/*
 * = Examples showing how to run crypt simulations with various mutations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * This tutorial assumes you have already read UserTutorials/RunningMeshBasedCryptSimulations
 * and UserTutorials/RunningVertexBasedCryptSimulations.
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to simulate a cylindrical model of an
 * intestinal crypt with mutations in both mesh and vertex based simulations.
 * Full details of the computational model can be found in the paper by
 * Osborne ''et al.'' (2010) [10.1098/rsta.2010.0173].
 *
 * The first thing to do is include the following header file, which allows us
 * to use certain methods in our test. This header file must be included
 * in any Chaste test.
 */
#include <cxxtest/TestSuite.h>
/*
 * Any test in which the {{{GetIdentifier()}}} method is used, even via the main
 * `cell_based` code (through calls to {{{AbstractCellPopulation}}} output methods),
 * must also include {{{CheckpointArchiveTypes.hpp}}} or {{{CellBasedSimulationArchiver.hpp}}}
 * as the first Chaste header file.
 */
#include "CheckpointArchiveTypes.hpp"

/* The next header includes the Boost shared_ptr smart pointer, and defines some useful
 * macros to save typing when using it. */
#include "SmartPointers.hpp"

/* The next header file defines a helper class for generating cells for crypt simulations. */
#include "CryptCellsGenerator.hpp"
/*
 * The next two header files define two different types of cell-cycle model.
 * In a {{{FixedDurationGenerationBasedCellCycleModel}}}, the duration of each phase
 * of the cell cycle is fixed. In a {{{WntCellCycleModel}}}, the duration of a cell's G1 phase
 * is determined by a system of nonlinear ODEs describing a cell's response to the local
 * concentration of Wnt,
 * a secreted cell–cell signalling molecule that is known to play a key role in cell
 * proliferation in the crypt. In our crypt simulations, we impose a fixed gradient of
 * Wnt up the axis of the crypt.
 */
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
/* The next header file defines a helper class for generating a suitable triangular mesh
 * for the crypt simulation, such that the cell corresponding to each node is initially
 * in mechanical equilibrium with its neighours and periodic boundary conditions are applied
 * at the left- and right-hand sides of the mesh (hence the "cylindrical"). */
#include "CylindricalHoneycombMeshGenerator.hpp"
/* The next header file defines a {{{CellPopulation}}} class that uses a triangular mesh, and allows
 * for the inclusion of 'ghost nodes'. These are nodes in the mesh that do not correspond
 * to cells; instead they help ensure that a sensible Delaunay triangulation is generated
 * at each timestep. This is because the triangulation algorithm requires a convex hull. */
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
/*
 * The next header file defines a force law, based on a linear spring, for describing
 * the mechanical interactions between neighbouring cells in the crypt.
 */
#include "GeneralisedLinearSpringForce.hpp"
/*
 * The next header file defines the class that simulates the evolution of a {{{CellPopulation}}},
 * specialized to deal with the cylindrical crypt geometry.
 */
#include "CryptSimulation2d.hpp"
/*
 * The next header file defines a Wnt singleton class, which (if used) deals with the
 * imposed Wnt gradient in our crypt model. This affects cell proliferation in the case
 * where we construct each cell with a {{{WntCellCycleModel}}}.
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
class TestRunningCryptSimulationsWithMutationsTutorial : public CxxTest::TestSuite
{
public:
    /* EMPTYLINE
     *
     * == Test 1: a mesh based crypt simulation with mutations ==
     *
     * EMPTYLINE
     *
     * In the first test, we demonstrate how to create a crypt simulation using a
     * cylindrical mesh, with each cell progressing through a simple wnt dependent cell-cycle model,
     * with the Wnt concentration varying within the crypt in a predefined manner.
     */
    void TestMeshBasedCryptWithMutations() throw(Exception)
    {
        /* We first re-initialize time to zero and reseed the random number generator. */
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(0);

        /* We then create a cylindrical mesh, and get the cell location indices, exactly as before. */
        CylindricalHoneycombMeshGenerator generator(6, 9, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        /* We create the cells, using the same method as before. Here, though, we use a {{{SimpleWntCellCycleModel}}}.*/
        std::vector<CellPtr> cells;
        CryptCellsGenerator<SimpleWntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        /* Make pointers to any mutations you want to use.
         * you need to do this before making the cell population or otherwise the numbers of
         * each type of mutation aren't tracked.
         */
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());

        /* We create the cell population, as before. */
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        /* In order to visualise mutant cells and to count how many cells there are of each type you need to use the following command.*/
        cell_population.SetOutputCellMutationStates(true);

        /*
         * We set the height of the crypt. As well as passing this variable into the {{{sloughingCellKiller}}},
         * we will pass it to the {{{WntConcentration}}} object (see below).
         */
        double crypt_height = 8.0;

        /*
         * When using a {{{SiumpleWntCellCycleModel}}}, we need a way of telling each cell what the Wnt concentration
         * is at its location. To do this, we set up a {{{WntConcentration}}} object. Like {{{SimulationTime}}},
         * {{{WntConcentration}}} is a singleton class, so when instantiated it is accessible from anywhere in
         * the code (and in particular, all cells and cell-cycle models can access it). We need to say what
         * the profile of the Wnt concentation should be up the crypt: here, we say it is {{{LINEAR}}} (linear
         * decreasing from 1 to 0 from the bottom of the crypt to the top). We also need to inform the
         * {{{WntConcentration}}} of the cell population and the height of the crypt.
         */
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_height);

        /* Create a simulator as before (except setting a different output directory). */
        CryptSimulation2d simulator(cell_population);
        simulator.SetOutputDirectory("MeshBasedCryptWithMutations");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10);

        /* As before, we create a force law and cell killer and pass these objects to the simulator, then call
         * Solve(). */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_height));
        simulator.AddCellKiller(p_killer);

        simulator.Solve();

        /*
         * Now we have run the simulation to a steady state we select a cell to become mutant.
         *
         *  EMPTYLINE
         *
         * First we re set the end time to some later time.
         */
        simulator.SetEndTime(20);

        /*
         * Now we select one of the cells and set the mutation state to {{{ApcTwoHitCellMutationState}}} (i.e. p_state).
         */
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                                                 cell_iter != cell_population.End();
                                                 ++cell_iter)
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);

            if (node_index == 74) // Chosen from looking at the results from steady state.
            {
                cell_iter->SetMutationState(p_state);
            }
        }

        // Change the Drag on the Mutant cell
        cell_population.SetDampingConstantMutant(10*cell_population.GetDampingConstantNormal());

        /* Solve to new end time */
        simulator.Solve();

        /*
         * Finally, we must tidy up by destroying the {{{SimulationTime}}} and the {{{WntConcentration}}}
         * singleton objects. This avoids memory leaks occurring. */
        WntConcentration<2>::Destroy();
        SimulationTime::Destroy();
    }
    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/MeshBasedCryptWithMutations/results_from_time_0}}}.
     *
     * These are the results before we add the mutations do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/MeshBasedCryptWithMutations/results_from_time_10}}}
     * to see the results from after the mutation has been added.
     *
     * We may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
     * java executable.
     *
     * In the results folderes there is also a file {{{cellmutationstates.dat}}} which tracks the numbers of each mutation type in the simulation.
     * These results are just tab separated columns so may be visualised by using gnuplot, Matlab or similar.
     *
     * EMPTYLINE
     */
};
#endif /*TESTRUNNINGCRYPTSIMULATIONSWITHMUTATIONSTUTORIAL_HPP_*/
