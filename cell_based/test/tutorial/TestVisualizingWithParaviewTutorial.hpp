/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef TESTVISUALIZINGWITHPARAVIEWTUTORIAL_HPP_
#define TESTVISUALIZINGWITHPARAVIEWTUTORIAL_HPP_
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
 * van Leeuwen ''et al'' (2009) [doi:10.1111/j.1365-2184.2009.00627.x].
 *
 * The first thing to do is include the following header, which allows us
 * to use certain methods in our test (this header file should be included
 * in any Chaste test).
 */
#include <cxxtest/TestSuite.h>

/*
 * The next header file defines a stochastic cell-cycle model.
 */
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
/* The next header file defines a helper class for generating a suitable mesh. */
#include "HoneycombMeshGenerator.hpp"
/* The next header file defines a {{{CellPopulation}}} class that uses a mesh. */
#include "MeshBasedCellPopulation.hpp"
/* The next header file defines a force law, based on a linear spring, for describing
 * the mechanical interactions between neighbouring cells in the cell population.
 */
#include "GeneralisedLinearSpringForce.hpp"
/* The next header file defines the class that simulates the evolution of a {{{CellPopulation}}}.
 */
#include "CellBasedSimulation.hpp"

/* Next, we define the test class, which inherits from {{{CxxTest::TestSuite}}}
 * and defines some test methods.
 */
class TestVisualizingWithParaviewTutorial : public CxxTest::TestSuite
{
public:
    /* EMPTYLINE
     *
     * == Test 1 - a cell-based monolayer simulation ==
     *
     * EMPTYLINE
     *
     * In the first test, we run a simple cell-based simulation, in which we use
     * a honeycomb mesh, give each cell a stochastic cell-cycle model, and enforce
     * random cell killing.
     */
	void Test2DMonolayerSimulationForVisualizing() throw (Exception)
    {
        /* As in '''all''' cell-based simulations, we must first set the start time.
         * In addition, it is advisable to reset the values of all model parameters.
         * {{{SimulationTime}}} and {{{CellBasedConfig}}} are ''singleton'' classes; this
         * means that one and only one of each of these objects is instantiated at
         * any time, and that that single object is accessible from anywhere in the
         * code. As a result, we do not need to keep passing round the current time or
         * model parameter values.
         */
        SimulationTime::Instance()->SetStartTime(0.0);

        /* Next, we generate a mesh we use the {{{HoneycombMeshGenerator}}}. This
         * generates a honeycomb-shaped mesh, in which all nodes are equidistant.
         * Here the first and second arguments define the size of the mesh - we have
         * chosen a mesh that is 5 nodes (i.e. cells) wide, and 5 nodes high.
         * The third argument indicates that we want 0 ghost nodes around the mesh.
         * The last boolean parameter indicates that we do not want cylindrical boundary
         * conditions.
         */
        HoneycombMeshGenerator generator(5, 5, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetCircularMesh(2.5);

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * Then we loop over the number of nodes in the mesh and assign a cell
         * to each node. Each cell will have a randomly chosen birth time. */
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            StochasticDurationGenerationBasedCellCycleModel* p_model = new StochasticDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(TRANSIT);
            p_model->SetMaxTransitGenerations(UINT_MAX);

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                (CellBasedConfig::Instance()->GetStemCellG1Duration()
                                    + CellBasedConfig::Instance()->GetSG2MDuration() );

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Now we have a mesh and a set of cells to go with it we can create a ''CellPopulation''.
         * In general, this class associates a collection of cells with a set of nodes or a mesh.
         * For this test we use a particular type of cell population called a
         * {{{MeshBasedCellPopulation}}}.
         */
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* In order to output the .vtu files required for paraview, we explicitly
         * instruct the simulation to output the data we need.
         */
        cell_population.SetOutputVoronoiData(true);

        /* We must now create one or more force laws, which determine the mechanics of
         * the cell population. For this test, we assume that a cell experiences a force from each
         * neighbour that can be represented as a linear overdamped spring. We put a pointer
         * to this force into a vector. We use a cut-off point which represents that cells farther
         * than 1.5 cell lengths apart, do not exert forces on one another.
         */
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        /* Now we define the cell-based simulation object, passing in the cell population and collection
         * of force laws: */
        CellBasedSimulation<2> simulator(cell_population, force_collection);

        /* Set the output directory on the simulator (relative to
         * "/tmp/<USER_NAME>/testoutput") and the end time (in hours).
         */
        simulator.SetOutputDirectory("Test2DMonolayerSimulationForVisualizing");
        simulator.SetEndTime(1.0);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* {{{SimulationTime::Destroy()}}} and {{{RandomNumberGenerator::Destroy()}}} '''must''' be called at the end of the test.
         * If not, when {{{SimulationTime::Instance()->SetStartTime(0.0);}}} is called
         * at the beginning of the next test in this file, an assertion will be triggered.
         */
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

};


#endif /* TESTVISUALIZINGWITHPARAVIEWTUTORIAL_HPP_ */
