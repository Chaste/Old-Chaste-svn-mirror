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
#ifndef TESTRUNNINGTUMOURSPHEROIDSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGTUMOURSPHEROIDSIMULATIONSTUTORIAL_HPP_
/*
 * = An example showing how to run tumour spheroid simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste is used to run discrete tumour
 * spheroid simulations. Like crypt simulations, tumour spheroid simulations
 * include cell cycle models and force laws to determine how cells divide and
 * move. In tumour spheroid simulations, however, these are also coupled to a
 * system of partial differential equations that determine the concentration
 * of specified nutrients (e.g. oxygen) throughout the cell population. Also, unlike
 * in crypt simulation, the cell population grows substantially as the cell-based simulation
 * progresses.
 *
 * In summary, the main differences between this tutorial and the crypt simulation
 * tutorials are
 *
 *  * a PDE is defined, to be used in the simulation, and
 *  * a non-periodic mesh is used.
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * The first thing to do is include the following header, which allows us
 * to use certain methods in our test (this header file should be included
 * in any Chaste test):
 */
#include <cxxtest/TestSuite.h>
/*
 * This header file defines a helper class for generating a suitable mesh:
 */
#include "HoneycombMeshGenerator.hpp"
/*
 * These are the classes that will be used in these tests (note that we use a
 * cell-based simulation subclass called {{{CellBasedSimulationWithPdes}}}):
 */
#include "CellBasedSimulationWithPdes.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OxygenBasedCellKiller.hpp"
#include "CellwiseSourcePde.hpp"
#include "WildTypeCellMutationState.hpp"
/*
 * !PetscSetupAndFinalize.hpp must be included in all tests which use Petsc. This is
 * a suite of data structures and routines that are used in the finite element
 * PDE solvers, which is how we solve the nutrient PDE(s).
 */
#include "PetscSetupAndFinalize.hpp"

/*
 * Next, we define the test class, which inherits from {{{CxxTest::TestSuite}}}.
 */
class TestRunningTumourSpheroidSimulationsTutorial : public CxxTest::TestSuite
{
public:
    void TestSpheroidTutorial() throw(Exception)
    {
        /*
         * This first line can be ignored, it's a macro which just says
         * don't run this test if in parallel.
         */
        EXIT_IF_PARALLEL; // defined in PetscTools.hpp

        /*
         * The first thing to do, as before, is to set up the start time and
         * reset the parameters.
         */
        SimulationTime::Instance()->SetStartTime(0.0);
        CellBasedConfig::Instance()->Reset();


        /*
         * Now we want to create a ''non-periodic'' 'honeycomb' mesh.
         * We use the honeycomb mesh generator, as before, saying 10 cells wide
         * and 10 cells high. Note that the thickness of the ghost nodes layer is
         * 0, i.e. no ghost nodes, and the {{{false}}} indicates not cylindrical.
         */
        HoneycombMeshGenerator generator(10, 10, 0, false);
        /*
         * Get the mesh. Note we call {{{GetMesh()}}} rather than {{{GetCyclindricalMesh}}},
         * and that a {{{MutableMesh}}} is returned.
         */
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /*
         * Next, we need to create some cells. Unlike before, we don't just use
         * a {{{CellsGenerator}}} class, but do it manually, in a loop. First,
         * define the cells vector.
         */
        std::vector<CellPtr> cells;

        /*
         * This line defines a mutation state to be used for all cells, of type
         * `WildTypeCellMutationState` (i.e. 'healthy'):
         */
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

        /*
         * Now loop over the nodes...
         */
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            /*
             * ...then create a cell, giving it a particular cell cycle model
             * - {{{SimpleOxygenBasedCellCycleModel}}}. The spatial dimension (1, 2 or 3) and
             * cell proliferative type (STEM, TRANSIT or DIFFERENTIATED) needs to be
             * set on the cell cycle model before being passed to the cell.
             */
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);
            CellPtr p_cell(new Cell(p_state, p_model));

            /*
			 * ...then alter the default cell cycle times
			 */
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            /*
             * We now define a random birth time, chosen from [-T,0], where
             * T = t,,1,, + t,,2,,, where t,,1,, is a parameter representing the G,,1,, duration
             * of a 'stem' cell, and t,,2,, is the basic S+G,,2,,+M phases duration.
             */
            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                 (  p_model->GetStemCellG1Duration()
                                  + p_model->GetSG2MDuration() );
            /*
             * ...then we set the birth time and push the cell back into the vector
             * of cells.
             */
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /*
         * Now that we have defined the cells, we can define the CellPopulation. This time it
         * is just a mesh-based cell population (i.e. not a {{{MeshBasedCellPopulationWithGhostNodes()}}}.
         * Again, the constructor takes in the mesh and the cells vector.
         */
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /*
         * Recall that in the Wnt based crypt simulation, we defined a singleton class
         * which cell-cycles used to get the wnt concentration. Here, we do the same kind
         * of thing, but using the singletom {{{CellwiseData}}} class, which stores the
         * value of the current nutrient concentration, for each cell. We have to
         * tell the {{{CellwiseData}}} object how many cells and variables per cell there
         * are (in this case, 1 variable per cell, i.e. the oxygen concentration), and
         * the cell population.
         */
        CellwiseData<2>::Instance()->SetNumCellsAndVars(cell_population.GetNumRealCells(),1);
        CellwiseData<2>::Instance()->SetCellPopulation(&cell_population);
        /*
         * Then we have to initialise the oxygen concentration for each node (to 1.0), by
         * calling {{{SetValue}}}. This takes in the concentration, and the location index
         * corresponding to the cell which this concentration is for.
         */
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            CellwiseData<2>::Instance()->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
        }

        /*
         * Next we instantiate an instance of the PDE class which we defined above.
         * This will be passed into the simulator. The !CellwiseSourcePde is
         * a Pde class which inherits from !AbstractLinearEllipticPde, and represents
         * the PDE: u_xx + u_yy = k(x) u, where k(x) = -0.03 (the coefficient below)
         * if x is in a live cell, and k(x)=0 if x is within a apoptotic cell
         */
        CellwiseSourcePde<2> pde(cell_population, -0.03);

        /*
         * To pass the PDE to our simulator, it first needs to be encapsulated in a
         * {{{PdeAndBoundaryConditions}}} object, together with the boundary condition for
         * the PDE. The latter is specified by the second and third arguments of the
         * {{{PdeAndBoundaryConditions}}} constructor below: the second argument defines the value
         * of the boundary condition and the third argument defines whether it is of Neumann type
         * (true) or Dirichlet type (false). Thus, in our case, we are a specifying no-flux
         * boundary condition. Note that we currently cannot impose more than one boundary
         * condition for each PDE (so that e.g. we cannot impose a zero-flux boundary condition
         * on some part of the boundary and a fixed-value boundary condition on the rest).
         */
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, 0.0, true);

        /*
         * After having created a {{{PdeAndBoundaryConditions}}} object, we then pass it
         * into a vector of pointers. This allows us to define any number of PDEs within
         * the cell-based simulation, in a similar way to how we create a vector of force laws
         * (see below).
         */
        std::vector<PdeAndBoundaryConditions<2>*> pde_and_bc_collection;
        pde_and_bc_collection.push_back(&pde_and_bc);

        /*
         * We must now create one or more force laws, which determine the mechanics of
         * the cell population. For this test, we assume that a cell experiences a force from each
         * neighbour that can be represented as a linear overdamped spring. Since this
         * model was first proposed in the context of crypt modelling by Meineke ''et al''
         * (Cell Prolif. 34:253-266, 2001), we call this object a
         * {{{GeneralisedLinearSpringForce}}}. We pass a pointer to this force into a vector.
         * Note that we have called the method {{{SetCutOffLength}}} on the
         * {{{GeneralisedLinearSpringForce}}} before passing it into the collection of force
         * laws - this modifies the force law so that two neighbouring cells do not impose
         * a force on each other if they are located more than 3 units (=3 cell widths)
         * away from each other. This modification is necessary when no ghost nodes are used,
         * for example to avoid artificially large forces between cells that lie close together
         * on the spheroid boundary.
         */
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(3);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        /*
         * The simulator object for these problems is
         * {{{CellBasedSimulationWithPdes}}}. We pass in the cell_population, the
         * mechanics system, and the PDE.
         */
        CellBasedSimulationWithPdes<2> simulator(cell_population, force_collection, pde_and_bc_collection);

        /*
         * As with {{{CryptSimulation2d}}} (which inherits from the same base class
         * as {{{CellBasedSimulationWithPdes}}}), we can set the output directory
         * and end time.
         */
        simulator.SetOutputDirectory("SpheroidTutorial");
        simulator.SetEndTime(10.0);

        /*
         * Solve.
         */
        simulator.Solve();

        /*
         * Finally, call {{{Destroy()}}} on the singleton classes. The results
         * can be visualised as in the previous test.
         */
        SimulationTime::Destroy();
        CellwiseData<2>::Destroy();
    }
};
#endif /*TESTRUNNINGTUMOURSPHEROIDSIMULATIONSTUTORIAL_HPP_*/
