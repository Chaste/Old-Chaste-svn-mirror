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
 * In this tutorial we show how Chaste can be used to simulate a growing cell monolayer culture or
 * multicellular tumour spheroid. Like the crypt simulations, tumour spheroid simulations
 * include cell-cycle models and force laws to determine how cells divide and
 * move. In tumour spheroid simulations, however, these are also coupled to a
 * system of partial differential equations (PDEs) that determine the concentration
 * of specified nutrients (e.g. oxygen) throughout the cell population. Also, unlike
 * in a crypt simulation, the cell population may grow substantially as the simulation
 * progresses.
 *
 * In summary, the main differences between this tutorial and the crypt simulation
 * tutorials are: a PDE is defined, to be used in the simulation; and a non-periodic mesh is used.
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * As in the crypt simulation tutorial, we begin by including the necessary header files. We have
 * encountered some of these files already. Recall that often {{{CheckpointArchiveTypes.hpp}}}
 * or {{{CellBasedSimulationArchiver.hpp}}} must be included the first Chaste header.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "RandomNumberGenerator.hpp"
#include "SmartPointers.hpp"
/*
 * The {{{SimpleOxygenBasedCellCycleModel}}} header file defines a cell-cycle model in which
 * a cell's rate of progress through G1 phase changes over time in a simple manner, according
 * to the local oxygen concentration. We also include the {{{WildTypeCellMutationState}}}
 * header file, which defines a wild type cell mutation state that we will use to construct
 * cells. A cell mutation state is always required when constructing a cell, however
 * in the crypt simulation tutorial we used a helper class ({{{CryptCellsGenerator}}}) that
 * allowed us to avoid having to construct cells directly.
 */
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
/*
 * The next three header files define: a PDE that describes how oxygen is transported via through the
 * domain via diffusion and is consumed by live cells; a constant-valued boundary condition to
 * associate with the PDE; and a PDE handler class, which is passed to the simulation object and
 * handles the numerical solution of any PDEs.
 */
#include "CellwiseSourcePde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "CellBasedPdeHandler.hpp"
/*
 * We also include a header file defining a cell killer, which implements the process of
 * hypoxia (low oxygen)-induced cell death.
 */
#include "OxygenBasedCellKiller.hpp"
/*
 * We use an {{{OffLatticeSimulation}}}.
 */
#include "OffLatticeSimulation.hpp"
/*
 * The header file {{{PetscSetupAndFinalize.hpp}}} must be included in all tests which use Petsc. This is
 * a suite of data structures and routines that are used in the finite element
 * PDE solvers, which is how we solve the oxygen transport PDE.
 */
#include "PetscSetupAndFinalize.hpp"

/*
 * Having included all the necessary header files, we proceed by defining the test class.
 */
class TestRunningTumourSpheroidSimulationsTutorial : public CxxTest::TestSuite
{
public:
    void TestSpheroidTutorial() throw(Exception)
    {
        /*
         * This first line can be ignored: it is a macro which just says
         * don't run this test if in parallel.
         */
        EXIT_IF_PARALLEL;

        /*
         * The first thing to do, as before, is to set up the start time.
         */
        SimulationTime::Instance()->SetStartTime(0.0);

        /*
         * Now we want to create a '''non-periodic''' 'honeycomb' mesh.
         * We use the honeycomb mesh generator, as before, saying 10 cells wide
         * and 10 cells high. Note that the thickness of the ghost nodes layer is
         * 0, i.e. there are no ghost nodes, and the {{{false}}} indicates that the
         * returned mesh is '''not''' cylindrical. In contrast to the crypt simulation
         * tutorial, here we call {{{GetMesh()}}} on the {{{HoneycombMeshGenerator}}}
         * object to return the mesh, which is of type {{{MutableMesh}}}.
         */
        HoneycombMeshGenerator generator(10, 10, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /*
         * Next, we need to create some cells. Unlike in the the crypt simulation
         * tutorial, we don't just use a {{{CellsGenerator}}} class, but do it manually,
         * in a loop. First, we define a {{{std::vector}}} of cell pointers.
         */
        std::vector<CellPtr> cells;

        /*
         * This line defines a mutation state to be used for all cells, of type
         * `WildTypeCellMutationState` (i.e. 'healthy'):
         */
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

        /*
         * Now we loop over the nodes...
         */
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            /*
             * ...then create a cell, giving it a {{{SimpleOxygenBasedCellCycleModel}}}.
             * The spatial dimension (1, 2 or 3) and
             * cell proliferative type (STEM, TRANSIT or DIFFERENTIATED) needs to be
             * set on the cell-cycle model before it is passed to the cell.
             */
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);
            CellPtr p_cell(new Cell(p_state, p_model));

            /*
             * We also alter the default cell-cycle times.
             */
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            /*
             * We now define a random birth time, chosen from [-T,0], where
             * T = t,,1,, + t,,2,,, where t,,1,, is a parameter representing the G,,1,, duration
             * of a 'stem' cell, and t,,2,, is the basic S+G,,2,,+M phases duration...
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
         * Now that we have defined the cells, we can define the {{{CellPopulation}}}. We use a
         * {{{MeshBasedCellPopulation}}} since although the cell population is mesh-based, it does
         * not include any ghost nodes. The constructor takes in the mesh and the cells vector.
         */
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /*
         * Recall that in the Wnt-based crypt simulation tutorial, we defined a singleton class
         * which cell-cycles used to get the Wnt concentration. Here, we do something similar
         * using the {{{CellwiseData}}} singleton class, which stores the
         * value of the current nutrient concentration for each cell. We have to
         * tell the {{{CellwiseData}}} object how many cells and variables per cell there
         * are (in this case, one variable per cell, namely the oxygen concentration), and
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
         * This will be passed into the {{{OffLatticeSimulationWithPdes}}} object. The
         * {{{CellwiseSourcePde}}} is a {{{PDE}}} class which inherits from
         * {{{AbstractLinearEllipticPde}}} and represents
         * the PDE ''u_xx'' + ''u_yy'' = ''k''(''x'',''y'') ''u'', where ''u''(''x'',''y'') denotes
         * the oxygen concentration at
         * position (''x'',''y'') and the function ''k''(''x'',''y'') specifies the rate of consumption by live cells
         * there. Here ''k''(''x'',''y'') takes the value -0.03 (the coefficient below) if
         * the cell located at (''x'',''y'') is a live cell, and zero if the cell has died due
         * to oxygen deprivation.
         */
        CellwiseSourcePde<2> pde(cell_population, -0.03);

        /*
         * We also create a constant-valued boundary condition to associate with the PDE.
         * This boundary condition object takes in a single argument in its constructor,
         * the value at the boundary. We also introduce a boolean to specify whether this value is the flux at the boundary
         * (a Neumann boundary condition) or the value of the state variable at the boundary
         * (a Dirichlet boundary condition) below.
         */
        ConstBoundaryCondition<2> bc(1.0);
        bool is_neumann_bc = true;

        /*
         * To pass the PDE to our simulator, it must first be encapsulated in a
         * {{{PdeAndBoundaryConditions}}} object, together with the boundary condition for
         * the PDE. The latter is specified by the second and third arguments of the
         * {{{PdeAndBoundaryConditions}}} constructor below: the second argument defines the value
         * of the boundary condition and the third argument defines whether it is of Neumann type
         * (true) or Dirichlet type (false). Thus, in our case, we are a specifying no-flux
         * boundary condition. Note that we currently cannot impose more than one boundary
         * condition for each PDE (so that e.g. we cannot impose a zero-flux boundary condition
         * on some part of the boundary and a fixed-value boundary condition on the rest), although
         * the boundary condition itself can be made spatially varying or time-dependent.
         */
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, is_neumann_bc);

        /*
         * After having created a {{{PdeAndBoundaryConditions}}} object, we then pass it
         * to a cell-based PDE handler object. This allows us to define any number of PDEs within
         * the cell-based simulation.
         */
        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);

        /*
         * We are now in a position to construct a {{{OffLatticeSimulationWithPdes}}} object,
         * using the cell population. We then pass the PDE handler object to the simulation.
         */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetCellBasedPdeHandler(&pde_handler);

        /*
         * We next set the output directory and end time.
         */
        simulator.SetOutputDirectory("SpheroidTutorial");
        simulator.SetEndTime(1.0);

        /*
         * We must now create one or more force laws, which determine the mechanics of
         * the cell population. As in the crypt simulation tutorial, we assume that a cell
         * experiences a force from each neighbour that can be represented as a linear overdamped
         * spring, so we use a {{{GeneralisedLinearSpringForce}}} object.
         * Note that we have called the method {{{SetCutOffLength}}} on the
         * {{{GeneralisedLinearSpringForce}}} before passing it to the simulator: this call
         * modifies the force law so that two neighbouring cells do not impose
         * a force on each other if they are located more than 3 units (=3 cell widths)
         * away from each other. This modification is necessary when no ghost nodes are used,
         * for example to avoid artificially large forces between cells that lie close together
         * on the spheroid boundary.
         */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(3);
        simulator.AddForce(p_linear_force);

        /*
         * We call {{{Solve()}}} on the simulator to run the simulation.
         */
        simulator.Solve();

        /*
         * Finally, we call {{{Destroy()}}} on the singleton classes. The results
         * can be visualized as in the crypt simulation tutorial.
         */
        SimulationTime::Destroy();
        CellwiseData<2>::Destroy();
    }
};
#endif /*TESTRUNNINGTUMOURSPHEROIDSIMULATIONSTUTORIAL_HPP_*/