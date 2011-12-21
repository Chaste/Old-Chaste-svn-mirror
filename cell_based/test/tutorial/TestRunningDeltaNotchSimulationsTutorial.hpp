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

#ifndef TESTRUNNINGDELTANOTCHSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGDELTANOTCHSIMULATIONSTUTORIAL_HPP_

/*
 * = An example showing how to run Delta/Notch simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to simulate a growing cell monolayer culture
 * into which a simple model of Delta/Notch signalling is incorporated. This model was developed
 * by Collier et al. ("Pattern formation by lateral inhibition with feedback: a mathematical
 * model of delta-notch intercellular signalling", J. Theor. Biol. 183:429-446) and comprises
 * two ODEs to describe the evolution in concentrations of Delta and Notch in each cell. The ODE
 * for Notch includes a reaction term that depends on the mean Delta concentration among neighbouring
 * cells. Thus in this simulation each cell needs to be able to access information about its
 * neighbours. We use the {{{CellwiseData}}} singleton to facilitate this, and introduce a subclass
 * of {{{OffLatticeSimulation}}} called {{{DeltaNotchOffLatticeSimulation}}} to handle the updating
 * of {{{CellwiseData}}} at each time step as cell neighbours change.
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * As in previous tutorials, we begin by including the necessary header files. We have
 * encountered these files already. Recall that often, either {{{CheckpointArchiveTypes.hpp}}}
 * or {{{CellBasedSimulationArchiver.hpp}}} must be included the first Chaste header.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SmartPointers.hpp"
/*
 * The next header file defines a simple stochastic cell-cycle model that includes the functionality
 * for solving each cell's Delta/Notch signalling ODE system at each time step, using information about neighbouring
 * cells through the {{{CellwiseData}}} singleton. We note that in this simple cell-cycle model, the
 * proliferative status of each cell is unaffected by its Delta/Notch activity; such dependence could
 * easily be introduced given an appropriate model of this coupling.
 */
#include "DeltaNotchCellCycleModel.hpp"
/*
 * The next header file defines the class that simulates the evolution of a {{{CellPopulation}}},
 * specialized to deal with updating of the {{{CellwiseData}}} singleton to deal with Delta-Notch
 * signalling between cells.
 */
#include "DeltaNotchOffLatticeSimulation.hpp"

/* Having included all the necessary header files, we proceed by defining the test class.
 */
class TestRunningDeltaNotchSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:

    /*
     * EMPTYLINE
     *
     * == Test 1: a vertex-based monolayer with Delta/Notch signalling ==
     *
     * EMPTYLINE
     *
     * In the first test, we demonstrate how to simulate a monolayer that incorporates
     * Delta/Notch signalling, using a vertex-based approach.
     */
    void TestVertexBasedMonolayerWithDeltaNotch() throw (Exception)
    {
        /* First we create a regular vertex mesh. */
        HoneycombVertexMeshGenerator generator(5, 5);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        /* We then create some cells, each with a cell-cycle model, {{{DeltaNotchCellCycleModel}}}, which
         * incorporates a Delta/Notch ODE system. In this example we choose to make each cell differentiated,
         * so that no cell division occurs. */
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);

        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
         * output to file. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellAges(true);
        cell_population.SetOutputCellVolumes(true);
        cell_population.SetOutputCellVariables(true);

        /* As we are using the {{{CellwiseData}}} singleton to store the information about each cell required to
         * solve the Delta/Notch ODE system, we must first instantiate this singleton and associate it with the
         * cell population. Note that we set the number of variables to 3. This is because each cell's ODE system
         * comprises two ODEs describing Delta/Notch activity, and an additional 'dummy' ODE with zero reaction term
         * that describes how the mean concentration of Delta among neighbouring cells changes. This latter quantity
         * remains constant for the purposes of solving the ODE system over each time step, but is updated at the
         * end of the time step by a method on {{{DeltaNotchOffLatticeSimulation}}}.
         */
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(p_mesh->GetNumElements(), 3);
        p_data->SetCellPopulation(&cell_population);
        /* We choose to initialise the concentrations to random levels in each cell. */
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            p_data->SetValue(RandomNumberGenerator::Instance()->ranf(), cell_population.GetLocationIndexUsingCell(*cell_iter), 0);
            p_data->SetValue(RandomNumberGenerator::Instance()->ranf(), cell_population.GetLocationIndexUsingCell(*cell_iter), 1);
            p_data->SetValue(RandomNumberGenerator::Instance()->ranf(), cell_population.GetLocationIndexUsingCell(*cell_iter), 2);
        }

        /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
         * and run the simulation. We can make the simulation run for longer to see more paterning by changing the end time. */
        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexBasedMonolayerWithDeltaNotch");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(1.0);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();

        /* Finally, as before, we call {{{Destroy()}}} on any singleton classes. */
        CellwiseData<2>::Destroy();
    }

    /*
     * EMPTYLINE
     *
     * To visualize the results, use Paraview. See the UserTutorials/VisualizingWithParaview tutorial for more information.
     *
     * Load the file {{{/tmp/$USER/testoutput/TestVertexBasedMonolayerWithDeltaNotch/results_from_time_0/results.pvd}}}.
     *
     * EMPTYLINE
     *
     * == Test 2 - a node-based monolayer with Delta/Notch signalling ==
     *
     * EMPTYLINE
     *
     * In the next test we run a similar simulation as before, but this time with node-based
     * 'overlapping spheres' model.
     */
    void TestNodeBasedMonolayerWithDeltaNotch() throw (Exception)
    {
        /*
         * Most of the code in this test is the same as in the previous test,
         * except we now create a 'nodes-only mesh' and {{{NodeBasedCellPopulation}}}.
         */
        HoneycombMeshGenerator generator(5, 5);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh);

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        /* The mechanics cut-off length is also used in this simulation to determine nearest
         * neighbours for the purpose of the Delta/Notch intercellular signalling model.
         */
        cell_population.SetMechanicsCutOffLength(1.5);

        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellIdData(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAges(true);

        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(p_mesh->GetNumNodes(), 3);
        p_data->SetCellPopulation(&cell_population);
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            p_data->SetValue(RandomNumberGenerator::Instance()->ranf(), cell_population.GetLocationIndexUsingCell(*cell_iter), 0);
            p_data->SetValue(RandomNumberGenerator::Instance()->ranf(), cell_population.GetLocationIndexUsingCell(*cell_iter), 1);
            p_data->SetValue(RandomNumberGenerator::Instance()->ranf(), cell_population.GetLocationIndexUsingCell(*cell_iter), 2);
        }

        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestNodeBasedMonolayerWithDeltaNotch");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(5.0);

        /* As we are using a node-based cell population, we use an appropriate force law. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        simulator.Solve();

        /* Finally, as before, we call {{{Destroy()}}} on any singleton classes. */
        CellwiseData<2>::Destroy();

        /* To avoid memory leaks, we also delete any pointers we created in the test. */
        delete p_mesh;
    }
    /*
     * EMPTYLINE
     *
     * To visualize the results, use Paraview. See the UserTutorials/VisualizingWithParaview tutorial for more information
     *
     * Load the file {{{/tmp/$USER/testoutput/TestNodeBasedMonolayerWithDeltaNotch/results_from_time_0/results.pvd}}},
     * add a spherical glyph.
     */
};

#endif /*TESTRUNNINGDELTANOTCHSIMULATIONSTUTORIAL_HPP_*/