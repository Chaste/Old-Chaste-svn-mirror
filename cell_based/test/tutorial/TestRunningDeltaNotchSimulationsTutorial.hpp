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
 * = An example showing how to run Delta-Notch simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 *  Still need to summarize the main points of this tutorial!
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
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SmartPointers.hpp"
/*
 * Talk about this header!
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
class TestRunningDeltaNotchSimulationsTutorial : public CxxTest::TestSuite
{
public:

    /*
     * Still need to document this test!
     */
    void TestVertexBasedMonolayerWithDeltaNotch() throw (Exception)
    {
        SimulationTime::Instance()->SetStartTime(0.0);

        /* First we create a regular vertex mesh. */
        HoneycombVertexMeshGenerator generator(5, 5);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        /* We then create some cells, each with a cell-cycle model, {{{DeltaNotchCellCycleModel}}}, which 
         * incorporates a delta-notch ODE system. */
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
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

        /* Using the vertex mesh and cells, we create a cell-based population object. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Configure output
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellAges(true);
        cell_population.SetOutputCellVolumes(true);
        cell_population.SetOutputCellVariables(true);

        // Create and initialize CellwiseData
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(p_mesh->GetNumElements(), 3);
        p_data->SetCellPopulation(&cell_population);
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            p_data->SetValue(RandomNumberGenerator::Instance()->ranf(), cell_population.GetLocationIndexUsingCell(*cell_iter), 0);
            p_data->SetValue(RandomNumberGenerator::Instance()->ranf(), cell_population.GetLocationIndexUsingCell(*cell_iter), 1);
            p_data->SetValue(RandomNumberGenerator::Instance()->ranf(), cell_population.GetLocationIndexUsingCell(*cell_iter), 2);
        }

        // Create and configure cell-based simulation
        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexBasedMonolayerWithDeltaNotch");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(5.0);

        // Create force law and add to simulation
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        /* Finally, as before, we call {{{Destroy()}}} on any singleton classes. */
        CellwiseData<2>::Destroy();
        SimulationTime::Destroy();
    }

    /*
     * Still need to document this test!
     */
    void TestNodeBasedMonolayerWithDeltaNotch() throw (Exception)
    {
        SimulationTime::Instance()->SetStartTime(0.0);

        // Create a 2D honeycomb mesh
        HoneycombMeshGenerator generator(5, 5);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create some cells, each with a cell-cycle model that incorporates a Delta-Notch ODE system
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            p_model->SetDimension(2);
            p_model->SetMaxTransitGenerations(UINT_MAX);
            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);
        cell_population.SetCellAncestorsToLocationIndices();

        // Create and initialize CellwiseData
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

        // Configure output
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellIdData(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAges(true);

        // Create and configure cell-based simulation
        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestNodeBasedMonolayerWithDeltaNotch");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(5.0);
        
        // Set up force law and add to simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        /* Finally, as before, we call {{{Destroy()}}} on any singleton classes. */
        CellwiseData<2>::Destroy();
        SimulationTime::Destroy();

        /* To avoid memory leaks, we also delete any pointers we created in the test. */
        delete p_mesh;
    }
};

#endif /*TESTRUNNINGDELTANOTCHSIMULATIONSTUTORIAL_HPP_*/
