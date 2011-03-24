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
#ifndef TESTCELLBASEDSIMULATIONWITHNODEBASEDCELLPOPULATION_HPP_
#define TESTCELLBASEDSIMULATIONWITHNODEBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "CellBasedSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "RandomCellKiller.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "LogFile.hpp"
#include "WildTypeCellMutationState.hpp"
#include "PlaneBoundaryCondition.hpp"

#include "Debug.hpp"

class TestCellBasedSimulationWithNodeBasedCellPopulation : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    /**
     * Create a simulation of a NodeBasedCellPopulation with a NodeBasedCellPopulationMechanicsSystem.
     * Test that no exceptions are thrown, and write the results to file.
     */
    void TestSimpleMonolayer() throw (Exception)
    {
        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(*p_mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell-based simulation
        CellBasedSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestCellBasedSimulationWithNodeBasedCellPopulation");
        simulator.SetEndTime(10.0);

        // Create a force law and pass it to the simulation
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        simulator.AddForce(&linear_force);

        simulator.Solve();

        // Check that nothing's gone badly wrong by testing that nodes aren't too close together
        double min_distance_between_cells = 1.0;

        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            for (unsigned j=i+1; j<simulator.rGetCellPopulation().GetNumNodes(); j++)
            {
                double distance = norm_2(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()-simulator.rGetCellPopulation().GetNode(j)->rGetLocation());
                if (distance < min_distance_between_cells)
                {
                    min_distance_between_cells = distance;
                }
            }
        }

        TS_ASSERT(min_distance_between_cells > 1e-3);
    }

    void TestSimulationWithBoxes() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(*p_mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell-based simulation
        CellBasedSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestCellBasedSimulationWithNodeBasedCellPopulation");
        simulator.SetEndTime(10.0);

        // Create a force law and pass it to the simulation
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        simulator.AddForce(&linear_force);

        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Check that nothing's gone badly wrong by testing that nodes aren't too close together
        double min_distance_between_cells = 1.0;

        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            for (unsigned j=i+1; j<simulator.rGetCellPopulation().GetNumNodes(); j++)
            {
                double distance = norm_2(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()-simulator.rGetCellPopulation().GetNode(j)->rGetLocation());
                if (distance < min_distance_between_cells)
                {
                    min_distance_between_cells = distance;
                }
            }
        }

        TS_ASSERT(min_distance_between_cells > 1e-3);
    }

    // results: with a few cells and small end times, Simple was twice as fast as meineke
    //          with 10000 cells, and t_end=0.05, (fixed cell cycle) takes 6.5 mins
    //          => 2 hours real time to do 1hr simulation time
    //   run commented test before to see how meineke does with 10000 cells
    //
    // \TODO this is out of date and needs to be updated with new force law system
//    void TestSimpleMonolayer2() throw (Exception)
//    {
//        // Create a simple mesh
//        int num_cells_depth = 100;
//        int num_cells_width = 100;
//        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
//        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
//
//        // Set up cells, one for each node. Give each a random birth time.
//        std::vector<CellPtr> cells = SetUpCells(p_mesh);
//
//        // Create a cell population
//        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
//
//        Meineke2001SpringSystem<2> mechanics_system(cell_population);
//
//        // Create a cell-based simulation
//        CellBasedSimulation<2> simulator(cell_population, &mechanics_system);
//        simulator.SetOutputDirectory("TestSimpleMonolayer2");
//        simulator.SetEndTime(0.05);
//
//        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
//    }

    /**
     * Create a simulation of a NodeBasedCellPopulation with a NodeBasedCellPopulationMechanicsSystem
     * and a CellKiller. Test that no exceptions are thrown, and write the results to file.
     */
    void TestCellDeath() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        // Create a node based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(*p_mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell-based simulation
        CellBasedSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestCellBasedSimulationWithNodeBasedCellPopulationCellPtrDeath");
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        simulator.AddForce(&linear_force);

        // Add cell killer
        RandomCellKiller<2> random_cell_killer(&node_based_cell_population, 0.997877574);
        simulator.AddCellKiller(&random_cell_killer);

        // Solve
        simulator.Solve();

        // Check some results
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 17u);

        std::vector<double> node_3_location = simulator.GetNodeLocation(3);
        TS_ASSERT_DELTA(node_3_location[0], 3.4931, 1e-4);
        TS_ASSERT_DELTA(node_3_location[1], 1.0062, 1e-4);

        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 1.0840, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 1.7208, 1e-4);
    }

    void TestStandardResultForArchivingTestsBelow() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        // Create a node based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(*p_mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell-based simulation
        CellBasedSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestCellBasedSimulationWithNodeBasedCellPopulationStandardResult");
        simulator.SetEndTime(2.5);

        // Create a force law and pass it to the simulation
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        simulator.AddForce(&linear_force);

        // Create some boundary conditions and pass them to the simulation
        PlaneBoundaryCondition<2> boundary_condition(&node_based_cell_population, zero_vector<double>(2), -1.0*unit_vector<double>(2,1));//y>0
        simulator.AddCellPopulationBoundaryCondition(&boundary_condition);

        // Solve
        simulator.Solve();

        // Check some results
        std::vector<double> node_3_location = simulator.GetNodeLocation(3);
        TS_ASSERT_DELTA(node_3_location[0], 2.9062, 1e-4);
        TS_ASSERT_DELTA(node_3_location[1], 0.0172, 1e-4);

        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 3.8102, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 0.000, 1e-4);
    }

    // Testing Save
    void TestSave() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        // Create a node based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(*p_mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell-based simulation
        CellBasedSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestCellBasedSimulationWithNodeBasedCellPopulationSaveAndLoad");
        simulator.SetEndTime(0.1);

        // Create a force law and pass it to the simulation
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        simulator.AddForce(&linear_force);

        // Create some boundary conditions and pass them to the simulation
        PlaneBoundaryCondition<2> boundary_condition(&node_based_cell_population, zero_vector<double>(2), -1.0*unit_vector<double>(2,1));//y>0
        simulator.AddCellPopulationBoundaryCondition(&boundary_condition);

        // Solve
        simulator.Solve();

        // Save the results
        CellBasedSimulationArchiver<2, CellBasedSimulation<2> >::Save(&simulator);
    }

    // Testing Load (based on previous two tests)
    void TestLoad() throw (Exception)
    {
        // Load the simulation from the TestSave method above and
        // run it from 0.1 to 1.0
        CellBasedSimulation<2>* p_simulator1;
        p_simulator1 = CellBasedSimulationArchiver<2, CellBasedSimulation<2> >::Load("TestCellBasedSimulationWithNodeBasedCellPopulationSaveAndLoad", 0.1);

        // Need to reset the Cut off length as not archived properly see #1496
        (dynamic_cast<NodeBasedCellPopulation<2>*>(&(p_simulator1->rGetCellPopulation())))->SetMechanicsCutOffLength(1.5);

        p_simulator1->SetEndTime(1.0);
        p_simulator1->Solve();


        // Save, then reload and run from 1.0 to 2.5
        CellBasedSimulationArchiver<2, CellBasedSimulation<2> >::Save(p_simulator1);
        CellBasedSimulation<2>* p_simulator2
            = CellBasedSimulationArchiver<2, CellBasedSimulation<2> >::Load("TestCellBasedSimulationWithNodeBasedCellPopulationSaveAndLoad", 1.0);

        // Need to reset the Cut off length as not archived properly see #1496
        (dynamic_cast<NodeBasedCellPopulation<2>*>(&(p_simulator2->rGetCellPopulation())))->SetMechanicsCutOffLength(1.5);

        p_simulator2->SetEndTime(2.5);
        p_simulator2->Solve();


        // These results are from time 2.5 in TestStandardResultForArchivingTestBelow()
        std::vector<double> node_3_location = p_simulator2->GetNodeLocation(3);
        TS_ASSERT_DELTA(node_3_location[0], 2.9062, 1e-4);
        TS_ASSERT_DELTA(node_3_location[1], 0.0172, 1e-4);

        std::vector<double> node_4_location = p_simulator2->GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 3.8102, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 0.000, 1e-4);

        // Tidy up
        delete p_simulator1;
        delete p_simulator2;
    }
};

#endif /*TESTCELLBASEDSIMULATIONWITHNODEBASEDCELLPOPULATION_HPP_*/
