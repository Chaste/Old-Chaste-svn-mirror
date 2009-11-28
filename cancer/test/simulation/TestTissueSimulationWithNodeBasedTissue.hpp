/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef TESTTISSUESIMULATIONWITHNODEBASEDTISSUE_HPP_
#define TESTTISSUESIMULATIONWITHNODEBASEDTISSUE_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "TissueSimulationArchiver.hpp"

#include "TissueSimulation.hpp"
#include "NodeBasedTissue.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "RandomCellKiller.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "LogFile.hpp"


class TestTissueSimulationWithNodeBasedTissue : public AbstractCellBasedTestSuite
{
private:
    template<unsigned DIM>
    std::vector<TissueCell> SetUpCells(TetrahedralMesh<DIM,DIM>* pMesh)
    {
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<pMesh->GetNumNodes(); i++)
        {
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                (TissueConfig::Instance()->GetStemCellG1Duration()
                                    + TissueConfig::Instance()->GetSG2MDuration() );
            TissueCell cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        return cells;
    }

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
     * Create a simulation of a NodeBasedTissue with a NodeBasedTissueMechanicsSystem.
     * Test that no exceptions are thrown, and write the results to file.
     */
    void TestSimpleMonolayer() throw (Exception)
    {
        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each node. Give each a random birth time.
        std::vector<TissueCell> cells = SetUpCells(p_mesh);

        // Create a node based tissue
        NodeBasedTissue<2> node_based_tissue(*p_mesh, cells);

        // Create a mechanics system
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(node_based_tissue, force_collection);
        simulator.SetOutputDirectory("TestTissueSimulationWithNodeBasedTissue");
        simulator.SetEndTime(10.0);

        simulator.Solve();

        // Check that nothing's gone badly wrong by testing that nodes aren't too close together
        double min_distance_between_cells = 1.0;

        for (unsigned i=0; i<simulator.rGetTissue().GetNumNodes(); i++)
        {
            for (unsigned j=i+1; j<simulator.rGetTissue().GetNumNodes(); j++)
            {
                double distance = norm_2(simulator.rGetTissue().GetNode(i)->rGetLocation()-simulator.rGetTissue().GetNode(j)->rGetLocation());
                if (distance < min_distance_between_cells)
                {
                    min_distance_between_cells = distance;
                }
            }
        }

        TS_ASSERT(min_distance_between_cells > 1e-3);
    }

    void TestSimulationWithNodeBoxes() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each node. Give each a random birth time.
        std::vector<TissueCell> cells = SetUpCells(p_mesh);

        // Create a node based tissue
        NodeBasedTissue<2> node_based_tissue(*p_mesh, cells);

        // Create a mechanics system
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(node_based_tissue, force_collection);
        simulator.SetOutputDirectory("TestTissueSimulationWithNodeBasedTissue");
        simulator.SetEndTime(10.0);

        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Check that nothing's gone badly wrong by testing that nodes aren't too close together
        double min_distance_between_cells = 1.0;

        for (unsigned i=0; i<simulator.rGetTissue().GetNumNodes(); i++)
        {
            for (unsigned j=i+1; j<simulator.rGetTissue().GetNumNodes(); j++)
            {
                double distance = norm_2(simulator.rGetTissue().GetNode(i)->rGetLocation()-simulator.rGetTissue().GetNode(j)->rGetLocation());
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
//    void TestSimpleMonolayer2() throw (Exception)
//    {
//        // Create a simple mesh
//        int num_cells_depth = 100;
//        int num_cells_width = 100;
//        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
//        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
//
//        // Set up cells, one for each node. Give each a random birth time.
//        std::vector<TissueCell> cells = SetUpCells(p_mesh);
//
//        // Create a tissue
//        MeshBasedTissue<2> tissue(*p_mesh, cells);
//
//        Meineke2001SpringSystem<2> mechanics_system(tissue);
//
//        // Create a tissue simulation
//        TissueSimulation<2> simulator(tissue, &mechanics_system);
//        simulator.SetOutputDirectory("TestSimpleMonolayer2");
//        simulator.SetEndTime(0.05);
//
//        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
//    }

    /**
     * Create a simulation of a NodeBasedTissue with a NodeBasedTissueMechanicsSystem
     * and a CellKiller. Test that no exceptions are thrown, and write the results to file.
     */
    void TestCellDeath() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each node. Give each a random birth time.
        std::vector<TissueCell> cells = SetUpCells(p_mesh);

        // Create a node based tissue
        NodeBasedTissue<2> node_based_tissue(*p_mesh, cells);

        // Create a mechanics system
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(node_based_tissue, force_collection);
        simulator.SetOutputDirectory("TestTissueSimulationWithNodeBasedTissueCellDeath");
        simulator.SetEndTime(0.5);

        // Add cell killer
        RandomCellKiller<2> random_cell_killer(&node_based_tissue, 0.997877574);
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
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each node. Give each a random birth time.
        std::vector<TissueCell> cells = SetUpCells(p_mesh);

        // Create a node based tissue
        NodeBasedTissue<2> node_based_tissue(*p_mesh, cells);

        // Create a mechanics system
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(node_based_tissue, force_collection);
        simulator.SetOutputDirectory("TestTissueSimulationWithNodeBasedTissueStandardResult");
        simulator.SetEndTime(2.5);

        // Solve
        simulator.Solve();

        // Check some results
        std::vector<double> node_3_location = simulator.GetNodeLocation(3);
        TS_ASSERT_DELTA(node_3_location[0], 2.9415, 1e-4);
        TS_ASSERT_DELTA(node_3_location[1], 0.0136, 1e-4);

        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 3.7813, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], -0.3702, 1e-4);
    }

    // Testing Save
    void TestSave() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        TetrahedralMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each node. Give each a random birth time.
        std::vector<TissueCell> cells = SetUpCells(p_mesh);

        // Create a node based tissue
        NodeBasedTissue<2> node_based_tissue(*p_mesh, cells);

        // Create a mechanics system
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(node_based_tissue, force_collection);
        simulator.SetOutputDirectory("TestTissueSimulationWithNodeBasedTissueSaveAndLoad");
        simulator.SetEndTime(0.1);
        // Solve
        simulator.Solve();

        // Save the results
        TissueSimulationArchiver<2, TissueSimulation<2> >::Save(&simulator);
    }

    // Testing Load (based on previous two tests)
    void TestLoad() throw (Exception)
    {
        // Load the simulation from the TestSave method above and
        // run it from 0.1 to 1.0
        TissueSimulation<2>* p_simulator1;
        p_simulator1 = TissueSimulationArchiver<2, TissueSimulation<2> >::Load("TestTissueSimulationWithNodeBasedTissueSaveAndLoad", 0.1);
        p_simulator1->SetEndTime(1.0);
        p_simulator1->Solve();

        // Save, then reload and run from 1.0 to 2.5

        TissueSimulationArchiver<2, TissueSimulation<2> >::Save(p_simulator1);
        TissueSimulation<2>* p_simulator2
            = TissueSimulationArchiver<2, TissueSimulation<2> >::Load("TestTissueSimulationWithNodeBasedTissueSaveAndLoad", 1.0);

        p_simulator2->SetEndTime(2.5);
        p_simulator2->Solve();
        // These results are from time 2.5 in TestStandardResultForArchivingTestBelow()
        std::vector<double> node_3_location = p_simulator2->GetNodeLocation(3);
        TS_ASSERT_DELTA(node_3_location[0], 2.9415, 1e-4);
        TS_ASSERT_DELTA(node_3_location[1], 0.0136, 1e-4);

        std::vector<double> node_4_location = p_simulator2->GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 3.7813, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], -0.3702, 1e-4);

        // Tidy up
        delete p_simulator1;
        delete p_simulator2;
    }

};

#endif /*TESTTISSUESIMULATIONWITHNODEBASEDTISSUE_HPP_*/

