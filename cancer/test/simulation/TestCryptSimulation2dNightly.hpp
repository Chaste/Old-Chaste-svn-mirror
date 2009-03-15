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
#ifndef TESTCRYPTSIMULATION2DNIGHTLY_HPP_
#define TESTCRYPTSIMULATION2DNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>

#include "CryptSimulation2d.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "WntConcentration.hpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "FixedDurationGenerationBasedCellCycleModelCellsGenerator.hpp"
#include "WntCellCycleModelCellsGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCancerTestSuite.hpp"
#include "CancerEventHandler.hpp"
#include "VoronoiTessellation.hpp"
#include "../../global/test/NumericFileComparison.hpp"


class TestCryptSimulation2dNightly : public AbstractCancerTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        CancerEventHandler::Disable(); // these tests fail with event-handling on
        mLastStartTime = std::clock();
        AbstractCancerTestSuite::setUp();
    }

    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCancerTestSuite::tearDown();
    }

public:

///////// NON-PERIODIC TESTS - these test the spring system and cell birth etc. /////////

    /**
     * Test the spring system.
     *
     * The cells in this test are given an intial age of 2.0 so that their
     * springs are at their natural length (i.e. we set birth time=-2.0).
     *
     * The mesh is initially a set of 10 by 10 squares, each square made up
     * of two triangles. The horizontal and vertical edges (springs) are at
     * rest length, the diagonals are longer, so this means the mesh skews
     * to a (sloughed) parallelogram, with each triangle trying to become
     * equilateral.
     *
     * If you want to view the results visually, set the end time to 24.0,
     * and the spring system will resemble a parallelogram. However we keep
     * the simulation time at 1.0 in order to keep the test short.
     */

    void Test2DSpringSystem() throw (Exception)
    {
        // Create mesh
        double crypt_length = 10;
        double crypt_width = 10;
        CancerParameters::Instance()->SetCryptLength(crypt_length);
        CancerParameters::Instance()->SetCryptWidth(crypt_width);

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, mesh, std::vector<unsigned>(), false, 0.0, 3.0, 6.5, 8.0);

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(mesh, cells);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);

        // Destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        simulator.SetEndTime(1.0);

        TS_ASSERT_THROWS_ANYTHING(simulator.Solve());// fails because output directory not set

        // Destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        simulator.SetOutputDirectory("Crypt2DSprings");
        simulator.SetEndTime(1.0);
        simulator.SetUpdateTissueRule(false);
        simulator.SetNoBirth(true);

        // Destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        simulator.Solve();
        std::vector<double> node_0_location = simulator.GetNodeLocation(0);
        TS_ASSERT_DELTA(node_0_location[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(node_0_location[1], 0.0, 1e-12);

        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt2DSprings",false);

        std::string node_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        TS_ASSERT_EQUALS(system(("diff " + node_results_file + " cancer/test/data/Crypt2DSpringsResults/results.viznodes").c_str()), 0);

        std::string cell_type_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizcelltypes";
        TS_ASSERT_EQUALS(system(("diff " + cell_type_results_file + " cancer/test/data/Crypt2DSpringsResults/results.vizcelltypes").c_str()), 0);

        std::string elem_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizelements";
        TS_ASSERT_EQUALS(system(("diff " + elem_results_file + " cancer/test/data/Crypt2DSpringsResults/results.vizelements").c_str()), 0);
    }

    /**
     * Provides a reasonable test for the ghost node system...
     */
    void Test2DHoneycombMeshNotPeriodic() throw (Exception)
    {
        // Create mesh
        int num_cells_depth = 11;
        int num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set crypt dimensions
        CancerParameters::Instance()->SetCryptLength(crypt_length);
        CancerParameters::Instance()->SetCryptWidth(crypt_width);

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);

        TS_ASSERT_EQUALS(cells.size(), location_indices.size());
        unsigned original_number_of_ghosts = p_mesh->GetNumNodes() - cells.size();
        TS_ASSERT_EQUALS(original_number_of_ghosts, 84u);

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt2DHoneycombMesh");
        simulator.SetOutputCellTypes(true);
        simulator.SetEndTime(12.0);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();

        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt2DHoneycombMesh", false);

        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = p_mesh->GetNumNodes();
        TS_ASSERT_EQUALS(number_of_cells, 62u);
        TS_ASSERT_EQUALS(number_of_nodes, 146u);

        std::set<unsigned> ghost_indices = crypt.GetGhostNodeIndices();
        TS_ASSERT_EQUALS(number_of_cells + ghost_indices.size(), number_of_nodes);

        c_vector<unsigned, NUM_CELL_TYPES> cell_type_count = crypt.GetCellTypeCount();
        TS_ASSERT_EQUALS(cell_type_count[0], 6u);   // Stem
        TS_ASSERT_EQUALS(cell_type_count[1], 21u);  // Transit
        TS_ASSERT_EQUALS(cell_type_count[2], 35u);  // Differentiated
        TS_ASSERT_EQUALS(cell_type_count[3], 0u);   // Apoptotic
    }

    void TestMonolayer() throw (Exception)
    {
        // Create mesh
        int num_cells_depth = 11;
        int num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        CancerParameters::Instance()->SetCryptLength(crypt_length);
        CancerParameters::Instance()->SetCryptWidth(crypt_width);

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true, -1.0);

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);
        crypt.SetWriteVoronoiData(true,true);

        // Set the first cell to be logged
        crypt.Begin()->SetLogged();

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Monolayer");
        simulator.SetOutputCellTypes(true);
        simulator.SetEndTime(1);

        // Run simulation
        simulator.Solve();

        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = p_mesh->GetNumNodes();
        TS_ASSERT_EQUALS(number_of_cells, 68u);
        TS_ASSERT_EQUALS(number_of_nodes, 152u);

        std::set<unsigned> ghost_indices = crypt.GetGhostNodeIndices();
        TS_ASSERT_EQUALS(number_of_cells + ghost_indices.size(), number_of_nodes);

        c_vector<unsigned, NUM_CELL_TYPES> cell_type_count = crypt.GetCellTypeCount();
        TS_ASSERT_EQUALS(cell_type_count[0], 0u);   // Stem
        TS_ASSERT_EQUALS(cell_type_count[1], 32u);  // Transit
        TS_ASSERT_EQUALS(cell_type_count[2], 36u);  // Differentiated
        TS_ASSERT_EQUALS(cell_type_count[3], 0u);   // Apoptotic

        VoronoiTessellation<2u>& r_tessellation = crypt.rGetVoronoiTessellation();

        /// \todo Voronoi is worked out for the ghost nodes too - is this necessary (probably not??)
        TS_ASSERT_EQUALS(r_tessellation.GetNumFaces(), number_of_cells + ghost_indices.size());
        TS_ASSERT_EQUALS(r_tessellation.GetNumVertices(), 273u);
    }

    /**
     * Starting with a small mesh with one stem cell and the rest
     * differentiated, check that the number of cells at the end
     * of the simulation is as expected.
     */
    void Test2DCorrectCellNumbers() throw (Exception)
    {
        // Set up singleton class
        CancerParameters* p_params = CancerParameters::Instance();

        // Check the stem cell cycle time is still 24 hrs, otherwise
        // this test might not pass
        TS_ASSERT_DELTA(p_params->GetStemCellG1Duration(), 14, 1e-12);
        TS_ASSERT_DELTA(p_params->GetTransitCellG1Duration(), 2, 1e-12);
        TS_ASSERT_DELTA(p_params->GetSG2MDuration(), 10, 1e-12);

        // Create mesh
        int num_cells_width = 7;
        int num_cells_depth = 5;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        double crypt_width = num_cells_width - 1.0;
        double crypt_length = num_cells_depth - 1.0;

        // Set crypt dimensions
        CancerParameters::Instance()->SetCryptLength(crypt_length);
        CancerParameters::Instance()->SetCryptWidth(crypt_width);

        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = location_indices.size();
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CellType cell_type;
            unsigned generation;
            double birth_time;

            if (location_indices[i]==27) // middle of bottom row of cells
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -1;
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = -1; //hours
            }

            TissueCell cell(cell_type, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            (static_cast<FixedDurationGenerationBasedCellCycleModel*>(cell.GetCellCycleModel()))->SetGeneration(generation);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt2DSpringsCorrectCellNumbers");
        simulator.SetEndTime(40); // hours

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();

        // Now count the number of each type of cell
        unsigned num_stem = 0;
        unsigned num_transit = 0;
        unsigned num_differentiated = 0;

        for (AbstractTissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            CellType type = cell_iter->GetCellType();

            if (type==STEM)
            {
                num_stem++;
            }
            else if (type==TRANSIT)
            {
                num_transit++;
            }
            else if (type==DIFFERENTIATED)
            {
                num_differentiated++;
            }
            else
            {
                // shouldn't get here
                TS_ASSERT(false);
            }
        }

        TS_ASSERT_EQUALS(num_stem, 1u);
        TS_ASSERT_EQUALS(num_transit, 2u);

        TS_ASSERT_LESS_THAN(num_differentiated, 25u);
        TS_ASSERT_LESS_THAN(15u, num_differentiated);
    }

///////// PERIODIC TESTS - These test the system as a whole /////////

    void Test2DPeriodicNightly() throw (Exception)
    {
        // Create mesh
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt2DPeriodicNightly");
        simulator.SetEndTime(12.0);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();

        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();

        TS_ASSERT_EQUALS(number_of_cells, 85u);
        TS_ASSERT_EQUALS(number_of_nodes, 133u);
    }

    void TestCrypt2DPeriodicWntNightly() throw (Exception)
    {
        CancerEventHandler::Enable();

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        WntCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create an instance of a Wnt concentration
        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(crypt);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt2DPeriodicWntNightly");
        simulator.SetEndTime(24.0);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();

        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)

        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
#ifdef CHASTE_CVODE
        // divisions occur marginally earlier with CVODE
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 95u);
        TS_ASSERT_EQUALS(number_of_nodes, 143u);
#else
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 94u);
        TS_ASSERT_EQUALS(number_of_nodes, 142u);
#endif //CHASTE_CVODE

        // Tidy up
        WntConcentration::Destroy();

        CancerEventHandler::Headings();
        CancerEventHandler::Report();
    }

    /**
     * This test is dontTest-ed out and not run every night as it
     * doesn't really test anything. It does show how to set up a
     * mutant simulation. Mutant viscosities are tested elsewhere
     * directly.
     */
    void dontRunTestWithMutantCellsUsingDifferentViscosities() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        WntCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);

        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            CellMutationState mutation_state;

            double x = p_mesh->GetNode(i)->GetPoint().rGetLocation()[0];
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            double dist_from_3_6 = sqrt((x-3)*(x-3)+(y-6)*(y-6));

            if (dist_from_3_6<1.1)
            {
                mutation_state = APC_TWO_HIT;
            }
            else
            {
                mutation_state = HEALTHY;
            }

            cells[i].SetMutationState(mutation_state);
        }

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create an instance of a Wnt concentration
        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(crypt);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt2DPeriodicMutant");
        simulator.SetEndTime(12.0);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();

        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<bool> ghost_cells = crypt.rGetGhostNodes();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();

        TS_ASSERT_EQUALS(number_of_nodes,ghost_cells.size());

        unsigned number_of_cells = 0;
        unsigned number_of_mutant_cells = 0;
        for (AbstractTissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            number_of_cells++;
            if (cell_iter->GetMutationState()==APC_TWO_HIT)
            {
                number_of_mutant_cells++;
            }
        }

        // Tidy up
        WntConcentration::Destroy();
    }

    /// \todo changed test below because once ghost nodes are not
    /// associated with cells, we will never be able to get to
    /// a situation where there are no real cells
    void TestRandomDeathWithPeriodicMesh() throw (Exception)
    {
        unsigned cells_across = 7;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt2DRandomDeathPeriodic");
        simulator.SetEndTime(4.0);

        // Create cell killer and pass in to crypt simulation
        RandomCellKiller<2> random_cell_killer(&crypt, 0.01);
        simulator.AddCellKiller(&random_cell_killer);

        // Run simulation
        simulator.Solve();

        // There should be no cells left after this amount of time
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 1u);
    }

    /**
     * Sloughing with a sloughing cell killer and not turning
     * into ghost nodes on a non-periodic mesh
     */
    void TestSloughingCellKillerOnNonPeriodicCrypt() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt2DSloughingDeathNonPeriodic");
        simulator.SetEndTime(4.0);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();
    }

    void TestSloughingDeathWithPeriodicMesh() throw (Exception)
    {
        unsigned cells_across = 7;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer,true,crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt2DSloughingDeathPeriodic");
        simulator.SetEndTime(4.0);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller cell_killer(&crypt);
        simulator.AddCellKiller(&cell_killer);

        // Run simulation
        simulator.Solve();

        std::vector<bool> ghost_node_indices_after = (static_cast<MeshBasedTissueWithGhostNodes<2>* >(&(simulator.rGetTissue())))->rGetGhostNodes();
        unsigned num_ghosts = 0;
        for (unsigned i=0; i<ghost_node_indices_after.size(); i++)
        {
            if (ghost_node_indices_after[i])
            {
                num_ghosts++;
            }
        }

        // Check no new ghost nodes have been created
        TS_ASSERT_EQUALS(num_ghosts, p_mesh->GetNumNodes() - crypt.GetNumRealCells());

        // There should be this number of cells left after this amount of time
        // (we have lost two rows of 7 but had a bit of birth too)
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 85u);
    }



    void TestMonolayerWithCutoffPointAndNoGhosts() throw (Exception)
    {
        int num_cells_depth = 11;
        int num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        CancerParameters::Instance()->SetCryptLength(crypt_length);
        CancerParameters::Instance()->SetCryptWidth(crypt_width);

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true, -1.0);

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(sqrt(2)); // root2 is a sensible choice
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("MonolayerCutoffPointNoGhosts");
        simulator.SetEndTime(12.0);

        // Run simulation
        simulator.Solve();
    }
};


#endif /*TESTCRYPTSIMULATION2DNIGHTLY_HPP_*/
