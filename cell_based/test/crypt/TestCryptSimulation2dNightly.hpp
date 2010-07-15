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
#ifndef TESTCRYPTSIMULATION2DNIGHTLY_HPP_
#define TESTCRYPTSIMULATION2DNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>

#include "CryptSimulation2d.hpp"
#include "CryptCellsGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "WntConcentration.hpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "NumericFileComparison.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"


#include "WntCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
class TestCryptSimulation2dNightly : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        CellBasedEventHandler::Disable(); // these tests fail with event-handling on
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

///////// NON-PERIODIC TESTS - these test the spring system and cell birth etc. /////////

    /**
     * Provides a reasonable test for the ghost node system...
     */
    void Test2DHoneycombMeshNotPeriodic() throw (Exception)
    {
        // Create mesh
        unsigned num_cells_depth = 11;
        unsigned num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set crypt dimensions
        TissueConfig::Instance()->SetCryptLength(crypt_length);
        TissueConfig::Instance()->SetCryptWidth(crypt_width);

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

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
        TissueConfig::Instance()->SetOutputCellProliferativeTypes(true);
        simulator.SetEndTime(12.0);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller<2> sloughing_cell_killer(&crypt, true);
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

        std::vector<unsigned> cell_type_count = crypt.rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_type_count.size(), 3u);
        TS_ASSERT_EQUALS(cell_type_count[0], 6u);   // Stem
        TS_ASSERT_EQUALS(cell_type_count[1], 21u);  // Transit
        TS_ASSERT_EQUALS(cell_type_count[2], 35u);  // Differentiated
    }

    void TestMonolayer() throw (Exception)
    {
        // Create mesh
        unsigned num_cells_depth = 11;
        unsigned num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        TissueConfig::Instance()->SetCryptLength(crypt_length);
        TissueConfig::Instance()->SetCryptWidth(crypt_width);

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true, -1.0);

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);
        TissueConfig::Instance()->SetOutputVoronoiData(true);

        // Set the first cell to be logged
        crypt.Begin()->SetLogged();

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Monolayer");
        TissueConfig::Instance()->SetOutputCellProliferativeTypes(true);
        simulator.SetEndTime(1);

        // Run simulation
        simulator.Solve();

        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = p_mesh->GetNumNodes();
        TS_ASSERT_EQUALS(number_of_cells, 68u);
        TS_ASSERT_EQUALS(number_of_nodes, 152u);

        std::set<unsigned> ghost_indices = crypt.GetGhostNodeIndices();
        TS_ASSERT_EQUALS(number_of_cells + ghost_indices.size(), number_of_nodes);

        std::vector<unsigned> cell_type_count = crypt.rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_type_count.size(), 3u);
        TS_ASSERT_EQUALS(cell_type_count[0], 0u);   // Stem
        TS_ASSERT_EQUALS(cell_type_count[1], 32u);  // Transit
        TS_ASSERT_EQUALS(cell_type_count[2], 36u);  // Differentiated

        TS_ASSERT_EQUALS(crypt.GetVoronoiTessellation()->GetNumElements(), number_of_nodes);
        TS_ASSERT_EQUALS(crypt.GetVoronoiTessellation()->GetNumNodes(), 273u);
    }

    /**
     * Starting with a small mesh with one stem cell and the rest
     * differentiated, check that the number of cells at the end
     * of the simulation is as expected.
     */
    void Test2DCorrectCellNumbers() throw (Exception)
    {
        // Set up singleton class
        TissueConfig* p_params = TissueConfig::Instance();

        // Check the stem cell cycle time is still 24 hrs, otherwise
        // this test might not pass
        TS_ASSERT_DELTA(p_params->GetStemCellG1Duration(), 14, 1e-12);
        TS_ASSERT_DELTA(p_params->GetTransitCellG1Duration(), 2, 1e-12);
        TS_ASSERT_DELTA(p_params->GetSG2MDuration(), 10, 1e-12);

        // Create mesh
        unsigned num_cells_width = 7;
        unsigned num_cells_depth = 5;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        double crypt_width = num_cells_width - 1.0;
        double crypt_length = num_cells_depth - 1.0;

        // Set crypt dimensions
        TissueConfig::Instance()->SetCryptLength(crypt_length);
        TissueConfig::Instance()->SetCryptWidth(crypt_width);

        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = location_indices.size();
        std::vector<TissueCellPtr> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CellProliferativeType cell_type;
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
            boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetGeneration(generation);
            p_model->SetCellProliferativeType(cell_type);

            TissueCellPtr p_cell(new TissueCell(p_state, p_model));
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
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
        SloughingCellKiller<2> sloughing_cell_killer(&crypt, true);
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
            CellProliferativeType type = cell_iter->GetCellCycleModel()->GetCellProliferativeType();

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
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

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
        SloughingCellKiller<2> sloughing_cell_killer(&crypt, true);
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
        CellBasedEventHandler::Enable();

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<WntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetTissue(crypt);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt2DPeriodicWntNightly");
        simulator.SetEndTime(24.0);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller<2> sloughing_cell_killer(&crypt, true);
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
        WntConcentration<2>::Destroy();

        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();
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
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<WntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            double x = p_mesh->GetNode(i)->GetPoint().rGetLocation()[0];
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            double dist_from_3_6 = sqrt((x-3)*(x-3)+(y-6)*(y-6));

            boost::shared_ptr<AbstractCellMutationState> p_healthy(new WildTypeCellMutationState);
            boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
            if (dist_from_3_6 < 1.1)
            {
                cells[i]->SetMutationState(p_apc2);
            }
            else
            {
                cells[i]->SetMutationState(p_healthy);
            }
        }

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetTissue(crypt);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt2DPeriodicMutant");
        simulator.SetEndTime(12.0);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller<2> sloughing_cell_killer(&crypt, true);
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
            if (cell_iter->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
            {
                number_of_mutant_cells++;
            }
        }

        // Tidy up
        WntConcentration<2>::Destroy();
    }

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
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

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
        RandomCellKiller<2> random_cell_killer(&crypt, 0.700619609);
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
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

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
        SloughingCellKiller<2> sloughing_cell_killer(&crypt, true);
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
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

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
        SloughingCellKiller<2> cell_killer(&crypt);
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
        unsigned num_cells_depth = 11;
        unsigned num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        TissueConfig::Instance()->SetCryptLength(crypt_length);
        TissueConfig::Instance()->SetCryptWidth(crypt_width);

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true, -1.0);

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

    /*
     * This tests that the results files are correct (added because of #1130).
     */
    void TestResultsFileForLongerCryptSimulation() throw(Exception)
    {
        // Set some model parameters
        TissueConfig* p_params = TissueConfig::Instance();
        p_params->SetSDuration(7.4);
        p_params->SetG2Duration(1.4);
        p_params->SetMDuration(0.72);
        p_params->SetTransitCellG1Duration(9.4);
        p_params->SetStemCellG1Duration(9.4);

        // Set output directory
        std::string output_directory = "TestResultsFileForLongerCryptSimulation";

        // Create cylindrical mesh
        HoneycombMeshGenerator generator(16, 19, 0, true, 14.1);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells in mesh
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up each cell with a simple Wnt-based cell cycle model
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<SimpleWntCellCycleModel> cell_generator;
        cell_generator.Generate(cells, p_mesh, location_indices, true);

        // Create crypt
        MeshBasedTissue<2> crypt(*p_mesh, cells);

        // Set up instance of WntConcentration singleton and associate it with crypt
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetTissue(crypt);

        // Set up force law
        GeneralisedLinearSpringForce<2> meineke_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&meineke_force);

        // Create crypt simulation
        CryptSimulation2d simulator(crypt, force_collection);

        // Set where to output simulation results
        simulator.SetOutputDirectory(output_directory);

        // Set simulation to output cell types
        TissueConfig::Instance()->SetOutputCellProliferativeTypes(true);

        // Set length of simulation
        simulator.SetEndTime(20.0);

        // Only save results every tenth time step
        simulator.SetSamplingTimestepMultiple(10);

        // Set up sloughing cell killer and pass in to simulation
        AbstractCellKiller<2>* p_cell_killer = new SloughingCellKiller<2>(&simulator.rGetTissue(), 0.01);
        simulator.AddCellKiller(p_cell_killer);

        // Unusual set-up here (corresponds to the Meineke crypt model parameters)
        p_params->SetDampingConstantNormal(1.0);
        p_params->SetDampingConstantMutant(p_params->GetDampingConstantNormal());
        p_params->SetMeinekeSpringStiffness(30.0);
        simulator.UseJiggledBottomCells();

        // Run simulation
        simulator.Solve();

        // Test that results files are correct
        OutputFileHandler handler(output_directory, false);
        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0";

        NumericFileComparison comp_ele(results_dir + "/results.vizelements", "cell_based/test/data/TestResultsFileForLongerCryptSimulation/results.vizelements");
        TS_ASSERT(comp_ele.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizelements cell_based/test/data/TestResultsFileForLongerCryptSimulation/results.vizelements").c_str()), 0);

        NumericFileComparison comp_nodes(results_dir + "/results.viznodes", "cell_based/test/data/TestResultsFileForLongerCryptSimulation/results.viznodes");
        TS_ASSERT(comp_nodes.CompareFiles(1e-15));

        NumericFileComparison comp_celltypes(results_dir + "/results.vizcelltypes", "cell_based/test/data/TestResultsFileForLongerCryptSimulation/results.vizcelltypes");
        TS_ASSERT(comp_celltypes.CompareFiles(1e-15));

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizsetup cell_based/test/data/TestResultsFileForLongerCryptSimulation/results.vizsetup").c_str()), 0);

        // Tidy up
        delete p_cell_killer;
        WntConcentration<2>::Destroy();
    }
};


#endif /*TESTCRYPTSIMULATION2DNIGHTLY_HPP_*/
