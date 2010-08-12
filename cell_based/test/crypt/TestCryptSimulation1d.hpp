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
#ifndef TESTCRYPTSIMULATION1D_HPP_
#define TESTCRYPTSIMULATION1D_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "TissueSimulationArchiver.hpp"

#include "CryptSimulation1d.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SloughingCellKiller.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "CellPropertyRegistry.hpp"

class TestCryptSimulation1d : public AbstractCellBasedTestSuite
{
public:

    /**
     * In this test, there is no birth as only differentiated cells are used.
     * There is also no death because the we have not passed a cell killer
     * into the simulation. We just perturb one of the nodes and check that
     * each spring relaxes to its rest length.
     */
    void Test1dCryptWithNoBirthOrDeath() throw(Exception)
    {
        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(21);

        ChastePoint<1> shifted_point;
        shifted_point.rGetLocation()[0] = 10.5;
        mesh.SetNode(10, shifted_point);

        // Set up cells
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_model));
            double birth_time = 0.0 - node_index;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create tissue
        MeshBasedTissue<1> crypt(mesh, cells);

        // Create force law
        GeneralisedLinearSpringForce<1> linear_force;
        std::vector<AbstractForce<1>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt1dWithNoBirthAndNoDeath");
        simulator.SetEndTime(1.0);

        // Run simulation
        simulator.Solve();

        // Note: converges very slowly so large tolerance of 0.1
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {
            Node<1>* p_node = mesh.GetNode(index);

            TS_ASSERT_EQUALS(p_node->IsDeleted(), false);

            double coord = p_node->rGetLocation()[0];
            TS_ASSERT_DELTA(coord, index, 0.05);
        }

        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt1dWithNoBirthAndNoDeath", false);

        std::string node_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        TS_ASSERT_EQUALS(system(("diff " + node_results_file + " cell_based/test/data/Crypt1dWithNoBirthAndNoDeath/results.viznodes").c_str()), 0);

        std::string cell_type_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizcelltypes";
        TS_ASSERT_EQUALS(system(("diff " + cell_type_results_file + " cell_based/test/data/Crypt1dWithNoBirthAndNoDeath/results.vizcelltypes").c_str()), 0);

        std::string elem_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizelements";
        TS_ASSERT_EQUALS(system(("diff " + elem_results_file + " cell_based/test/data/Crypt1dWithNoBirthAndNoDeath/results.vizelements").c_str()), 0);
    }


    /**
     * In this test, we pass a sloughing cell killer into the simulation, and
     * check that a cell starting at the end of the crypt is sloughed off.
     */
    void Test1dCryptWithDeathButNoBirth() throw(Exception)
    {
        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(24);

        // Set up cells
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_model));
            double birth_time = 0.0 - node_index;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create tissue
        MeshBasedTissue<1> crypt(mesh, cells);

        // Create force law
        GeneralisedLinearSpringForce<1> linear_force;
        std::vector<AbstractForce<1>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt1dWithDeathButNoBirth");
        simulator.SetEndTime(1.0);

        // Add sloughing cell killer to simulation
        SloughingCellKiller<1> sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        TS_ASSERT_EQUALS(simulator.rGetTissue().GetNumRealCells(), 25u);

        // Run simulation
        simulator.Solve();

        // Since the default crypt length is 22, two cells should have been sloughed
        TS_ASSERT_EQUALS(simulator.rGetTissue().GetNumRealCells(), 23u);

        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt1dWithDeathButNoBirth", false);

        std::string node_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        TS_ASSERT_EQUALS(system(("diff " + node_results_file + " cell_based/test/data/Crypt1dWithDeathButNoBirth/results.viznodes").c_str()), 0);

        std::string cell_type_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizcelltypes";
        TS_ASSERT_EQUALS(system(("diff " + cell_type_results_file + " cell_based/test/data/Crypt1dWithDeathButNoBirth/results.vizcelltypes").c_str()), 0);

        std::string elem_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizelements";
        TS_ASSERT_EQUALS(system(("diff " + elem_results_file + " cell_based/test/data/Crypt1dWithDeathButNoBirth/results.vizelements").c_str()), 0);
    }


    /**
     * In this test, we allow cells to proliferate.
     */
    void Test1dCryptWithBirthButNoDeath() throw (Exception)
    {
        // Get pointers to singleton objects
        RandomNumberGenerator* p_rand_gen = RandomNumberGenerator::Instance();
        TissueConfig* p_params = TissueConfig::Instance();

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(22);

        // Set up cells by iterating through the nodes
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            CellProliferativeType cell_type;
            unsigned generation;
            double birth_time;
            if (i == 0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -p_rand_gen->ranf()*(p_params->GetStemCellG1Duration()
                                                  + p_params->GetSG2MDuration()); // hours - doesn't matter for stem cell
            }
            else if (i < 15)
            {
                cell_type = TRANSIT;
                generation = 1 + (i - 1) / 5;
                birth_time = -p_rand_gen->ranf()*(p_params->GetTransitCellG1Duration()
                                                    + p_params->GetSG2MDuration()); // hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0; // hours
            }

            StochasticDurationGenerationBasedCellCycleModel* p_model = new StochasticDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(cell_type);
            p_model->SetGeneration(generation);

            TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_model));
            p_cell->InitialiseCellCycleModel();
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Create tissue
        MeshBasedTissue<1> crypt(mesh, cells);

        // Create force law
        GeneralisedLinearSpringForce<1> linear_force;
        std::vector<AbstractForce<1>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt1dWithCells");
        simulator.SetEndTime(10.0);

        // Run simulation
        simulator.Solve();

        // There should be 34 cells at the end of the simulation
        TS_ASSERT_EQUALS(simulator.rGetTissue().GetNumRealCells(), 34u);

        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt1dWithCells", false);

        std::string node_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        TS_ASSERT_EQUALS(system(("diff " + node_results_file + " cell_based/test/data/Crypt1dWithCells/results.viznodes").c_str()), 0);

        std::string cell_type_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizcelltypes";
        TS_ASSERT_EQUALS(system(("diff " + cell_type_results_file + " cell_based/test/data/Crypt1dWithCells/results.vizcelltypes").c_str()), 0);

        std::string elem_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizelements";
        TS_ASSERT_EQUALS(system(("diff " + elem_results_file + " cell_based/test/data/Crypt1dWithCells/results.vizelements").c_str()), 0);
    }


    /**
     * In this test, we check that the daughters of a cell that has just divided
     * are put in the correct positions.
     */
    void TestCalculateCellDivisionVector() throw (Exception)
    {
        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(3);

        // Set up cells by iterating through the nodes
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double birth_time;
            if (i == 1)
            {
                birth_time = -23.0;
            }
            else
            {
                birth_time = -1.0;
            }

            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_model));
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Create tissue
        MeshBasedTissue<1> crypt(mesh, cells);

        // Create force law
        GeneralisedLinearSpringForce<1> linear_force;
        std::vector<AbstractForce<1>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("TestCalculateCellDivisionVector");
        simulator.SetEndTime(1.001);

        // Run simulation
        simulator.Solve();

        // There should be 5 cells at the end of the simulation
        TS_ASSERT_EQUALS(simulator.rGetTissue().GetNumRealCells(), 5u);

        // The second cell should have just divided, so be a distance 0.15 away from its initial location
        TS_ASSERT_DELTA(simulator.rGetTissue().GetNode(1)->rGetLocation()[0], 1.15, 1e-3);

        // The new cell should also be a distance 0.15 away from the second cell's initial location
        TS_ASSERT_DELTA(simulator.rGetTissue().GetNode(4)->rGetLocation()[0], 0.85, 1e-3);

        // The other cells should still be at their initial locations
        TS_ASSERT_DELTA(simulator.rGetTissue().GetNode(0)->rGetLocation()[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(simulator.rGetTissue().GetNode(2)->rGetLocation()[0], 2.0, 1e-3);
        TS_ASSERT_DELTA(simulator.rGetTissue().GetNode(3)->rGetLocation()[0], 3.0, 1e-3);
    }


    /**
     * In this test, we include cell birth and cell death.
     */
    void Test1dCryptWithBirthAndDeath() throw (Exception)
    {
        // Get pointers to singleton objects
        RandomNumberGenerator* p_rand_gen = RandomNumberGenerator::Instance();
        TissueConfig* p_params = TissueConfig::Instance();

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(22);

        // Set up cells by iterating through the nodes
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            CellProliferativeType cell_type;
            unsigned generation;
            double birth_time;
            if (i == 0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -p_rand_gen->ranf()*(p_params->GetStemCellG1Duration()
                                                  + p_params->GetSG2MDuration()); // hours - doesn't matter for stem cell
            }
            else if (i < 15)
            {
                cell_type = TRANSIT;
                generation = 1 + (i - 1) / 5;
                birth_time = -p_rand_gen->ranf()*(p_params->GetTransitCellG1Duration()
                                                    + p_params->GetSG2MDuration()); // hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0; // hours
            }

            StochasticDurationGenerationBasedCellCycleModel* p_model = new StochasticDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(cell_type);
            p_model->SetGeneration(generation);

            TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_model));
            p_cell->InitialiseCellCycleModel();
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Create tissue
        MeshBasedTissue<1> crypt(mesh, cells);

        // Create force law
        GeneralisedLinearSpringForce<1> linear_force;
        std::vector<AbstractForce<1>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt1dWithCellsAndGrowth");
        simulator.SetEndTime(10.0);

        // Add sloughing cell killer to simulation
        SloughingCellKiller<1> sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();

        // There should be 34 cells at the end of the simulation
        TS_ASSERT_EQUALS(simulator.rGetTissue().GetNumRealCells(), 30u);

        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt1dWithCellsAndGrowth", false);

        std::string node_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        TS_ASSERT_EQUALS(system(("diff " + node_results_file + " cell_based/test/data/Crypt1dWithCellsAndGrowth/results.viznodes").c_str()), 0);

        std::string cell_type_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizcelltypes";
        TS_ASSERT_EQUALS(system(("diff " + cell_type_results_file + " cell_based/test/data/Crypt1dWithCellsAndGrowth/results.vizcelltypes").c_str()), 0);

        std::string elem_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizelements";
        TS_ASSERT_EQUALS(system(("diff " + elem_results_file + " cell_based/test/data/Crypt1dWithCellsAndGrowth/results.vizelements").c_str()), 0);
    }


    void Test1DChainWithTysonNovakCellsAndNoDeath() throw (Exception)
    {
        // Get pointers to singleton objects
        RandomNumberGenerator* p_rand_gen = RandomNumberGenerator::Instance();
        TissueConfig* p_params = TissueConfig::Instance();

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(22);

        // For Tyson-Novak Cells
        double temp_stem = p_params->GetStemCellG1Duration() + p_params->GetSG2MDuration();
        double temp_transit = p_params->GetTransitCellG1Duration() + p_params->GetSG2MDuration();
        p_params->SetStemCellG1Duration(0.12);
        p_params->SetTransitCellG1Duration(0.12);
        p_params->SetSDuration(0.01);
        p_params->SetG2Duration(0.01);
        p_params->SetMDuration(0.01);

        // Set up cells by iterating through the nodes
        unsigned num_cells_at_start = mesh.GetNumNodes();
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            CellProliferativeType cell_type;
            double birth_time;
            if (i == 0)
            {
                cell_type = STEM;
                birth_time = -p_rand_gen->ranf()*(p_params->GetStemCellG1Duration()
                                                  + p_params->GetSG2MDuration());
            }
            else if (i < 15)
            {
                cell_type = TRANSIT;
                birth_time = -p_rand_gen->ranf()*(p_params->GetTransitCellG1Duration()
                                                    + p_params->GetSG2MDuration());
            }
            else
            {
                cell_type = DIFFERENTIATED;
                birth_time = 0;
            }

            TysonNovakCellCycleModel* p_model = new TysonNovakCellCycleModel();
            p_model->SetCellProliferativeType(cell_type);

            TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_model));
            p_cell->InitialiseCellCycleModel();
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Create tissue
        MeshBasedTissue<1> crypt(mesh, cells);

        // Create force law
        GeneralisedLinearSpringForce<1> linear_force;
        std::vector<AbstractForce<1>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("CryptWithTysonNovakCells");
        simulator.SetEndTime(1.35);

        // Run simulation
        simulator.Solve();

        /*
         * Check that several cell divisions have occurred. The Tyson-Novak cell
         * cycle has period 1.25, so by the end of the simulation the number of
         * cells should have doubled.
         *
         * Note that this test will fail if it is run for much longer, because the
         * Tyson-Novak cell-cycle period is so short that the cells because too
         * squashed together.
         */
        TS_ASSERT_EQUALS(simulator.rGetTissue().GetNumRealCells(), num_cells_at_start + 23u);
        TS_ASSERT_EQUALS(simulator.rGetTissue().GetNumRealCells(), 2*num_cells_at_start);

        p_params->SetStemCellG1Duration(temp_stem - 10.0);
        p_params->SetTransitCellG1Duration(temp_transit - 10.0);
    }


    /**
     * Create a crypt containing a single stem cell and all other cells differentiated.
     * Check that there is the correct number of cells at the end of a simulation, and
     * that they are in the correct order.
     */
    void Test1dChainCorrectCellNumbers()
    {
        // Get pointers to singleton object
        TissueConfig* p_params = TissueConfig::Instance();

        // The stem cell cycle time must still be 24 h, otherwise this test may not pass
        TS_ASSERT_DELTA(p_params->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetSG2MDuration(), 10.0, 1e-12);

        p_params->SetCryptLength(5.0);

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(5);

        // For coverage, shift a node to a negative position
        ChastePoint<1> shifted_point;
        shifted_point.rGetLocation()[0] = -0.5;
        mesh.SetNode(0, shifted_point);

        // Set up cells by iterating through the nodes
        unsigned num_cells = mesh.GetNumNodes();
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        for (unsigned i=0; i<num_cells; i++)
        {
            CellProliferativeType cell_type = DIFFERENTIATED;
            unsigned generation = 4;
            double birth_time = 0;
            if (i == 0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = 0;
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0;
            }

            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(cell_type);
            p_model->SetGeneration(generation);

            TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_model));
            p_cell->InitialiseCellCycleModel();
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Create tissue
        MeshBasedTissue<1> crypt(mesh, cells);
        crypt.SetOutputCellProliferativeTypes(true);

        // Create force law
        GeneralisedLinearSpringForce<1> linear_force;
        std::vector<AbstractForce<1>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt1dTestCorrectCellNumbers");
        simulator.SetEndTime(40);

        // Add sloughing cell killer to simulation
        SloughingCellKiller<1> sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();

        unsigned num_cells_at_end = simulator.rGetTissue().GetNumRealCells();

        /*
         * At the end of the simulation, there should be 6 live cells and one sloughed cell.
         * Note that throughout the simulation 2 cells are sloughed off, but one place is reused.
         */
        TS_ASSERT_EQUALS(num_cells_at_end, 6u);

        // Check we have the correct number of each cell type
        std::vector<unsigned> cell_type_count = simulator.rGetTissue().rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_type_count.size(), 3u);
        TS_ASSERT_EQUALS(cell_type_count[0], 1u);
        TS_ASSERT_EQUALS(cell_type_count[1], 2u);
        TS_ASSERT_EQUALS(cell_type_count[2], 3u);

        for (AbstractTissue<1>::Iterator cell_iter = simulator.rGetTissue().Begin();
             cell_iter != simulator.rGetTissue().End();
             ++cell_iter)
        {
            c_vector<double, 1> cell_location = simulator.rGetTissue().GetLocationOfCellCentre(*cell_iter);
            double x = cell_location[0];

            if (fabs(x) < 1e-2)
            {
                TS_ASSERT_EQUALS(cell_iter->GetCellCycleModel()->GetCellProliferativeType(), STEM);
            }
            else if ( (fabs(x-1) < 1e-2) || (fabs(x-2) < 1e-2) )
            {
                TS_ASSERT_EQUALS(cell_iter->GetCellCycleModel()->GetCellProliferativeType(), TRANSIT);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->GetCellCycleModel()->GetCellProliferativeType(), DIFFERENTIATED);
            }
        }
    }

    /**
     * Test with Wnt-dependent cells.
     */
    void TestWntCellsCannotMoveAcrossYEqualsZero() throw (Exception)
    {
        // Get pointers to singleton object
        TissueConfig* p_params = TissueConfig::Instance();

        // The stem cell cycle time must still be 24 h, otherwise this test may not pass
        TS_ASSERT_DELTA(p_params->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetSG2MDuration(), 10.0, 1e-12);

        p_params->SetCryptLength(5.0);

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(3);

        // Set up cells by iterating through the nodes
        // (don't use any stem cells as we want to test the jiggling)
        unsigned num_cells = mesh.GetNumNodes();
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        for (unsigned i=0; i<num_cells; i++)
        {
            WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
            p_cell_cycle_model1->SetDimension(1);
            p_cell_cycle_model1->SetCellProliferativeType(TRANSIT);
            TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_cell_cycle_model1));
            p_cell->SetBirthTime(0.0);
            cells.push_back(p_cell);
        }

        // Create tissue
        MeshBasedTissue<1> crypt(mesh, cells);

        crypt.SetOutputCellMutationStates(true);
        crypt.SetOutputCellProliferativeTypes(true);

        AbstractTissue<1>::Iterator cell_iterator = crypt.Begin();
        cell_iterator->SetBirthTime(-1.0);   // Make cell cycle models do minimum work
        ++cell_iterator;
        cell_iterator->SetBirthTime(-1.0);

        boost::shared_ptr<AbstractCellProperty> p_apc1(crypt.GetCellPropertyRegistry()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_bcat1(crypt.GetCellPropertyRegistry()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_label(crypt.GetCellPropertyRegistry()->Get<CellLabel>());

        cell_iterator->AddCellProperty(p_label);
        ++cell_iterator;
        cell_iterator->SetBirthTime(-1.0);
        cell_iterator->SetMutationState(p_apc1);
        ++cell_iterator;
        cell_iterator->SetBirthTime(-1.0);
        cell_iterator->SetMutationState(p_bcat1);

        // Create an instance of a Wnt concentration
        WntConcentration<1>::Instance()->SetType(LINEAR);
        WntConcentration<1>::Instance()->SetTissue(crypt);

        // Create force law
        GeneralisedLinearSpringForce<1> linear_force;
        std::vector<AbstractForce<1>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        CryptSimulation1d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt1DWntMatureCells");
        simulator.SetEndTime(0.01);

        // Run simulation
        simulator.Solve();

        // Check that nothing has moved below y=0
        for (AbstractTissue<1>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            TS_ASSERT_LESS_THAN(-1e-15, crypt.GetLocationOfCellCentre(*cell_iter)[0]);
        }

        std::vector<unsigned> cell_mutation_state_count = simulator.rGetTissue().GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_state_count.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_state_count[0], 2u);
        TS_ASSERT_EQUALS(cell_mutation_state_count[1], 1u);
        TS_ASSERT_EQUALS(cell_mutation_state_count[2], 0u); // No APC two hit
        TS_ASSERT_EQUALS(cell_mutation_state_count[3], 1u);

        std::vector<unsigned> cell_type_count = simulator.rGetTissue().rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_type_count.size(), 3u);
        TS_ASSERT_EQUALS(cell_type_count[0], 0u);
        TS_ASSERT_EQUALS(cell_type_count[1], 4u);
        TS_ASSERT_EQUALS(cell_type_count[2], 0u);

        // Tidy up
        WntConcentration<1>::Destroy();
    }


    /**
     * Test saving a CryptSimulation1d object.
     */
    void TestSave() throw (Exception)
    {
        // Get pointers to singleton object
        TissueConfig* p_params = TissueConfig::Instance();
        RandomNumberGenerator* p_rand_gen = RandomNumberGenerator::Instance();

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(22);

        // Set up cells by iterating through the nodes
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            CellProliferativeType cell_type;
            unsigned generation;
            double birth_time;
            if (i == 0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -p_rand_gen->ranf()*(p_params->GetStemCellG1Duration()
                                                  + p_params->GetSG2MDuration());
            }
            else if (i < 15)
            {
                cell_type = TRANSIT;
                generation = 1 + (i - 1)/5;
                birth_time = -p_rand_gen->ranf()*(p_params->GetTransitCellG1Duration()
                                                    + p_params->GetSG2MDuration());
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0;
            }

            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(cell_type);
            p_model->SetGeneration(generation);

            TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_model));
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Create tissue
        MeshBasedTissue<1> crypt(mesh, cells);

        // Create force law
        GeneralisedLinearSpringForce<1> linear_force;
        std::vector<AbstractForce<1>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt1DSaveAndLoad");

        // Our full end time is 0.25, here we run until 0.1 then load and run more below.
        simulator.SetEndTime(0.1);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller<1> sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();

        // Save the results
        TissueSimulationArchiver<1, CryptSimulation1d>::Save(&simulator);
    }


    /**
     * Test loading a CryptSimulation1d object.
     */
    void TestLoad() throw (Exception)
    {
        // Load the simulation from the TestSave method above and
        // run it from 0.1 to 0.2
        CryptSimulation1d* p_simulator1;

        p_simulator1 = TissueSimulationArchiver<1, CryptSimulation1d>::Load("Crypt1DSaveAndLoad", 0.1);

        p_simulator1->SetEndTime(0.2);
        p_simulator1->Solve();

        // Save then reload and run from 0.2 to 0.25
        MutableMesh<1,1>& r_mesh1 = (static_cast<MeshBasedTissue<1>*>(&(p_simulator1->rGetTissue())))->rGetMesh();
        TissueSimulationArchiver<1, CryptSimulation1d>::Save(p_simulator1);

        CryptSimulation1d* p_simulator2 = TissueSimulationArchiver<1, CryptSimulation1d>::Load("Crypt1DSaveAndLoad", 0.2);
        MutableMesh<1,1>& r_mesh2 = (static_cast<MeshBasedTissue<1>*>(&(p_simulator2->rGetTissue())))->rGetMesh();

        TS_ASSERT_EQUALS(r_mesh1.GetNumAllNodes(), r_mesh2.GetNumAllNodes());
        TS_ASSERT_EQUALS(r_mesh1.GetNumNodes(), r_mesh2.GetNumNodes());
        TS_ASSERT_EQUALS(r_mesh1.GetNumBoundaryNodes(), r_mesh2.GetNumBoundaryNodes());

        for (unsigned i=0; i<r_mesh1.GetNumAllNodes(); i++)
        {
            Node<1>* p_node = r_mesh1.GetNode(i);
            Node<1>* p_node2 = r_mesh2.GetNode(i);
            TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
            TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());
            TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());
            TS_ASSERT_DELTA(p_node->rGetLocation()[0], p_node2->rGetLocation()[0], 1e-16);
        }

        TS_ASSERT_EQUALS(r_mesh1.GetNumElements(), r_mesh2.GetNumElements());
        TS_ASSERT_EQUALS(r_mesh1.GetNumAllElements(), r_mesh2.GetNumAllElements());
        TS_ASSERT_EQUALS(r_mesh1.GetNumBoundaryElements(), r_mesh2.GetNumBoundaryElements());
        TS_ASSERT_EQUALS(r_mesh1.GetNumAllBoundaryElements(), r_mesh2.GetNumAllBoundaryElements());

        MutableMesh<1,1>::ElementIterator it = r_mesh1.GetElementIteratorBegin();
        MutableMesh<1,1>::ElementIterator it2 = r_mesh2.GetElementIteratorBegin();
        for (;
             it != r_mesh1.GetElementIteratorEnd();
             ++it, ++it2)
        {
            TS_ASSERT_EQUALS(it->GetNumNodes(), it2->GetNumNodes());
            for (unsigned i=0; i<it->GetNumNodes(); i++)
            {
                TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(i), it2->GetNodeGlobalIndex(i));
            }
        }

        p_simulator2->SetEndTime(0.25);

        // Run simulation
        p_simulator2->Solve();

        std::vector<double> node_4_location = p_simulator2->GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 4.000, 1e-4);

        std::vector<double> node_11_location = p_simulator2->GetNodeLocation(11);
        TS_ASSERT_DELTA(node_11_location[0], 11.000, 1e-4);

        // Tidy up
        delete p_simulator1;
        delete p_simulator2;
    }

};

#endif /*TESTCRYPTSIMULATION1D_HPP_*/

