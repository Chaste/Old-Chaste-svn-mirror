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
#ifndef TESTCRYPTSTATISTICS_HPP_
#define TESTCRYPTSTATISTICS_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "TissueSimulationArchiver.hpp"

#include "CryptStatistics.hpp"
#include "CryptSimulation2d.hpp"
#include "CryptCellsGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"

/**
 * Note that all these tests call setUp() and tearDown() before running,
 * so if you copy them into a new test suite be sure to copy these methods
 * too.
 */
class TestCryptStatistics : public CxxTest::TestSuite
{
private:

    void setUp()
    {
        // Initialise singleton classes
        SimulationTime::Instance()->SetStartTime(0.0);
        TissueConfig::Instance()->Reset();
    }
    void tearDown()
    {
        // Clear up singleton classes
        SimulationTime::Destroy();
    }

public:

    void TestGetSection() throw (Exception)
    {
        // Create mesh
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 0;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);// true = mature cells

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);
        crypt.InitialiseCells(); // must be called explicitly as there is no simulation

        CryptStatistics crypt_statistics(crypt);

        std::vector<TissueCellPtr> test_section = crypt_statistics.GetCryptSection(0.5,1.5,sqrt(3));

        // Test the cells are correct
        TS_ASSERT_EQUALS(test_section.size(), 6u);

        unsigned expected_indices[6] = {0,1,3,4,7,8};

        for (unsigned i=0; i<test_section.size(); i++)
        {
            TS_ASSERT_EQUALS(crypt.GetLocationIndexUsingCell(test_section[i]), expected_indices[i]);
        }

        // Test that we get a valid section when the x-values are the same
        std::vector<TissueCellPtr> test_section_vertical = crypt_statistics.GetCryptSection(0.5,0.5,sqrt(3));

        // Test the cells are correct
        TS_ASSERT_EQUALS(test_section_vertical.size(), 5u);

        unsigned expected_indices_vertical[6] = {0,1,3,6,7};

        for (unsigned i=0; i<test_section_vertical.size(); i++)
        {
            TS_ASSERT_EQUALS(crypt.GetLocationIndexUsingCell(test_section_vertical[i]), expected_indices_vertical[i]);
        }

        std::vector<TissueCellPtr> test_section_periodic = crypt_statistics.GetCryptSectionPeriodic(0.5,2.5,sqrt(3));

        // Test the cells are correct
        TS_ASSERT_EQUALS(test_section_periodic.size(), 6u);

        unsigned expected_indices_periodic[6] = {0,1,3,5,6,8};

        for (unsigned i=0; i<test_section_periodic.size(); i++)
        {
            TS_ASSERT_EQUALS(crypt.GetLocationIndexUsingCell(test_section_periodic[i]), expected_indices_periodic[i]);
        }

        std::vector<TissueCellPtr> test_section_periodic_2 = crypt_statistics.GetCryptSectionPeriodic(2.5,0.5,sqrt(3));

        // Test the cells are correct
        TS_ASSERT_EQUALS(test_section_periodic_2.size(), 6u);

        unsigned expected_indices_periodic_2[6] = {0,2,3,5,6,7};

        for (unsigned i=0; i<test_section_periodic_2.size(); i++)
        {
            TS_ASSERT_EQUALS(crypt.GetLocationIndexUsingCell(test_section_periodic_2[i]), expected_indices_periodic_2[i]);
        }

        // Test an overwritten method
        std::vector<TissueCellPtr> test_section_periodic_3 = crypt_statistics.GetCryptSectionPeriodic();

        // Test the cells are correct
        TS_ASSERT_EQUALS(test_section_periodic_3.size(), 3u);
        unsigned expected_indices_periodic_3[6] = {2,4,8};

        for (unsigned i=0; i<test_section_periodic_3.size(); i++)
        {
            TS_ASSERT_EQUALS(crypt.GetLocationIndexUsingCell(test_section_periodic_3[i]), expected_indices_periodic_3[i]);
        }
    }


    void TestMakeMeinekeGraphs() throw (Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();

        std::string output_directory = "MakeMeinekeGraphs";

        // Create mesh
        unsigned cells_across = 13;
        unsigned cells_up = 25;
        double crypt_width = 12.1;
        unsigned thickness_of_ghost_layer = 3;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCellPtr> temp_cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(temp_cells, p_mesh, std::vector<unsigned>(), true, 0.3, 2.0, 3.0, 4.0, true);

        // This awkward way of setting up the cells is a result of #430
        std::vector<TissueCellPtr> cells;
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            cells.push_back(temp_cells[location_indices[i]]);
        }

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        CryptSimulation2d simulator(crypt, force_collection, false, false);

        simulator.SetOutputDirectory(output_directory);
        double time_of_each_run = simulator.GetDt(); // for each run

        // Set simulation to output cell types
        TissueConfig::Instance()->SetOutputCellMutationStates(true);

        // Set length of simulation here
        simulator.SetEndTime(time_of_each_run);
        SloughingCellKiller<2> cell_killer(&simulator.rGetTissue(),0.01);
        simulator.AddCellKiller(&cell_killer);

        // UNUSUAL SET UP HERE /////////////////////////////////////

        p_params->SetDampingConstantNormal(1.0);    // normally 1

        // Do not give mutant cells any different movement properties to normal ones
        p_params->SetDampingConstantMutant(p_params->GetDampingConstantNormal());

        p_params->SetMeinekeSpringStiffness(30.0); //normally 15.0;
        // 0.3/30 = 0.01 (i.e. Meineke's values)

        simulator.UseJiggledBottomCells();

        // END OF UNUSUAL SET UP! //////////////////////////////////


        // TEST CryptStatistics::GetCryptSectionPeriodic by labelling a column of cells...
        CryptStatistics crypt_statistics(crypt);
        std::vector<TissueCellPtr> test_section = crypt_statistics.GetCryptSectionPeriodic(8.0,8.0);

        boost::shared_ptr<AbstractCellProperty> p_label(new CellLabel);
        for (unsigned i=0; i<test_section.size(); i++)
        {
            test_section[i]->AddCellProperty(p_label);
        }

        simulator.Solve();
        TissueSimulationArchiver<2, CryptSimulation2d>::Save(&simulator);

        // ... and checking visualization of labelled cells against previous run
        OutputFileHandler handler(output_directory, false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        TS_ASSERT_EQUALS(system(("diff " + results_file + " cell_based/test/data/MakeMeinekeGraphs/results.viznodes").c_str()), 0);

        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizcelltypes";
        TS_ASSERT_EQUALS(system(("diff " + results_file2 + " cell_based/test/data/MakeMeinekeGraphs/results.vizcelltypes").c_str()), 0);

        // TEST crypt_statistics::LabelSPhaseCells

        // First remove labels
        for (AbstractTissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            cell_iter->RemoveCellProperty<CellLabel>();
        }

        crypt_statistics.LabelSPhaseCells();

        // Iterate over cells checking for correct labels
        unsigned counter = 0;
        for (AbstractTissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            bool is_labelled = cell_iter->HasCellProperty<CellLabel>();
            bool in_s_phase = (cell_iter->GetCellCycleModel()->GetCurrentCellCyclePhase() == S_PHASE);

            TS_ASSERT_EQUALS(is_labelled, in_s_phase);

            if (in_s_phase)
            {
                counter++;
            }
        }

        TS_ASSERT_EQUALS(counter, 15u);

        // Test that LabelAllCellsAsHealthy sets all cells back to be UNLABELLED wild type cells
        crypt_statistics.LabelAllCellsAsHealthy();

        // Iterate over cells checking for correct labels
        counter = 0;
        for (AbstractTissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>(), true);
            TS_ASSERT_EQUALS(cell_iter->HasCellProperty<CellLabel>(), false);
            counter++;
        }

        TS_ASSERT_EQUALS(counter, simulator.rGetTissue().GetNumRealCells());

        crypt_statistics.LabelSPhaseCells();

        simulator.SetEndTime(2*time_of_each_run);
        simulator.Solve();

        // TEST CryptStatistics::AreCryptSectionCellsLabelled

        // Set cells which are not in the crypt section to be in state APC +/-, so that we can
        // see the section
        test_section = crypt_statistics.GetCryptSectionPeriodic(8.0,8.0);

        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);

        for (AbstractTissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            bool in_section = false;
            for (unsigned vector_index=0; vector_index<test_section.size(); vector_index++)
            {
                if (test_section[vector_index] == *cell_iter)
                {
                    in_section = true;
                }
            }
            if (!in_section)
            {
                cell_iter->SetMutationState(p_apc1);
            }
        }
        simulator.SetEndTime(3*time_of_each_run);
        simulator.Solve();

        std::vector<TissueCellPtr> crypt_section = crypt_statistics.GetCryptSection(8.0,8.0);
        std::vector<bool> labelled = crypt_statistics.AreCryptSectionCellsLabelled(crypt_section);

        // Test that the vector of booleans corresponds with a visualisation of the data -
        // only the third cell had been labelled
        for (unsigned vector_index=0; vector_index<labelled.size(); vector_index++)
        {
            if (vector_index == 2u)
            {
                TS_ASSERT_EQUALS(labelled[vector_index], true);
            }
            else
            {
                TS_ASSERT_EQUALS(labelled[vector_index], false);
            }
        }
        RandomNumberGenerator::Destroy();
    }


    /**
     * This test runs multiple crypt simulations and records whether
     * or not labelled cells are in a randomly chosen crypt section.
     */
    void TestMultipleCryptSimulations() throw (Exception)
    {
        std::string output_directory = "MakeMoreMeinekeGraphs";

        // Create mesh
        unsigned cells_across = 13;
        unsigned cells_up = 25;
        double crypt_width = 12.1;
        unsigned thickness_of_ghost_layer = 3;

        unsigned num_simulations = 2;

        // Guess of maximum number of cells a crypt section may contain
        unsigned max_length_of_crypt_section = 5 * (unsigned)sqrt(pow(cells_across/2.0+1,2) + pow(cells_up,2));

        std::vector<unsigned> labelled_cells_counter(max_length_of_crypt_section);

        for (unsigned i=0; i<max_length_of_crypt_section; i++)
        {
            labelled_cells_counter[i] = 0;
        }

        TissueConfig* p_params = TissueConfig::Instance();

        p_params->SetDampingConstantNormal(1.0);    // normally 1
        // Do not give mutant cells any different movement properties to normal ones
        p_params->SetDampingConstantMutant(p_params->GetDampingConstantNormal());
        p_params->SetMeinekeSpringStiffness(30.0); //normally 15.0;

        double time_of_each_run;
        AbstractCellKiller<2>* p_cell_killer;
        std::vector<bool> labelled;

        CryptStatistics* p_crypt_statistics;

        // Create tissue
        MeshBasedTissueWithGhostNodes<2>* p_crypt;

        HoneycombMeshGenerator generator = HoneycombMeshGenerator(cells_across, cells_up, thickness_of_ghost_layer, true, crypt_width/cells_across);
        std::vector<unsigned> location_indices;

        Cylindrical2dMesh* p_mesh;
        SimulationTime* p_simulation_time;

        // Loop over the number of simulations
        for (unsigned simulation_index=0; simulation_index<num_simulations; simulation_index++)
        {
            // Create new structures for each simulation
            p_mesh = generator.GetCylindricalMesh();
            location_indices = generator.GetCellLocationIndices();

            // Reset start time
            SimulationTime::Destroy();
            p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);

            // Set up cells
            std::vector<TissueCellPtr> temp_cells;
            CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
            cells_generator.Generate(temp_cells, p_mesh, std::vector<unsigned>(), true, 0.3, 2.0, 3.0, 4.0, true);

            // This awkward way of setting up the cells is a result of #430
            std::vector<TissueCellPtr> cells;
            for (unsigned i=0; i<location_indices.size(); i++)
            {
                cells.push_back(temp_cells[location_indices[i]]);
            }

            // Set up crypt
            p_crypt = new MeshBasedTissueWithGhostNodes<2>(*p_mesh, cells, location_indices);

            // Set up force law
            GeneralisedLinearSpringForce<2> linear_force;
            std::vector<AbstractForce<2>*> force_collection;
            force_collection.push_back(&linear_force);

            // Set up crypt simulation
            CryptSimulation2d simulator(*p_crypt, force_collection, false, false);
            simulator.SetOutputDirectory(output_directory);

            // Set simulation to output cell types
            TissueConfig::Instance()->SetOutputCellMutationStates(true);

            // Set length of simulation here
            time_of_each_run = 10.0*simulator.GetDt(); // for each run
            simulator.SetEndTime(time_of_each_run);

            // Set up cell killer
            p_cell_killer = new SloughingCellKiller<2>(&simulator.rGetTissue(), 0.01);
            simulator.AddCellKiller(p_cell_killer);

            simulator.UseJiggledBottomCells();

            // set up crypt statistics
            p_crypt_statistics = new CryptStatistics(*p_crypt);

            // run for a bit
            simulator.Solve();
            p_crypt_statistics->LabelSPhaseCells();

            simulator.SetEndTime(2.0*time_of_each_run);
            simulator.Solve();

            std::vector<TissueCellPtr> crypt_section = p_crypt_statistics->GetCryptSection(8.0, 8.0);
            labelled = p_crypt_statistics->AreCryptSectionCellsLabelled(crypt_section);

            // Store information from this simulation in a global vector
            for (unsigned cell_index=0; cell_index < labelled.size(); cell_index++)
            {
                TS_ASSERT(cell_index < labelled_cells_counter.size());

                if (labelled[cell_index])
                {
                    labelled_cells_counter[cell_index]++;
                }
            }

            // Tidy up
            cells.clear();
            labelled.clear();
            WntConcentration<2>::Destroy();

            delete p_crypt_statistics;
            delete p_crypt;
            delete p_cell_killer;
        }

        // Calculate percentage of labelled cells at each position in 'labelled_cells_counter'
        std::vector<double> percentage_of_labelled_cells(max_length_of_crypt_section);
        for (unsigned index=0; index < max_length_of_crypt_section; index ++)
        {
            percentage_of_labelled_cells[index] = (double) labelled_cells_counter[index]/num_simulations;
        }

        // Write data to file
        SimpleDataWriter writer1(output_directory, "percentage_of_labelled_cells.dat", percentage_of_labelled_cells, false);

        // Test against previous run
        // ... and checking visualization of labelled cells against previous run
        OutputFileHandler handler(output_directory, false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "percentage_of_labelled_cells.dat";
        TS_ASSERT_EQUALS(system(("diff " + results_file + " cell_based/test/data/MakeMoreMeinekeGraphs/percentage_of_labelled_cells.dat").c_str()), 0);

        delete p_params;
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTCRYPTSTATISTICS_HPP_*/
