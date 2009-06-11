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
#ifndef TESTMEINEKELABELLINGEXPERIMENTSSUNTERDATA_HPP_
#define TESTMEINEKELABELLINGEXPERIMENTSSUNTERDATA_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cancer headers
#include "TissueSimulationArchiver.hpp"

#include "CryptStatistics.hpp"
#include "CryptSimulation2d.hpp"
#include "MeshBasedTissue.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "SloughingCellKiller.hpp"
#include "StochasticCellCycleModelCellsGenerator.hpp"
#include "SimpleWntCellCycleModelCellsGenerator.hpp"
#include "StochasticWntCellCycleModelCellsGenerator.hpp"
#include "IngeWntSwatCellCycleModelCellsGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "LogFile.hpp"
#include "AbstractCancerTestSuite.hpp"


/*
 * Simulations of crypt labelling and sections experiments for a 2D cylindrical crypt.
 *
 * To run this set of experiments:
 * 1. Run the TestGenerateSteadyStateCrypt.hpp file to generate steady-state archive.
 *
 * 2. Copy the resulting archive folders into the relevant data folder in
 *    projects/CellProliferation09/test/data/<output_directory>/archive
 *
 * 3. Run the test
 *
 * 4. (Edit and) Run the script CompileLabellingResults_40_9.sh
 *
 * 5. run plot_40_9_data.m
 */
class TestMeinekeLabellingExperimentsSunterData : public AbstractCancerTestSuite
{
private:

    // Store information from this section in a big vector
    void LogResults(std::vector<bool> labelledThisTime, std::vector<unsigned>& rLabelledCellsCounter)
    {
        for (unsigned cell_index=0; cell_index<labelledThisTime.size(); cell_index++)
        {
            if (cell_index>=rLabelledCellsCounter.size())
            {
                EXCEPTION("labelled_cells_counter vector is not big enough\n");
            }

            if (labelledThisTime[cell_index])
            {
                rLabelledCellsCounter[cell_index]++;
            }
        }
    }

public:

    void TestLoadSteadyStateResultsAndRunSimulations_40_9() throw (Exception)
    {
        // Set experimental parameters
        unsigned num_experiments = 50;          // number of successive experiments to simulate
        double load_time = 300;                 // time at which simulation must be first loaded
        double time_of_each_run = 10;           // write experiment data to files every 10 hours
        double randomise_time = 50.0;           // time to run simulation after loading prior to starting experiments

        // Change the number below to re-seed the random number generator in different simulations
        unsigned experiment_run = 1;

        // Change directory_to_copy_from according to which archive (i.e. choice of geometry and cell cycle model) is used
        std::stringstream directory_to_copy_from;
        directory_to_copy_from << "projects/CellProliferation09/test/data/SteadyStateSimpleWnt/sunter3_archive";

        // Set output directory
        std::string output_directory = "MeinekeLabellingExperiment";

        // Where the archive will be copied to
        std::string directory_to_copy_to = OutputFileHandler::GetChasteTestOutputDirectory() + output_directory;

        // Set up log file
        LogFile::Instance()->Set(1, output_directory, "log.dat");

        // Copy archive to output directory and rename so that it is recognised by TissueSimulationArchiver
        // (comment out the following six lines for running on a cluster and copy the archives there manually)
        std::string command = "cp -rf --remove-destination " + directory_to_copy_from.str() + " " + directory_to_copy_to + "/";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);

        command = "mkdir -p " + directory_to_copy_to + "/archive";
        return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);

        command = "mv " + directory_to_copy_to + "/sunter3_archive/* " + directory_to_copy_to + "/archive";
        return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);

        // Create results file handler
        OutputFileHandler results_handler(output_directory, false);

        // Create overall results file
        out_stream overall_results_file = results_handler.OpenOutputFile("overall_results.dat");
        *overall_results_file << "Starting experiments\n" << std::flush;

        // Calculate when the set of experiments will start and end
        double experiments_start_time = load_time + randomise_time;
        double experiments_end_time = experiments_start_time + num_experiments*time_of_each_run;
        if (experiments_end_time < load_time + randomise_time + time_of_each_run)
        {
            EXCEPTION("End time of simulation must be increased, crypt was not allowed to settle to a random state\n");
        }

        // Load simulation
        CryptSimulation2d* p_simulator = TissueSimulationArchiver<2, CryptSimulation2d>::Load(output_directory, load_time);

        // Set output directory for simulation
        p_simulator->SetOutputDirectory(output_directory);

        // Run simulation for a given duration using the new random number,
        // to get a different starting point for each set of experiments
        RandomNumberGenerator::Instance()->Reseed(experiment_run);
        p_simulator->SetEndTime(load_time + randomise_time);
        p_simulator->Solve();

        // Save simulation and tidy up
        TissueSimulationArchiver<2, CryptSimulation2d>::Save(p_simulator);
        delete p_simulator;

        *overall_results_file << "Finished a " << randomise_time << "hr run to randomise start archive\n" << std::flush;

        // Set up some data structures to record results
        unsigned max_length_of_crypt_section = 100;
        std::vector<unsigned> labelled_cells_counter_t40(max_length_of_crypt_section);
        std::vector<unsigned> labelled_cells_counter_t9(max_length_of_crypt_section);

        for (unsigned i=0; i<max_length_of_crypt_section; i++)
        {
            labelled_cells_counter_t40[i] = 0u;
            labelled_cells_counter_t9[i] = 0u;
        }
        unsigned experiment_counter = 0u;

        // Run experiments
        for (double t=experiments_start_time; t<experiments_end_time; t+=time_of_each_run)
        {
            experiment_counter++;

            // Load simulation
            CryptSimulation2d* p_simulator = TissueSimulationArchiver<2, CryptSimulation2d>::Load(output_directory, t);

            *overall_results_file << "Experiment " << experiment_counter << ", time = " << SimulationTime::Instance()->GetTime() << " running..." << std::flush;

            // Label cells in S phase
            MeshBasedTissue<2>* p_crypt = static_cast<MeshBasedTissue<2>*>(&p_simulator->rGetTissue());
            CryptStatistics crypt_statistics(*p_crypt);
            crypt_statistics.LabelSPhaseCells();

            // Run for 40 minutes
            p_simulator->SetEndTime(t + 0.66666666);  // run for just short of target time (9 hours)
            p_simulator->Solve();

            // Take two t'=0.66666666 (40 minute) sections
            std::vector<TissueCell*> crypt_section = crypt_statistics.GetCryptSection();
            std::vector<bool> labelled_t40_a = crypt_statistics.GetWhetherCryptSectionCellsAreLabelled(crypt_section);
            crypt_section = crypt_statistics.GetCryptSection();
            std::vector<bool> labelled_t40_b = crypt_statistics.GetWhetherCryptSectionCellsAreLabelled(crypt_section);
            LogResults(labelled_t40_a, labelled_cells_counter_t40);
            LogResults(labelled_t40_b, labelled_cells_counter_t40);

            // Run for 9 hours
            p_simulator->SetEndTime(t + time_of_each_run - 1);  // run for just short of target time (9 hours)
            p_simulator->Solve();

            // Take two t'=9 sections
            crypt_section = crypt_statistics.GetCryptSection();
            std::vector<bool> labelled_t9_a = crypt_statistics.GetWhetherCryptSectionCellsAreLabelled(crypt_section);
            crypt_section = crypt_statistics.GetCryptSection();
            std::vector<bool> labelled_t9_b = crypt_statistics.GetWhetherCryptSectionCellsAreLabelled(crypt_section);
            LogResults(labelled_t9_a, labelled_cells_counter_t9);
            LogResults(labelled_t9_b, labelled_cells_counter_t9);

            // Un-label all cells
            crypt_statistics.LabelAllCellsAsHealthy();

            // Run for another hour
            p_simulator->SetEndTime(t + time_of_each_run);
            p_simulator->Solve();

            // Save simulation and tidy up
            TissueSimulationArchiver<2, CryptSimulation2d>::Save(p_simulator);
            delete p_simulator;

            // Write results to file
            std::stringstream string_stream1;
            string_stream1 << "labelled_t40_totals" << experiment_counter << "_run_" << experiment_run << ".dat";
            OutputFileHandler results_handler(output_directory + "/results/", false);
            out_stream file_t40 = results_handler.OpenOutputFile(string_stream1.str());
            for (unsigned i=0; i<labelled_cells_counter_t40.size(); i++)
            {
                *file_t40 << labelled_cells_counter_t40[i] << std::endl;
            }
            file_t40->close();

            std::stringstream string_stream2;
            string_stream2 << "labelled_t9_totals" << experiment_counter << "_run_" << experiment_run << ".dat";
            out_stream file_t9 = results_handler.OpenOutputFile(string_stream2.str());
            for (unsigned i=0 ; i< labelled_cells_counter_t9.size() ; i++)
            {
                *file_t9 << labelled_cells_counter_t9[i] << std::endl;
            }
            file_t9->close();

            // Finish this time loop
            *overall_results_file << "... to " << SimulationTime::Instance()->GetTime() << " DONE \n" << std::flush;
        }

        // Close results files and tidy up
        *overall_results_file << "Closing file\n" << std::flush;
        overall_results_file->close();
        LogFile::Close();
        WntConcentration::Destroy();
    }
};

#endif /*TESTMEINEKELABELLINGEXPERIMENTSSUNTERDATA_HPP_*/
