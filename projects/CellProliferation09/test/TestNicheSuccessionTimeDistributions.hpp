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
#ifndef TESTNICHESUCCESSIONTIMEDISTRIBUTIONS_HPP_
#define TESTNICHESUCCESSIONTIMEDISTRIBUTIONS_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cancer headers
#include "TissueSimulationArchiver.hpp"

#include "CryptSimulation2d.hpp"
#include "MeinekeInteractionWithVariableSpringConstantsForce.hpp"
#include "SloughingCellKiller.hpp"
#include "StochasticCellCycleModelCellsGenerator.hpp"
#include "SimpleWntCellCycleModelCellsGenerator.hpp"
#include "StochasticWntCellCycleModelCellsGenerator.hpp"
#include "IngeWntSwatCellCycleModelCellsGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "LogFile.hpp"
#include "AbstractCancerTestSuite.hpp"


/*
 * Simulations to generate data on the distribution of the time it takes for
 * the 2D crypt to reach a monoclonal population.
 *
 * To run this set of experiments:
 * 1. Run the TestGenerateSteadyStateCrypt.hpp file to generate steady-state archive.
 *
 * 2. Copy the resulting archive folders into the relevant data folder in
 *    projects/CellProliferation09/test/data/<output_directory>/archive
 *
 * 3. Run the test
 *
 * 4. Run the matlab file monoclonality_graphs.m to generate graphs.
 */
class TestNicheSuccessionTimeDistributions : public AbstractCancerTestSuite
{

public:

    /**
     * Here we load a crypt simulation in dynamic equilibrium and label each of the
     * cells as having an ancestor, which is equivalent to their current node position.
     *
     * We then evolve the crypt in a sequence of short runs until it becomes monoclonal
     * (when all cells share a common ancestor).
     *
     * This experiment is repeated many times, in order to build up an estimate for the
     * distribution of niche succession times (the time from labelling to monoclonality)
     * for the Meineke (unwrapped cylinder) crypt geometry.
     *
     */
    void TestGetNicheSuccessionTimeDistribution() throw (Exception)
    {
        // Set experimental parameters
        unsigned num_experiments = 50;          // number of successive experiments to simulate
        double load_time = 300;                 // time at which simulation must be first loaded
        double max_experiment_duration = 6000;  // maximum duration of each experiment
        double time_of_each_run = 10;           // write experiment data to files every 10 hours

        // Change directory_to_copy_from according to which archive (i.e. choice of geometry and cell cycle model) is used
        std::stringstream directory_to_copy_from;
        directory_to_copy_from << "projects/CellProliferation09/test/data/SteadyStateSimpleWnt/sunter3_archive";

        // Set output directory
        std::string output_directory = "NicheSuccessionTime";

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
        out_stream overall_results_file = results_handler.OpenOutputFile("overall_results_2d_.dat");

        // Run experiments
        for (unsigned i=0; i<num_experiments; i++)
        {
            // Set up results files
            std::stringstream ss;
            ss << i;
            out_stream ancestor_results_file = results_handler.OpenOutputFile("ancestor_results_2d_" + ss.str() + ".dat");
            out_stream cell_type_results_file = results_handler.OpenOutputFile("cell_type_results_2d_" + ss.str() + ".dat");

            // Load simulation
            CryptSimulation2d* p_simulator = TissueSimulationArchiver<2, CryptSimulation2d>::Load(output_directory, load_time);

            // Set output directory for simulation
            p_simulator->SetOutputDirectory(output_directory);

            // Set simulation to output cell ancestors too
            p_simulator->SetOutputCellAncestors(true);

            // Label each cell according to its current node index
            p_simulator->rGetTissue().SetCellAncestorsToNodeIndices();

            // Run the experiment for several short runs
            for (double t=load_time; t<max_experiment_duration+load_time+0.5; t+=time_of_each_run)
            {
                // Reset end time and run simulation
                p_simulator->SetEndTime(t + time_of_each_run);
                p_simulator->Solve();

                // Get ancestors
                std::set<unsigned> ancestor_set = p_simulator->rGetTissue().GetCellAncestors();

                // Output time, ancestor indices and ancestor population numbers to ancestor results file

                *ancestor_results_file << SimulationTime::Instance()->GetTime();

                for (std::set<unsigned>::iterator ancestor_iter = ancestor_set.begin();
                     ancestor_iter != ancestor_set.end();
                     ++ancestor_iter)
                {
                    *ancestor_results_file << "\t" << *ancestor_iter;

                    unsigned num_cells_with_this_ancestor = 0;
                    for (AbstractTissue<2>::Iterator cell_iter = p_simulator->rGetTissue().Begin();
                         cell_iter != p_simulator->rGetTissue().End();
                         ++cell_iter)
                    {
                         if (cell_iter->GetAncestor() == *ancestor_iter)
                         {
                            num_cells_with_this_ancestor++;
                         }
                    }
                    *ancestor_results_file << "\t" << num_cells_with_this_ancestor;
                }
                *ancestor_results_file << "\n" << std::flush;

                // Output time and cell type numbers to cell type results file

                *cell_type_results_file << SimulationTime::Instance()->GetTime();

                c_vector<unsigned, NUM_CELL_TYPES> cell_type_count = p_simulator->GetCellTypeCount();

                for (unsigned cell_type_index=0; cell_type_index<NUM_CELL_TYPES-1; cell_type_index++) // ignore necrotic cells
                {
                    *cell_type_results_file << "\t" << cell_type_count[cell_type_index];
                }
                *cell_type_results_file << "\n" << std::flush;

                // If crypt has become monoclonal, then print the current time and stop
                if (ancestor_set.size() == 1u)
                {
                    std::cout << SimulationTime::Instance()->GetTime() - load_time << std::endl;
                    *overall_results_file << SimulationTime::Instance()->GetTime() - load_time << "\n" << std::flush;
                    break;
                }
            }

            // Output trace
            std::cout << "\n Monoclonality experiment finished \n" << std::flush;

            // Save simulation and tidy up
            TissueSimulationArchiver<2, CryptSimulation2d>::Save(p_simulator);
            delete p_simulator;

            // Remove results directories to free up space - comment this out if you want to generate movies
            std::string results_directory_to_remove = OutputFileHandler::GetChasteTestOutputDirectory() + "/" + output_directory + "/results_from_time_*";
            system(("rm -r " + results_directory_to_remove).c_str());

            // Reset the time at which to load the simulation
            load_time = SimulationTime::Instance()->GetTime();

            // Close results files
            ancestor_results_file->close();
            cell_type_results_file->close();
        }

        // Close results files and tidy up
        overall_results_file->close();
        LogFile::Instance()->Close();
        WntConcentration::Destroy();
    }

};

#endif /*TESTNICHESUCCESSIONTIMEDISTRIBUTIONS_HPP_*/
