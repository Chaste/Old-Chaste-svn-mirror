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
#ifndef TESTGENERATEPLOTSOFASINGLECELL_HPP_
#define TESTGENERATEPLOTSOFASINGLECELL_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cancer headers
#include "TissueSimulationArchiver.hpp"

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
 * Simulations to track the progeny of a single cell in the crypt over several experiments.
 *
 * To run this set of experiments:
 * 1. Run the TestGenerateSteadyStateCrypt.hpp file to generate steady-state archive.
 *
 * 2. Copy the resulting archive folders into the relevant data folder in
 *    projects/CellProliferation09/test/data/<output_directory>/archive
 *
 * 3. Run the test
 *
 * 4. (Edit and) Run single_cell_compile_results script
 */
class TestGeneratePlotsOfASingleCell : public AbstractCancerTestSuite
{
public:

    void TestLoadSteadyStateAndFollowASingleCellFromBottomToTop() throw (Exception)
    {
        // Set experimental parameters
        unsigned num_experiments = 40;  // number of successive experiments to simulate
        double load_time = 300;         // time at which simulation must be first loaded
        double time_of_each_run = 10;   // write experiment data to files every 10 hours
        bool first_run = true;          // helper flag

        // Change directory_to_copy_from according to which archive (i.e. choice of geometry and cell cycle model) is used
        std::stringstream directory_to_copy_from;
        directory_to_copy_from << "projects/CellProliferation09/test/data/SteadyStateSimpleWnt/sunter3_archive";

        // Set output directory
        std::string output_directory = "CryptTrackSingleCell";

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

        // Here we do all the setup necessary
        if (first_run)
        {
            // Load and set up simulation
            CryptSimulation2d* p_simulator = TissueSimulationArchiver<2, CryptSimulation2d>::Load(output_directory, load_time);

            // Set output directory for simulation
            p_simulator->SetOutputDirectory(output_directory);

            // Set simulation to output cell mutation states too
            p_simulator->SetOutputCellMutationStates(true);

            // Write results to file for every time step
            p_simulator->SetSamplingTimestepMultiple(1);

            // Follow a cell at the base of the crypt
            MeshBasedTissue<2>* p_crypt = static_cast<MeshBasedTissue<2>*>(&p_simulator->rGetTissue());
            p_crypt->SetWriteVoronoiData(true, true);

            // Set the cell to be logged and label it so that we can see it in visualizer.
            // Note that we should not just pick an index and use rGetCellUsingLocationIndex(),
            // since this may correspond to a ghost node. A safer way of picking a cell is to
            // use the tissue iterator, as below.
            unsigned some_number = 71;
            AbstractTissue<2>::Iterator cell_iter = p_simulator->rGetTissue().Begin();

            for (unsigned i=0; i<some_number; i++)
            {
                ++cell_iter;
            }
            cell_iter->SetLogged();
            cell_iter->SetMutationState(LABELLED);

            *overall_results_file << "Finished setup and altering start archive\n" << std::flush;

            // Save simulation and tidy up
            TissueSimulationArchiver<2, CryptSimulation2d>::Save(p_simulator);
            delete p_simulator;
        }

        // Run the experiment for several short runs
        double t = load_time;
        for (unsigned i=0; i<num_experiments; i++)
        {
            // Load simulation
            CryptSimulation2d* p_simulator = TissueSimulationArchiver<2, CryptSimulation2d>::Load(output_directory, t);

            // Reset end time and run simulation
            t += time_of_each_run;
            p_simulator->SetEndTime(t);
            *overall_results_file << "Experiment " << i << ", time = " << SimulationTime::Instance()->GetTime() << " running..." << std::flush;
            p_simulator->Solve();

            c_vector<unsigned, 5> mutations = p_simulator->GetCellMutationStateCount();

            // Finish this time loop
            *overall_results_file << "... to " << SimulationTime::Instance()->GetTime() << " DONE\n" << std::flush;

            if (mutations[1]==0)
            {
                std::cout << "Healthy cells = " << mutations[0] << ", labelled cells = " << mutations[1] << "\n" << std::flush;
                *overall_results_file << "Labelled population is lost from crypt\n" << std::flush;
                break;
            }

            // Save simulation and tidy up
            TissueSimulationArchiver<2, CryptSimulation2d>::Save(p_simulator);
            delete p_simulator;
        }

        // Close results files and tidy up
        *overall_results_file << "Closing file\n" << std::flush;
        overall_results_file->close();
        LogFile::Instance()->Close();
        WntConcentration::Destroy();
    }
};

#endif /*TESTGENERATEPLOTSOFASINGLECELL_HPP_*/
