/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#ifndef TESTWILDTYPEMONOCLONALITY_HPP_
#define TESTWILDTYPEMONOCLONALITY_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CryptSimulation2dWithAncestorStoppingEvent.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "SloughingCellKiller.hpp"
#include "CryptCellsGenerator.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "Version.hpp"
#include "FileFinder.hpp"

/*
 * A crypt invasion experiment to determine the "wild-type" behaviour...
 */
class TestWildTypeMonoclonality : public AbstractCellBasedTestSuite
{
public:

    /*
     * Here we take a crypt in dynamic equilibrium and track the ancestors of the different cells.
     */
    void TestForWildTypeMonoclonalityBehaviour() throw (Exception)
    {
        // Set parameters for this compile
        const unsigned num_experiments = 1000;  // Stop simulation after this number of experiments.
        double load_time = 300;                 // time at which each simulation must be first loaded
        double maximum_duration = 10000;        // time we are prepared to wait for each simulation to end

        // Compilation information
        std::cout << "Compiled from Chaste revision number: " << ChasteBuildInfo::GetVersionString() << "\n\n";

        // Change directory_to_copy_from according to which archive
        // (i.e. choice of geometry and cell cycle model) is used
        FileFinder directory_to_copy_from("projects/CryptInvasion/test/data/SteadyStateWT/sunter3_archive", RelativeTo::ChasteSourceRoot);
        TS_ASSERT(directory_to_copy_from.IsDir());

        CommandLineArguments* p_args = CommandLineArguments::Instance();
        unsigned argc = *(p_args->p_argc); // has the number of arguments, and
        char **argv = *(p_args->p_argv); // is a char** of them.
        std::cout << "#" << argc-1 << " arguments supplied.\n" << std::flush;
        if (argc != 2)
        {            std::cerr << "TestWildTypeMonoclonality::Please input one arguments\n"
                         "* output suffix (simulation run number)\n" << std::flush;
            return;
        }
        unsigned sim_run = (unsigned)atoi(argv[1]);

        // Set output directory
	    std::string arg1_string(argv[1]);
        std::string output_directory = "CryptInvasion_WT_" + arg1_string;

        // Create results file handler
        OutputFileHandler results_handler(output_directory, false);

        // Where the archive will be copied to
        std::string directory_to_copy_to = results_handler.GetChasteTestOutputDirectory() + output_directory;

        // Make a folder to put archives into
        std::string command = "mkdir -p " + directory_to_copy_to + "/archive";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);

        // Remove existing archives to avoid confusion!
        command = "rm -rf " + directory_to_copy_to + "/archive/*";

        // Copy archive to output directory and rename so that it is recognised by CellBasedSimulationArchiver
        command = "cp -rf --remove-destination " + directory_to_copy_from.GetAbsolutePath() + "/* " + directory_to_copy_to + "/archive";
        return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);

        // Create overall results file
        std::string overall_results_filename = "overall_results_WT.dat";
        out_stream overall_results_file = results_handler.OpenOutputFile(overall_results_filename);

        // Run experiments
        double previous_simulation_end_time = load_time;
        for (unsigned i=(sim_run-1)*num_experiments; i<sim_run*num_experiments; i++)
        {
            std::cout << "EXPERIMENT = " << i << "\n" << std::flush;
            *overall_results_file << "EXPERIMENT = " << i << "\n" << std::flush;

            // Load new initial condition for simulation
            CryptSimulation2dWithAncestorStoppingEvent* p_simulator;

            p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2dWithAncestorStoppingEvent>::Load(output_directory, previous_simulation_end_time);

            if (i==(sim_run-1)*num_experiments)
            {
                RandomNumberGenerator::Instance()->Reseed(sim_run);
            }

            p_simulator->SetOutputDirectory(output_directory);
            p_simulator->rGetCellPopulation().SetOutputCellProliferativeTypes(false);
            p_simulator->rGetCellPopulation().SetOutputCellMutationStates(false);
            p_simulator->rGetCellPopulation().SetOutputCellAncestors(false);

            // Run simulation
            p_simulator->SetCheckForStoppingEvent();
            p_simulator->SetEndTime(previous_simulation_end_time+maximum_duration);
            p_simulator->Solve(); // Every time this class calls a solve it re-initialises the ancestors and starts tracking them.

            if (fabs(SimulationTime::Instance()->GetTime() - (previous_simulation_end_time+maximum_duration))<1)
            {
                EXCEPTION("Experiment did not complete - increase maximum simulation duration");
            }

            TS_ASSERT_EQUALS(p_simulator->rGetCellPopulation().GetCellAncestors().size(), 1u);
            unsigned ancestor_cell = *(p_simulator->rGetCellPopulation().GetCellAncestors().begin());
            *overall_results_file << "Ancestor height = " << p_simulator->GetAncestorHeight(ancestor_cell) << "\n" << std::flush;
            *overall_results_file << "Monoclonality time by = " << p_simulator->GetMonoclonalityTime() - previous_simulation_end_time << "\n" << std::flush;
            *overall_results_file << "Cells lost time by = " << p_simulator->GetOriginalCellsLostTime() - previous_simulation_end_time << "\n" << std::flush;

            // Save simulation and tidy up
            CellBasedSimulationArchiver<2, CryptSimulation2dWithAncestorStoppingEvent>::Save(p_simulator);
            delete p_simulator;

            previous_simulation_end_time = SimulationTime::Instance()->GetTime();
        }

        // Close results files and tidy up
        *overall_results_file << "EXPERIMENTS COMPLETE\n" << std::flush;
        overall_results_file->close();
        WntConcentration<2>::Destroy();
    }

};

#endif /*TESTWILDTYPEMONOCLONALITY_HPP_*/
