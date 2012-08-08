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

#ifndef TESTPROLIFERATIONONLY_HPP_
#define TESTPROLIFERATIONONLY_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CryptSimulation2dWithCryptInvasionStoppingEvent.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "SloughingCellKiller.hpp"
#include "CryptCellsGenerator.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "Version.hpp"
#include "CellwiseData.hpp"
#include "FileFinder.hpp"

/*
 * A single crypt invasion experiment.
 *
 * \todo Use SimpleAreaBasedWntCellCycleModel
 */
class TestProliferationOnly : public AbstractCellBasedTestSuite
{
public:

    /*
     * Here we take a crypt in dynamic equilibrium and bestow a labelled mutation on a single cell
     * somewhere in the bottom 5% of the crypt. We assume that the mutation affects only the height
     * up the crypt at which proliferation ceases.
     *
     * We must choose at what height up the crypt (as a fraction of the crypt length) proliferation
     * ceases for mutant cells.
     *
     * We then simulate the crypt, keeping track of the mutant cell and its progeny.
     *
     * The simulation is stopped as soon as a stopping event occurs; this is implemented in the
     * simulation method StoppingEventHasOccurred(), which is overridden in the new class
     * CryptSimulation2dWithCryptInvasionStoppingEvent.
     *
     * We then record the initial conditions of mutation, which of the stopping events occurred
     * and how long the simulation ran for.
     */
    void TestForCryptInvasion() throw (Exception)
    {
        // Set parameters for this compile
        const unsigned max_num_experiments_per_success = 1000; // Maximum number of successive experiments to simulate per success
        const unsigned max_experiments_possible = 15000; // Never do more than this number of experiments.
        const unsigned num_successes_aim = 500; // Stop simulation after this number of experiments.
        unsigned num_successes = 0;             // No successes yet...
        const double load_time = 300;                 // time at which each simulation must be first loaded
        const double maximum_duration = 10000;        // time we are prepared to wait for each simulation to end

        // Compilation information
        std::cout << "Compiled from Chaste revision number: " << ChasteBuildInfo::GetVersionString() << "\n\n";

        // THIS IS A HACK WHICH HIJACKS PETSC TO GET ARGUMENTS INTO A TEST!
        CommandLineArguments* p_args = CommandLineArguments::Instance();
        unsigned argc = *(p_args->p_argc); // has the number of arguments, and
        char **argv = *(p_args->p_argv); // is a char** of them.
        std::cout << "#" << argc-1 << " arguments supplied.\n" << std::flush;

        if (argc != 2)
        {
            std::cerr << "TestProliferationOnly::Please input one argument\n"
                         "* fraction of crypt height where proliferation should cease for mutant cells\n" << std::flush;
            return;
        }

        double bucket_size = 0.025; // Make the mutant be introduced only very close to the base (so force averaging is near the bottom).
        double mutation_height_fraction = 0.0*bucket_size;

        // Sanity check
        assert(mutation_height_fraction >= 0.0);
        assert(mutation_height_fraction + bucket_size <= 1.0);

        double proliferation_ceiling_height_fraction = atof(argv[1]);

        // Sanity check
        assert(proliferation_ceiling_height_fraction >= 0.0);
        assert(proliferation_ceiling_height_fraction <= 1.0);

        double wnt_labelled_threshold = 1.0 - proliferation_ceiling_height_fraction;

        std::cout << "Proliferation ceiling height fraction for labelled cells = " << proliferation_ceiling_height_fraction << "\n" << std::flush;

        /*
         * Change directory_to_copy_from according to which archive
         * (i.e. choice of geometry and cell cycle model) is used.
         */
        FileFinder directory_to_copy_from("projects/CryptInvasion/test/data/SteadyStateSimpleWnt/sunter3_archive", RelativeTo::ChasteSourceRoot);
        TS_ASSERT(directory_to_copy_from.IsDir());

        // Set output directory
        std::string arg1_string(argv[1]);
        std::string output_directory = "ProliferationOnly_" + arg1_string;

        // Create results file handler
        OutputFileHandler results_handler(output_directory, false);

        // Where the archive will be copied to
        std::string directory_to_copy_to = results_handler.GetChasteTestOutputDirectory() + output_directory;

        // Make a folder to put archives into.
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
        std::string overall_results_filename = "proliferation_only_overall_results" + arg1_string + ".dat";
        out_stream overall_results_file = results_handler.OpenOutputFile(overall_results_filename);

        // Run experiments
        bool can_reuse_previous_simulation = false;
        double previous_simulation_end_time = load_time;
        unsigned max_num_experiments = max_num_experiments_per_success; // More experiments allowed after a success below.
        for (unsigned i=0; i<max_num_experiments; i++)
        {
            std::cout << "EXPERIMENT = " << i << "\n" << std::flush;
            *overall_results_file << "EXPERIMENT = " << i << "\n" << std::flush;

            // Load new initial condition for simulation
            CryptSimulation2dWithCryptInvasionStoppingEvent* p_simulator;

            /*
             * If the previous simulation ended with the mutant population being washed
             * out of the crypt, then we can re-use this as our next initial condition...
             */
            if (!can_reuse_previous_simulation)
            {
                /*
                 * ...otherwise, we need to load the original archive and run for a random
                 * length of time in order to get a new initial condition.
                 */
                p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Load(output_directory, load_time);

                // Don't bother recording data for now
                p_simulator->rGetCellPopulation().SetOutputCellVolumes(false);
                p_simulator->SetOutputNodeVelocities(false);

                // We must reseed the random number generator, which was archived along with the simulation
                RandomNumberGenerator::Instance()->Reseed(i);

                previous_simulation_end_time = load_time + 50.0*RandomNumberGenerator::Instance()->ranf();

                p_simulator->SetEndTime(previous_simulation_end_time);
                p_simulator->SetOutputDirectory(output_directory);
                p_simulator->Solve();

                CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Save(p_simulator);
            }

            p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Load(output_directory, previous_simulation_end_time);
            p_simulator->SetOutputDirectory(output_directory);

            // Record cell areas and node velocities
            p_simulator->rGetCellPopulation().SetOutputCellVolumes(true);
            p_simulator->SetOutputNodeVelocities(true);

            // Configure simulation
            p_simulator->SetEndTime(previous_simulation_end_time + maximum_duration);
            p_simulator->SetSamplingTimestepMultiple(1200);
            p_simulator->SetCheckForStoppingEvent(true);

            // Set dependence of damping constants on cell areas
            static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation()))->SetAreaBasedDampingConstant(false);

            double crypt_length = p_simulator->GetCryptHeight();

            // Compute minimum and maximum height of initial mutation
            double minimum_mutation_height = mutation_height_fraction*crypt_length;
            double maximum_mutation_height = (mutation_height_fraction + bucket_size)*crypt_length;

            // Iterate over the tissue and generate a set of real cells located in the specified bucket
            std::vector<boost::shared_ptr<Cell> > cells_in_bucket;
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_simulator->rGetCellPopulation().Begin();
                 cell_iter != p_simulator->rGetCellPopulation().End();
                 ++cell_iter)
            {
                // Set effect of mutation on Wnt threshold
                AbstractCellCycleModel* p_model = (*cell_iter)->GetCellCycleModel();
                if (dynamic_cast<SimpleWntCellCycleModel*>(p_model))
                {
                    SimpleWntCellCycleModel* p_ccm = static_cast<SimpleWntCellCycleModel*>(p_model);
                    p_ccm->SetWntLabelledThreshold(wnt_labelled_threshold);
                }
                else
                {
                    EXCEPTION("Wrong kind of cell cycle model for setting Wnt threshold.");
                }
                // Decide where the `buckets' are.
                double cell_height = p_simulator->rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[1];
                if ( cell_height >= minimum_mutation_height && cell_height < maximum_mutation_height )
                {
                    cells_in_bucket.push_back(*cell_iter);
                }
            }

            // Now choose a random cell from this set and give it a mutation
            unsigned random_index = RandomNumberGenerator::Instance()->randMod(cells_in_bucket.size());
            boost::shared_ptr<Cell> p_mutant_cell = cells_in_bucket[random_index];

            // Sanity check
            if (   p_simulator->rGetCellPopulation().GetLocationOfCellCentre(p_mutant_cell)[1] < minimum_mutation_height
                || p_simulator->rGetCellPopulation().GetLocationOfCellCentre(p_mutant_cell)[1] > maximum_mutation_height)
            {
                EXCEPTION("The height of the initial mutation is not within the prescribed limits");
            }
            // Add the 'Labelled' cell property to track this 'mutant' clone.
            p_mutant_cell->AddCellProperty(p_simulator->rGetCellPopulation().GetCellPropertyRegistry()->Get<CellLabel>());
            p_simulator->ResetForceAveraging();
            p_simulator->SetForceAveragingThreshold(maximum_mutation_height);

            // Store the initial height and corresponding location index of the mutant cell
            unsigned ancestor_cell = p_simulator->rGetCellPopulation().GetLocationIndexUsingCell(p_mutant_cell);
            double mutant_cell_height = p_simulator->rGetCellPopulation().GetLocationOfCellCentre(p_mutant_cell)[1];

            // Record details of initial mutation
            *overall_results_file << "Initial mutation height = " << mutant_cell_height << "\n" << std::flush;
            *overall_results_file << "Initial mutation proliferation ceiling height fraction = " << proliferation_ceiling_height_fraction << "\n" << std::flush;
            *overall_results_file << "Initial mutation occurred at node " << ancestor_cell << "\n" << std::flush;
            *overall_results_file << "Experiment start time = " << SimulationTime::Instance()->GetTime() << "\n" << std::flush;

            // Run simulation
            p_simulator->Solve();

            // Record actual duration of experiment (a stopping event may have occurred)
            double actual_experiment_duration = SimulationTime::Instance()->GetTime() - previous_simulation_end_time;
            *overall_results_file << "Experiment duration = " << actual_experiment_duration << "\n" << std::flush;

            // Record how many cells are in the crypt at the end of the experiment
            unsigned num_cells_at_end = p_simulator->rGetCellPopulation().GetNumRealCells();
            *overall_results_file << "Number of cells in crypt at end = " << num_cells_at_end << "\n" << std::flush;

            c_vector<double, 2> average_forces = p_simulator->GetAverageVerticalForces();
            *overall_results_file << "Average force on WT = " << average_forces[0] << "\n" << std::flush;
            *overall_results_file << "Average force on mutant = " << average_forces[1] << "\n" << std::flush;
            *overall_results_file << "Time averaged = " << p_simulator->GetTimeOverWhichForcesAveraged() << "\n" << std::flush;

            /*
             * Record what happened to the mutant population (note that this is only
             * updated every mSamplingTimestepMultiple in the Solve() method).
             */
            boost::shared_ptr<CellPropertyRegistry> p_registry = p_simulator->rGetCellPopulation().GetCellPropertyRegistry();
            unsigned num_cells = p_simulator->rGetCellPopulation().GetNumRealCells();
            unsigned num_wild_type_cells = p_registry->Get<WildTypeCellMutationState>()->GetCellCount();
            unsigned num_labelled_cells = p_registry->Get<CellLabel>()->GetCellCount();

            // Save simulation and tidy up
            CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Save(p_simulator);
            delete p_simulator;

            if (num_wild_type_cells == 0 || num_labelled_cells==num_cells)
            {
                /*
                 * If there are no healthy cells left, then the mutant progeny must have taken
                 * over the entire crypt, so we must go back to the original initial condition
                 * for the next simulation.
                 */
                *overall_results_file << "Mutant population taken over crypt\n" << std::flush;
                can_reuse_previous_simulation = false;
                previous_simulation_end_time = load_time;
                num_successes++;  // Keep track of how many successes we have had.
                max_num_experiments += max_num_experiments_per_success;   // Allow more experiments after this success...
                if (max_num_experiments > max_experiments_possible)
                {
                    max_num_experiments = max_experiments_possible;
                }
                if (num_successes >= num_successes_aim)
                {
                    // End experiments
                    break;
                }
            }
            else if (num_wild_type_cells == num_cells && num_labelled_cells==0)
            {
                /*
                 * If there are no mutant cells left, then the mutant progeny must have been
                 * washed out of the crypt, so we can use the end of this simulation as the
                 * initial condition for the next simulation.
                 */
                *overall_results_file << "Mutant population washed out of crypt\n" << std::flush;
                can_reuse_previous_simulation = true;
                previous_simulation_end_time = SimulationTime::Instance()->GetTime();
            }
            else
            {
                /*
                 * If both healthy and mutant cells are still present in the crypt at the end
                 * of the maximum allowed duration, then we must go back to the original initial
                 * condition for the next simulation.
                 */
                *overall_results_file << "No result\n" << std::flush;
                can_reuse_previous_simulation = false;
                previous_simulation_end_time = load_time;
            }
        }

        // Close results files and tidy up
        *overall_results_file << "EXPERIMENTS COMPLETE\n" << std::flush;
        overall_results_file->close();
        WntConcentration<2>::Destroy();
    }
};

#endif /*TESTPROLIFERATIONONLY_HPP_*/
