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
#ifndef TESTMUTATIONSPREAD_HPP_
#define TESTMUTATIONSPREAD_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cancer headers
#include "TissueSimulationArchiver.hpp"

#include "CryptSimulation2d.hpp"
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
 * A single monoclonality experiment, this generates results to make movie files.
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
class TestMutationSpread : public AbstractCancerTestSuite
{
public:

    /*
     * Here we take a crypt in a steady state and label each of the cells
     * as having an ancestor which is equivalent to their current node position.
     * We then run the crypt until it becomes monoclonal (all cells share a
     * common ancestor). We print out which cell this is and can then label
     * accordingly in a second run to generate movies.
     */
    void TestWhetherMutationsSpread() throw (Exception)
    {
        // Set experimental parameters
        double load_time = 300;          // time at which simulation must be first loaded
        double time_of_each_run = 10;    // write experiment data to files every 10 hours
        double end_of_simulation = 2000; // end time for experiment

        // Change directory_to_copy_from according to which archive (i.e. choice of geometry and cell cycle model) is used
        std::stringstream directory_to_copy_from;
        directory_to_copy_from << "projects/CellProliferation09/test/data/SteadyStateSimpleWnt/sunter3_archive";

        // Set output directory
        std::string output_directory = "MutationSpread";

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

        // Load simulation
        CryptSimulation2d* p_simulator = TissueSimulationArchiver<2, CryptSimulation2d>::Load(output_directory, load_time);

        // Set output directory for simulation
        p_simulator->SetOutputDirectory(output_directory);

        // Set simulation to output cell ancestors too
        p_simulator->SetOutputCellAncestors(true);

        MeshBasedTissue<2>* p_crypt = static_cast<MeshBasedTissue<2>*>(&p_simulator->rGetTissue());

        // Set a cell to be labelled so that we can see it in visualizer.
        // Note that we should not just pick an index and use rGetCellUsingLocationIndex(),
        // since this may correspond to a ghost node. A safer way of picking a cell is to
        // use the tissue iterator, as below.
        unsigned some_number = 52;
        AbstractTissue<2>::Iterator cell_iter = p_simulator->rGetTissue().Begin();
        for (unsigned i=0; i<some_number; i++)
        {
            ++cell_iter;
        }
        cell_iter->SetMutationState(LABELLED);
        unsigned ancestor_cell = p_crypt->GetNodeCorrespondingToCell(&(*cell_iter))->GetIndex();

        p_crypt->SetCellAncestorsToNodeIndices();

        // Save simulation and tidy up
        TissueSimulationArchiver<2, CryptSimulation2d>::Save(p_simulator);
        delete p_simulator;

        // Create results file handler
        OutputFileHandler results_handler(output_directory, false);

        // Create overall results file
        out_stream overall_results_file = results_handler.OpenOutputFile("overall_results.dat");
        *overall_results_file << "Node = " << ancestor_cell << " has been labelled\n" << std::flush;

        // Run simulation
        for (double t=load_time; t<end_of_simulation+0.5; t+=time_of_each_run)
        {
            // Load simulation
            CryptSimulation2d* p_simulator = TissueSimulationArchiver<2, CryptSimulation2d>::Load(output_directory, t);

            // Reset end time and run simulation
            p_simulator->SetEndTime(t + time_of_each_run);
            p_simulator->Solve();

            // Write results to file

            std::set<unsigned> ancestor_set = p_simulator->rGetTissue().GetCellAncestors();

            *overall_results_file << "Time = " << SimulationTime::Instance()->GetTime()
                                  << "\t" << "#Ancestors = " << ancestor_set.size() << ", Ancestor list =";

            std::set<unsigned>::iterator iter = ancestor_set.begin();
            while (iter!=ancestor_set.end())
            {
                *overall_results_file << "\t" << *iter;
                iter++;
            }
            *overall_results_file << "\n" << std::flush;

            // If crypt has become monoclonal then stop
            if (ancestor_set.size()==1u)
            {
                *overall_results_file << "Crypt is monoclonal - finished simulation\n" << std::flush;
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

#endif /*TESTMUTATIONSPREAD_HPP_*/
