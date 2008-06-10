/*

Copyright (C) University of Oxford, 2008

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
#ifndef TESTREPRESENTATIVESIMULATION_HPP_
#define TESTREPRESENTATIVESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CryptSimulation2d.hpp"
#include "OutputFileHandler.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"


class TestRepresentativeSimulation : public CxxTest::TestSuite
{
public:
void TestRepresentativeSimulationForProfiling() throw (Exception)
    {
        SimulationTime::Instance()->SetStartTime(0.0);

        std::string test_to_load = "SteadyStateCrypt";
        std::string test_to_profile = "CryptProfiling";
        double t = 150;   // this is the folder and time that the stored results were archived (needed to know foldernames)
        double run_for = 10; // run for 10 hours.

        // create a new clean directory...
        OutputFileHandler file_handler(test_to_profile,true);

        // The archive needs to be copied from cancer/test/data/<test_to_profile>
        // to the testoutput directory to continue running the simulation.
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string test_data_directory = "cancer/test/data/" + test_to_load +"/";
        std::string command = "cp -Rf --remove-destination " + test_data_directory +"* "+ test_output_directory +"/" + test_to_profile + "/";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);

        CryptSimulation2d* p_simulator = CryptSimulation2d::Load(test_to_profile,t);
        p_simulator->SetEndTime(t+run_for); // start time + duration
        p_simulator->Solve();
        delete p_simulator;

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTREPRESENTATIVESIMULATION_HPP_*/
