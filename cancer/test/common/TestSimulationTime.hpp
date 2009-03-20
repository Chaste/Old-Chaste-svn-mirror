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
#ifndef TESTSIMULATIONTIME_HPP_
#define TESTSIMULATIONTIME_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "OutputFileHandler.hpp"
#include "SimulationTime.hpp"

class TestSimulationTime : public CxxTest::TestSuite
{
public:
    void TestTime()
    {
        // create the simulation time object
        // set the simulation length and number of time steps
        SimulationTime* p_simulation_time = SimulationTime :: Instance();

        TS_ASSERT_EQUALS(p_simulation_time->IsStartTimeSetUp(), false);

        p_simulation_time->SetStartTime(0.0);

        TS_ASSERT_EQUALS(p_simulation_time->IsStartTimeSetUp(), true);

        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 3);
        // get the time step
        TS_ASSERT_DELTA(p_simulation_time->GetTimeStep(), 3.33333333, 1e-6);

        // get a second instance
        // check that the time step is set correctly
        SimulationTime* p_simulation_time2 = SimulationTime :: Instance();
        TS_ASSERT_DELTA(p_simulation_time2->GetTimeStep(), 3.33333333, 1e-6);


        // check that number of time steps starts at 0
        TS_ASSERT_EQUALS(p_simulation_time->GetTimeStepsElapsed(), 0u);

        // increment the time
        p_simulation_time->IncrementTimeOneStep();

        // check the number of time steps
        TS_ASSERT_EQUALS(p_simulation_time->GetTimeStepsElapsed(), 1u);

        // check the simulation time from the second instance
        TS_ASSERT_DELTA(p_simulation_time2->GetTime(), 3.33333333, 1e-6);

        // increment the time twice
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        // check the simulation time from the first instance
        TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 10.0);

        SimulationTime::Destroy();

        SimulationTime *p_simulation_time3 = SimulationTime :: Instance();
        p_simulation_time3->SetStartTime(0.0);
        p_simulation_time3->SetEndTimeAndNumberOfTimeSteps(10.0,5);
        TS_ASSERT_DELTA(p_simulation_time3->GetTimeStep(), 2.0, 1e-6);

        SimulationTime::Destroy();

        p_simulation_time3 = SimulationTime :: Instance();
        p_simulation_time3->SetStartTime(5.0);
        p_simulation_time3->SetEndTimeAndNumberOfTimeSteps(10.0,5);
        TS_ASSERT_DELTA(p_simulation_time3->GetTimeStep(), 1.0, 1e-6);

        SimulationTime::Destroy();
    }

    void TestResetTime()
    {
        // create the simulation time object
        // set the simulation length and number of time steps
        SimulationTime* p_simulation_time = SimulationTime :: Instance();
        p_simulation_time->SetStartTime(0.0);
        unsigned num_steps = 4;
        double first_end = 10.0;
        double second_end = 20.0;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(first_end, num_steps);

        for (unsigned i=0; i<num_steps; i++)
        {
            double time_should_be = i*first_end/(double)num_steps;
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), time_should_be, 1e-9);
            p_simulation_time->IncrementTimeOneStep();
        }
        //Reset the end time and number of steps
        num_steps = 20;
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(second_end, num_steps);

        for (unsigned i=0; i<num_steps; i++)
        {
            double time_should_be = first_end + i*(second_end-first_end)/(double)num_steps;
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), time_should_be, 1e-9);
            p_simulation_time->IncrementTimeOneStep();
        }

        TS_ASSERT_DELTA(p_simulation_time->GetTime(), second_end, 1e-9);


        SimulationTime::Destroy();
    }

    void TestArchiveSimulationTime()
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "time.arch";

        // Create and archive simulation time
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);
            p_simulation_time->IncrementTimeOneStep();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(2.0,6);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), 0.75, 1e-9);

            output_arch << static_cast<const SimulationTime&> (*p_simulation_time);
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), 0.75, 1e-9);

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), 1.0, 1e-9);

            SimulationTime::Destroy();
        }

        // Restore
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(-100.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(5.0, 5);
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), -100.0, 1e-9);

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> *p_simulation_time;

            TS_ASSERT_DELTA(p_simulation_time->GetTime(), 0.75,1e-9);
            TS_ASSERT_DELTA(p_simulation_time->GetTimeStep(), 0.25, 1e-9);

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), 1.0, 1e-9);

            SimulationTime::Destroy();
        }
    }


};
#endif /*TESTSIMULATIONTIME_HPP_*/
