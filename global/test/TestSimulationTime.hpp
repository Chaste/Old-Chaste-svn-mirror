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
        SimulationTime *p_simulation_time = SimulationTime :: Instance();
        
        TS_ASSERT(p_simulation_time->IsStartTimeSetUp()==false);
        
        p_simulation_time->SetStartTime(0.0);
        
        TS_ASSERT(p_simulation_time->IsStartTimeSetUp()==true);
        
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 3);
        // get the time step
        TS_ASSERT_DELTA(p_simulation_time->GetTimeStep(), 3.33333333, 1e-6);
        
        // get a second instance
        // check that the time step is set correctly
        SimulationTime *p_simulation_time2 = SimulationTime :: Instance();
        TS_ASSERT_DELTA(p_simulation_time2->GetTimeStep(), 3.33333333, 1e-6);
        
        
        // check that number of time steps starts at 0
        TS_ASSERT_EQUALS(p_simulation_time->GetTimeStepsElapsed(), 0U);
        
        // increment the time
        p_simulation_time->IncrementTimeOneStep();
        
        // check the number of time steps
        TS_ASSERT_EQUALS(p_simulation_time->GetTimeStepsElapsed(), 1U);
        
        // check the simulation time from the second instance
        TS_ASSERT_DELTA(p_simulation_time2->GetDimensionalisedTime(), 3.33333333, 1e-6);
        
        // increment the time twice
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();
        
        // check the simulation time from the first instance
        TS_ASSERT_EQUALS(p_simulation_time->GetDimensionalisedTime(), 10.0);
        
        SimulationTime::Destroy();
        
        SimulationTime *p_simulation_time3 = SimulationTime :: Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time3->SetEndTimeAndNumberOfTimeSteps(10.0,5);
        TS_ASSERT_DELTA(p_simulation_time3->GetTimeStep(), 2.0, 1e-6);
        
        SimulationTime::Destroy();
        
        p_simulation_time3 = SimulationTime :: Instance();
        p_simulation_time3->SetStartTime(5.0);
        p_simulation_time3->SetEndTimeAndNumberOfTimeSteps(10.0,5);
        TS_ASSERT_DELTA(p_simulation_time3->GetTimeStep(), 1.0, 1e-6);
        
        SimulationTime::Destroy();
    }

    void TestArchiveSimulationTime()
    {
        OutputFileHandler handler("archive");
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "time.arch";
        
        // Create and archive simulation time
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);
            p_simulation_time->IncrementTimeOneStep();
            
            std::ofstream ofs(archive_filename.c_str());       
            boost::archive::text_oarchive output_arch(ofs);
            
            output_arch << static_cast<const SimulationTime&>(*p_simulation_time);
            TS_ASSERT_EQUALS(p_simulation_time->GetDimensionalisedTime(), 0.5);
            
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_simulation_time->GetDimensionalisedTime(), 1.0);
            
            SimulationTime::Destroy();
        }
        
        // Restore
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(-100.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(5.0, 5);
            TS_ASSERT_EQUALS(p_simulation_time->GetDimensionalisedTime(), -100.0);
            
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);       
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> *p_simulation_time;
            
            TS_ASSERT_EQUALS(p_simulation_time->GetDimensionalisedTime(), 0.5);
            TS_ASSERT_EQUALS(p_simulation_time->GetTimeStep(), 0.5);

            SimulationTime::Destroy();
        }
    }

    
};
#endif /*TESTSIMULATIONTIME_HPP_*/
