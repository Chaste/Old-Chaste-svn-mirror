#ifndef TESTSIMULATIONTIME_HPP_
#define TESTSIMULATIONTIME_HPP_
#include <cxxtest/TestSuite.h>

#include "SimulationTime.hpp"

class TestSimulationTime : public CxxTest::TestSuite
{
public:
    void TestTime()
    {
        // create the simulation time object
        // set the simulation length and number of time steps
        SimulationTime *p_simulation_time = SimulationTime :: Instance();
        
        TS_ASSERT(p_simulation_time->IsSimulationTimeSetUp()==false);
        
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 3);
        
        TS_ASSERT(p_simulation_time->IsSimulationTimeSetUp()==true);
        
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
        p_simulation_time3->SetEndTimeAndNumberOfTimeSteps(10.0,5);
        TS_ASSERT_DELTA(p_simulation_time3->GetTimeStep(), 2.0, 1e-6);
        
        SimulationTime::Destroy();
    }
    
};
#endif /*TESTSIMULATIONTIME_HPP_*/
