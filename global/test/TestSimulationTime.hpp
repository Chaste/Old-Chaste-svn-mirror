#ifndef TESTSIMULATIONTIME_HPP_
#define TESTSIMULATIONTIME_HPP_
#include <cxxtest/TestSuite.h>

#include "SimulationTime.hpp"

class TestSimulationTime : public CxxTest::TestSuite
{
public:
    void TestTime()
    {
    	SimulationTime *p_simulation_time = SimulationTime :: Instance();
    	TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.0);
    	
    	p_simulation_time->IncrementTime(0.3);
    	TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.3);
    	
    	// Second call to instance should return the same object as the first call.
    	SimulationTime *p_simulation_time2 = SimulationTime :: Instance();
    	TS_ASSERT_EQUALS(p_simulation_time2->GetTime(), 0.3);
    }
    
};
#endif /*TESTSIMULATIONTIME_HPP_*/
