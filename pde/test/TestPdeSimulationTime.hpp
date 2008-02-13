#ifndef TESTPDESIMULATIONTIME_HPP_
#define TESTPDESIMULATIONTIME_HPP_

#include <cxxtest/TestSuite.h>
#include "PdeSimulationTime.hpp"

class TestPdeSimulationTime : public CxxTest::TestSuite
{
public:
    void TestTime()
    {
        PdeSimulationTime::SetTime(1.0);
        TS_ASSERT_EQUALS(PdeSimulationTime::GetTime(), 1.0);
    }
};
#endif /*TESTPDESIMULATIONTIME_HPP_*/
