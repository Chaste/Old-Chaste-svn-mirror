#ifndef TESTTIMESTEPPER_HPP_
#define TESTTIMESTEPPER_HPP_

#include <cxxtest/TestSuite.h>

#include "TimeStepper.hpp"

#include <cfloat>

class TestTimeStepper : public CxxTest::TestSuite
{
public:
    void TestConstruction()
    {
        // The provided dt should divide the simulation interval
        TS_ASSERT_THROWS_ANYTHING(TimeStepper stepper(0.0, 10.0, 1.5));
        // The end time should be greater than the start time
        TS_ASSERT_THROWS_ANYTHING(TimeStepper stepper(10.0, -10.0, 1.0));
        TS_ASSERT_THROWS_ANYTHING(TimeStepper stepper(-10.0, -10.0, 0.0));

        TS_ASSERT_THROWS_NOTHING(TimeStepper stepper(0.0, 10.0, 1.0));
        
        //edge cases
        TS_ASSERT_THROWS_NOTHING(TimeStepper stepper(0.0, 10.0, 10.0/3.0));
        TS_ASSERT_THROWS_NOTHING(TimeStepper stepper(0.0, 10.0, 10.0/3.0+DBL_EPSILON));
        TS_ASSERT_THROWS_NOTHING(TimeStepper stepper(0.0, 10.0, 10.0/3.0-DBL_EPSILON));
        TS_ASSERT_THROWS_NOTHING(TimeStepper stepper(0.0, 10.0, 10.0/3.0+8*DBL_EPSILON));
        TS_ASSERT_THROWS_NOTHING(TimeStepper stepper(0.0, 10.0, 10.0/3.0-8*DBL_EPSILON));
        // Whether factors 9--11 throw depends on the compiler...
        TS_ASSERT_THROWS_ANYTHING(TimeStepper stepper(0.0, 10.0, 10.0/3.0+12*DBL_EPSILON));
        TS_ASSERT_THROWS_ANYTHING(TimeStepper stepper(0.0, 10.0, 10.0/3.0-12*DBL_EPSILON));
    }
    
    void TestGetTotalSteps()
    {
        TimeStepper stepper(0.0, 10.0, 1.0);
        TS_ASSERT_EQUALS(stepper.GetTotalSteps(), 10u);
        
        TimeStepper stepper2(0.0, 10.0, 10.0/3.0+8*DBL_EPSILON);
        TS_ASSERT_EQUALS(stepper2.GetTotalSteps(), 3u);
        TimeStepper stepper3(0.0, 10.0, 10.0/3.0-8*DBL_EPSILON);
        TS_ASSERT_EQUALS(stepper3.GetTotalSteps(), 3u);
    }
    
    void TestGetTime()
    {
        TimeStepper stepper(1.0, 10.0, 0.01);
        TS_ASSERT_DELTA(stepper.GetTime(), 1.0, DBL_EPSILON);
        
        stepper.AdvanceOneTimeStep();
        TS_ASSERT_DELTA(stepper.GetTime(), 1.01, DBL_EPSILON);
    }
};

#endif /*TESTTIMESTEPPER_HPP_*/
