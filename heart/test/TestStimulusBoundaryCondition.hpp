#ifndef TESTSTIMULUSBOUNDARYCONDITION_HPP_
#define TESTSTIMULUSBOUNDARYCONDITION_HPP_

#include <cxxtest/TestSuite.h>
#include "StimulusBoundaryCondition.hpp"
#include "InitialStimulus.hpp"

class TestStimulusBoundaryCondition : public CxxTest::TestSuite
{
public:
    void TestStimulusBoundaryConditionMethod()
    {
        ChastePoint<1> zero(0);
        InitialStimulus sq_wave(23.0, 2.0, 1.0); // magnitude, duration, start time
        StimulusBoundaryCondition<1> stim_bc(&sq_wave);
        
        PdeSimulationTime::SetTime(0.5);
        TS_ASSERT_EQUALS(0.0, stim_bc.GetValue(zero));
       
        PdeSimulationTime::SetTime(1.5);
        TS_ASSERT_EQUALS(23.0, stim_bc.GetValue(zero));
        
        PdeSimulationTime::SetTime(5.0);
        TS_ASSERT_EQUALS(0.0, stim_bc.GetValue(zero));
    }
};

#endif /*TESTSTIMULUSBOUNDARYCONDITION_HPP_*/
