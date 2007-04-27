#ifndef TESTSUMSTIMULUS_HPP_
#define TESTSUMSTIMULUS_HPP_


#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"
#include "TimeStepper.hpp"
#include "SumStimulus.hpp"

class TestSumStimulus : public CxxTest::TestSuite
{
    public:
    void TestSum()
    {
        InitialStimulus r1(2,1,0);
        InitialStimulus r2(3,1,3); 
        SumStimulus s(r1,r2);
        
        TimeStepper t(0,10,1);
        while (!t.IsTimeAtEnd())
        {
            TS_ASSERT_EQUALS(  s.GetStimulus(t.GetTime()),
                              r1.GetStimulus(t.GetTime())+
                              r2.GetStimulus(t.GetTime()) );
            t.AdvanceOneTimeStep();
        }
    }
};

#endif /*TESTSUMSTIMULUS_HPP_*/
