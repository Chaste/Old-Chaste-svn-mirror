#ifndef TESTSTIMULUS_HPP_
#define TESTSTIMULUS_HPP_

#include <float.h>

#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"
#include "RegularStimulus.hpp"
#include "TimeStepper.hpp"
#include "SumStimulus.hpp"

class TestStimulus : public CxxTest::TestSuite
{
    public:
    void TestInitialStimulus()
    {
        double magnitude_of_stimulus = 1.0;
        double duration_of_stimulus  = 0.5 ;  // ms
        double when = 100.0;

    	InitialStimulus initial_at_zero(magnitude_of_stimulus,
             duration_of_stimulus);
                
        TS_ASSERT_EQUALS(initial_at_zero.GetStimulus(0.0), magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_at_zero.GetStimulus(duration_of_stimulus*(1-DBL_EPSILON)), 
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_at_zero.GetStimulus(duration_of_stimulus), magnitude_of_stimulus);
        
        //Made more sloppy
        TS_ASSERT_EQUALS(initial_at_zero.GetStimulus(duration_of_stimulus*(1+3*DBL_EPSILON)), 
            0.0);
        
        InitialStimulus initial_later(magnitude_of_stimulus,
             duration_of_stimulus, when);
        
        TS_ASSERT_EQUALS(initial_later.GetStimulus(when*(1-DBL_EPSILON)), 
            0.0);            
        TS_ASSERT_EQUALS(initial_later.GetStimulus(when), magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_later.GetStimulus(when*(1+DBL_EPSILON)), 
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_later.GetStimulus((when+duration_of_stimulus)*(1-DBL_EPSILON)), 
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_later.GetStimulus(when+duration_of_stimulus), magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_later.GetStimulus((when+duration_of_stimulus)*(1+2*DBL_EPSILON)), 
            0.0);
        
        
    	 
    }
    void TestDelayedInitialStimulus()
    {
        
        //TestIonicModels::TestNoble98WithDelayedInitialStimulus() used to fail
         double magnitude_of_stimulus = 1.0;
        double duration_of_stimulus  = 0.002 ;  // ms
        double when = 0.06;
        
        InitialStimulus initial(magnitude_of_stimulus,
             duration_of_stimulus, when);
        TS_ASSERT_EQUALS(initial.GetStimulus(0.062),  magnitude_of_stimulus);

    }
    void TestRegularStimulus()
    {
        double magnitude_of_stimulus = 1.0;
        double duration_of_stimulus  = 0.5 ;  // ms
        double frequency = 1.0/1000.0; // 1Hz
        double when = 100.0;
        RegularStimulus regular_stimulus(magnitude_of_stimulus,
                                 duration_of_stimulus,
                                 frequency,
                                 when);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(0.0), 
            0.0); 
            
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.0*(1-DBL_EPSILON)), 
            0.0); 
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.0), 
            magnitude_of_stimulus); 
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.0*(1+DBL_EPSILON)), 
            magnitude_of_stimulus); 
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.5*(1-DBL_EPSILON)), 
            magnitude_of_stimulus); 
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.5), 
            magnitude_of_stimulus); 
        //Made more sloppy
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.5+(1000*DBL_EPSILON)), 
            0.0); 

        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.0*(1-DBL_EPSILON)), 
            0.0); 
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.0), 
            magnitude_of_stimulus); 
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.0*(1+DBL_EPSILON)), 
            magnitude_of_stimulus); 
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.5*(1-DBL_EPSILON)), 
            magnitude_of_stimulus); 
        //An overeager floating point optimiser may fail the following test
        //by turning the stimulus off too early
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.5), 
            magnitude_of_stimulus); 
        //Made more sloppy
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.5*(1+2*DBL_EPSILON)), 
            0.0); 
  
    }
    
    //void TestBasicFmod() removed since the exact floating point behaviour
    //is too difficult to reproduce
    
    void TestSumStimulus()
    {
        InitialStimulus r1(2,1,0);
        InitialStimulus r2(3,1,3); 
        SumStimulus s(&r1,&r2);
        
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

#endif /*TESTSTIMULUS_HPP_*/
