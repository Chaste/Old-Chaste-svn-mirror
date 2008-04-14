/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TESTSTIMULUS_HPP_
#define TESTSTIMULUS_HPP_

#include <float.h>

#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"
#include "RegularStimulus.hpp"
#include "TimeStepper.hpp"
#include "SumStimulus.hpp"
#include "ZeroStimulus.hpp"
#include "MultiStimulus.hpp"

class TestStimulus : public CxxTest::TestSuite
{
public:

    void TestInitialStimulus()
    {
        double magnitude_of_stimulus = 1.0;
        double duration_of_stimulus  = 0.5;  // ms
        double when = 100.0;

    	InitialStimulus initial_at_zero(magnitude_of_stimulus, duration_of_stimulus);
                
        TS_ASSERT_EQUALS(initial_at_zero.GetStimulus(0.0), magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_at_zero.GetStimulus(duration_of_stimulus*(1-DBL_EPSILON)), 
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_at_zero.GetStimulus(duration_of_stimulus), magnitude_of_stimulus);
        
        //Made more sloppy
        TS_ASSERT_EQUALS(initial_at_zero.GetStimulus(duration_of_stimulus*(1+3*DBL_EPSILON)), 
            0.0);
        
        InitialStimulus initial_later(magnitude_of_stimulus, duration_of_stimulus, when);
        
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
        double magnitude_of_stimulus = 1.0;
        double duration_of_stimulus  = 0.002;  // ms
        double when = 0.06;
        
        InitialStimulus initial(magnitude_of_stimulus, duration_of_stimulus, when);
        TS_ASSERT_EQUALS(initial.GetStimulus(0.062),  magnitude_of_stimulus);
    }

    void TestRegularStimulus()
    {
        double magnitude_of_stimulus = 1.0;
        double duration_of_stimulus  = 0.5;  // ms
        double period = 1000.0; // 1s
        double when = 100.0;

        RegularStimulus regular_stimulus(magnitude_of_stimulus,
                                         duration_of_stimulus,
                                         period,
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
    
    void TestZeroStimulus()
    {
        ZeroStimulus zero_stim;
        TS_ASSERT_EQUALS( zero_stim.GetStimulus(1), 0);
    }
    
    void TestMultiStimulus()
    {
        MultiStimulus multi_stim;    

        // No stimulus after creation.
        TS_ASSERT_EQUALS( multi_stim.GetStimulus(1.0), 0.0);

        InitialStimulus init_stim_a(2,1,0);
        InitialStimulus init_stim_b(3,1,30);
        // RegularStimulus result at a boundary point isn't the same for Default and IntelProduction build
        // (in fact it gives different answers to "fmod" within the IntelProduction test) 
        RegularStimulus regular_stim(2.0, 1.0, 1.0/0.15, 1);
        
        multi_stim.AddStimulus(&init_stim_a);
        multi_stim.AddStimulus(&init_stim_b);
        multi_stim.AddStimulus(&regular_stim);

        TimeStepper t(0,100,1);
        while (!t.IsTimeAtEnd())
        {
            double time=t.GetTime();
            // Stimulus equals to the sum of the individual stimuli
            TS_ASSERT_EQUALS( multi_stim.GetStimulus(time),
                              init_stim_a.GetStimulus(time)+
                              init_stim_b.GetStimulus(time)+
                              regular_stim.GetStimulus(time)
                            );
            t.AdvanceOneTimeStep();
        }
        
    }
};

#endif /*TESTSTIMULUS_HPP_*/
