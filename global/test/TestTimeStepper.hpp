#ifndef TESTTIMESTEPPER_HPP_
#define TESTTIMESTEPPER_HPP_

#include <cxxtest/TestSuite.h>

#include "TimeStepper.hpp"
#include <assert.h>

#include <cfloat>

class TestTimeStepper : public CxxTest::TestSuite
{
public:
    void TestAdvance()
    {
        const double smidge=1e-10;
        
        double startTime=0.0;
        double endTime=2.0;
        double timeStep=3.7e-05;

        TimeStepper stepper(startTime, endTime, timeStep);
        
        
        double real_time_step = timeStep;
        unsigned time_step_number = 0;
        double current_time = startTime;
        
        /* We'll trap for stopping times that are close to the end time
         * in order to avoid having a timestep of 1e-14 (or whatever) at
         * the end in the case of rounding errors.
         */
        double close_to_end_time = endTime - smidge*timeStep;
        
        while (current_time < endTime)
        {
            assert(!stepper.IsTimeAtEnd());
            
            time_step_number++;
            // Determine what the value time step should really be like
            double to_time = startTime+time_step_number*timeStep;
    
            if (to_time >= close_to_end_time)
            {
                real_time_step = endTime - current_time;
                // std::cout<<"InternalSolve "<<timeStep<<" "<<real_time_step<<"\n";
                to_time = endTime;
            }
            
            //std::cout << stepper.GetNextTimeStep()-real_time_step << std::endl;
    
            assert(stepper.GetNextTimeStep()==real_time_step);
            assert(stepper.GetTime()==current_time);
            assert(stepper.GetNextTime()==to_time);
                                
            // Determine the new current time
            current_time = to_time;
            stepper.AdvanceOneTimeStep();
    
            assert(current_time == stepper.GetTime());
        }
        
        assert(stepper.IsTimeAtEnd());

    }
};

#endif /*TESTTIMESTEPPER_HPP_*/
