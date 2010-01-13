/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef TESTTIMESTEPPER_HPP_
#define TESTTIMESTEPPER_HPP_

#include <cxxtest/TestSuite.h>

#include "TimeStepper.hpp"
#include <cassert>
#include <cmath>
#include <cfloat>

class TestTimeStepper : public CxxTest::TestSuite
{
public:
    void TestOverflow()
    {
        TimeStepper stepper(0.0, DBL_MAX, DBL_EPSILON);
        stepper.mTimeStep = (unsigned)(-1);
        TS_ASSERT(!stepper.IsTimeAtEnd());
        TS_ASSERT_THROWS_THIS(stepper.AdvanceOneTimeStep(),"Time step counter has overflowed.");
    }

    void TestAdvance()
    {
        const double smidge=1e-10;

        double start_time = 0.0;
        double end_time = 2.0;
        double timestep = 3.7e-05;

        ////////////////////////////////////////////////////////
        // This is how a time stepper is normally used
        ////////////////////////////////////////////////////////
        TimeStepper my_stepper(start_time, end_time, timestep);
        while ( !my_stepper.IsTimeAtEnd() )
        {
            // do something

            my_stepper.AdvanceOneTimeStep();
        }


        ////////////////////////////////////////////////////////
        // tests
        ////////////////////////////////////////////////////////
        TS_ASSERT_THROWS_THIS(TimeStepper(end_time, start_time, timestep),"The simulation duration must be positive");

        TimeStepper stepper(start_time, end_time, timestep);

        TS_ASSERT_EQUALS( stepper.EstimateTimeSteps(),
                          (unsigned) floor((end_time - start_time)/timestep) );

        double real_time_step = timestep;
        unsigned time_step_number = 0;
        double current_time = start_time;

        /* We'll trap for stopping times that are close to the end time
         * in order to avoid having a timestep of 1e-14 (or whatever) at
         * the end in the case of rounding errors.
         */
        double close_to_end_time = end_time - smidge*timestep;

        while (current_time < end_time)
        {
            TS_ASSERT(!stepper.IsTimeAtEnd());

            time_step_number++;
            // Determine what the value time step should really be like
            double to_time = start_time+time_step_number*timestep;

            if (to_time >= close_to_end_time)
            {
                real_time_step = end_time - current_time;
                // std::cout<<"InternalSolve "<<timestep<<" "<<real_time_step<<"\n";
                to_time = end_time;
            }

            //std::cout << stepper.GetNextTimeStep()-real_time_step << std::endl;

            TS_ASSERT_EQUALS(stepper.GetNextTimeStep(), real_time_step);
            TS_ASSERT_EQUALS(stepper.GetTime(), current_time);
            TS_ASSERT_EQUALS(stepper.GetNextTime(), to_time);

            // Determine the new current time
            current_time = to_time;
            stepper.AdvanceOneTimeStep();

            TS_ASSERT_EQUALS(current_time, stepper.GetTime());
        }

        TS_ASSERT(stepper.IsTimeAtEnd());
        TS_ASSERT(stepper.GetTimeStepsElapsed()==time_step_number);
    }

    void TestEnforceConstantTimeStep() throw(Exception)
    {
        TimeStepper stepper(0.0, 1.0, 0.3); // timestep does not divide, but no checking..

        TS_ASSERT_THROWS_THIS( TimeStepper bad_const_dt_stepper(0.0, 1.0, 0.3, true),
                "TimeStepper estimate non-constant timesteps will need to be used: "
                "check timestep divides (end_time-start_time) (or divides printing timestep)" );

        TS_ASSERT_THROWS_THIS( TimeStepper bad_const_dt_stepper2(0.0, 1.0, 0.99999999, true),
                "TimeStepper estimate non-constant timesteps will need to be used: "
                "check timestep divides (end_time-start_time) (or divides printing timestep)" );

        TimeStepper const_dt_stepper(0.0, 1.0, 0.1, true);
        unsigned counter = 0;
        while (!const_dt_stepper.IsTimeAtEnd())
        {
            counter++;
            const_dt_stepper.AdvanceOneTimeStep();
        }
        TS_ASSERT_EQUALS(counter,10u);
    }

};

#endif /*TESTTIMESTEPPER_HPP_*/
