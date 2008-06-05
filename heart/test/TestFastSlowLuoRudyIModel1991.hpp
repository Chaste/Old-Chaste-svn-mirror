/*

Copyright (C) University of Oxford, 2008

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
#ifndef TESTFASTSLOWLUORUDYIMODEL1991_HPP_
#define TESTFASTSLOWLUORUDYIMODEL1991_HPP_

#include <cxxtest/TestSuite.h>
#include "FastSlowLuoRudyIModel1991.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "ZeroStimulus.hpp"

class TestFastSlowLuoRudyIModel1991 : public CxxTest::TestSuite
{
public:
    void TestFastLuoRudyCalculateDerivativesFastMode(void) throw(Exception)
    {
        //Set up normal cell model and evaluate Y derivatives
        ZeroStimulus stimulus;
        EulerIvpOdeSolver solver;
        double time_step = 0.01;
        
        LuoRudyIModel1991OdeSystem luo_rudy(&solver, time_step, &stimulus);
        
        std::vector<double> DY_normal(8);

        luo_rudy.EvaluateYDerivatives(0.0, luo_rudy.rGetStateVariables(), DY_normal);

        //Set up fast cell model and evaluate Y derivatives
        FastSlowLuoRudyIModel1991 fast_luo_rudy(&solver, time_step, &stimulus);
        fast_luo_rudy.SetState(FAST);
        TS_ASSERT_EQUALS(fast_luo_rudy.IsFast(), true);
        TS_ASSERT_EQUALS(fast_luo_rudy.GetNumberOfStateVariables(), 6u);
        TS_ASSERT_EQUALS(fast_luo_rudy.GetNumSlowValues(), 2u);


        std::vector<double> DY_fast(6);
        std::vector<double> slow_values(2);
        slow_values[0] = luo_rudy.rGetStateVariables()[5];
        slow_values[1] = luo_rudy.rGetStateVariables()[6];
        
        fast_luo_rudy.SetSlowValues(slow_values);
                
        fast_luo_rudy.EvaluateYDerivatives(0.0, fast_luo_rudy.rGetStateVariables(), DY_fast);
        
        //Compare the resulting Y derivatives
        for (unsigned i = 0; i < 5; ++i)
        {
            TS_ASSERT_DELTA(DY_normal[i], DY_fast[i], 1e-5);
        }
        
        TS_ASSERT_DELTA(DY_normal[7], DY_fast[5], 1e-5);
    }


    void TestFastLuoRudyCalculateDerivativesSlowMode(void) throw(Exception)
    {
        //Set up normal cell model and evaluate Y derivatives
        ZeroStimulus stimulus;
        EulerIvpOdeSolver solver;
        double time_step = 0.01;
        
        LuoRudyIModel1991OdeSystem luo_rudy(&solver, time_step, &stimulus);
        std::vector<double> DY_normal(8);
        luo_rudy.EvaluateYDerivatives(0.0, luo_rudy.rGetStateVariables(), DY_normal);

        //Set up fast cell model in slow mode and evaluate Y derivatives
        FastSlowLuoRudyIModel1991 slow_luo_rudy(&solver, time_step, &stimulus);
        slow_luo_rudy.SetState(SLOW);
        TS_ASSERT_EQUALS(slow_luo_rudy.IsFast(), false);
        TS_ASSERT_EQUALS(slow_luo_rudy.GetNumberOfStateVariables(), 8u);

        std::vector<double> DY_fast(8);
        slow_luo_rudy.EvaluateYDerivatives(0.0, slow_luo_rudy.rGetStateVariables(), DY_fast);
        
        //Compare the resulting Y derivatives
        for (unsigned i = 0; i < 7; ++i)
        {
            TS_ASSERT_DELTA(DY_normal[i], DY_fast[i], 1e-5);
        }

        std::vector<double> slow_values;
        slow_luo_rudy.GetSlowValues(slow_values);
        TS_ASSERT_DELTA(slow_values[0], luo_rudy.rGetStateVariables()[5], 1e-5);
        TS_ASSERT_DELTA(slow_values[1], luo_rudy.rGetStateVariables()[6], 1e-5);
    }
};

#endif /*TESTFASTSLOWLUORUDYIMODEL1991_HPP_*/
