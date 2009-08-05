/*

Copyright (C) University of Oxford, 2005-2009

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
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus);
        
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        LuoRudyIModel1991OdeSystem luo_rudy(p_solver, p_stimulus);

        std::vector<double> DY_normal(8);

        luo_rudy.EvaluateYDerivatives(0.0, luo_rudy.rGetStateVariables(), DY_normal);

        //Set up fast cell model and evaluate Y derivatives
        FastSlowLuoRudyIModel1991 fast_luo_rudy(p_solver, p_stimulus);
        fast_luo_rudy.SetState(FAST_VARS_ONLY);
        TS_ASSERT_EQUALS(fast_luo_rudy.IsFastOnly(), true);
        TS_ASSERT_EQUALS(fast_luo_rudy.GetNumberOfStateVariables(), 4u);
        TS_ASSERT_EQUALS(fast_luo_rudy.GetNumSlowValues(), 4u);

        TS_ASSERT_DELTA(fast_luo_rudy.GetVoltage(), -84, 1.0);

        std::vector<double> DY_fast(4);
        std::vector<double> slow_values(4);
        slow_values[0] = luo_rudy.rGetStateVariables()[3];
        slow_values[1] = luo_rudy.rGetStateVariables()[5];
        slow_values[2] = luo_rudy.rGetStateVariables()[6];
        slow_values[3] = luo_rudy.rGetStateVariables()[7];

        fast_luo_rudy.SetSlowValues(slow_values);

        fast_luo_rudy.EvaluateYDerivatives(0.0, fast_luo_rudy.rGetStateVariables(), DY_fast);

        TS_ASSERT_DELTA(DY_normal[0], DY_fast[0], 1e-5);
        TS_ASSERT_DELTA(DY_normal[1], DY_fast[1], 1e-5);
        TS_ASSERT_DELTA(DY_normal[2], DY_fast[2], 1e-5);
        TS_ASSERT_DELTA(DY_normal[4], DY_fast[3], 1e-5);

        std::vector<double> slow_value_out_of_range(4);
        slow_value_out_of_range[0] = -0.01;
        slow_value_out_of_range[1] = -0.01;
        slow_value_out_of_range[2] = -0.01;
        slow_value_out_of_range[3] = -0.01;
        fast_luo_rudy.AdjustOutOfRangeSlowValues(slow_value_out_of_range);
        for(unsigned i=0; i<4; i++)
        {
            TS_ASSERT_DELTA(slow_value_out_of_range[i], 0.0, 1e-9);
        }

        slow_value_out_of_range[1] = 1.01;
        slow_value_out_of_range[2] = 1.01;
        slow_value_out_of_range[3] = 1.01;
        fast_luo_rudy.AdjustOutOfRangeSlowValues(slow_value_out_of_range);
        for(unsigned i=1; i<4; i++)
        {
            TS_ASSERT_DELTA(slow_value_out_of_range[i], 1.0, 1e-9);
        }

    }


    void TestFastLuoRudyCalculateDerivativesSlowMode(void) throw(Exception)
    {
        //Set up normal cell model and evaluate Y derivatives
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus);        
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        LuoRudyIModel1991OdeSystem luo_rudy(p_solver, p_stimulus);
        std::vector<double> DY_normal(8);
        luo_rudy.EvaluateYDerivatives(0.0, luo_rudy.rGetStateVariables(), DY_normal);

        //Set up fast cell model in slow mode and evaluate Y derivatives
        FastSlowLuoRudyIModel1991 slow_luo_rudy(p_solver, p_stimulus);
        slow_luo_rudy.SetState(ALL_VARS);
        TS_ASSERT_EQUALS(slow_luo_rudy.IsFastOnly(), false);
        TS_ASSERT_EQUALS(slow_luo_rudy.GetNumberOfStateVariables(), 8u);

        TS_ASSERT_DELTA(slow_luo_rudy.GetVoltage(), -84, 1.0);

        std::vector<double> DY_fast(8);
        slow_luo_rudy.EvaluateYDerivatives(0.0, slow_luo_rudy.rGetStateVariables(), DY_fast);

        //Compare the resulting Y derivatives
        for (unsigned i = 0; i < 7; ++i)
        {
            TS_ASSERT_DELTA(DY_normal[i], DY_fast[i], 1e-5);
        }

        std::vector<double> slow_values;
        slow_luo_rudy.GetSlowValues(slow_values);
        TS_ASSERT_DELTA(slow_values[0], luo_rudy.rGetStateVariables()[3], 1e-5);
        TS_ASSERT_DELTA(slow_values[1], luo_rudy.rGetStateVariables()[5], 1e-5);
        TS_ASSERT_DELTA(slow_values[2], luo_rudy.rGetStateVariables()[6], 1e-5);
        TS_ASSERT_DELTA(slow_values[3], luo_rudy.rGetStateVariables()[7], 1e-5);

        TS_ASSERT_DELTA(slow_values[0], slow_luo_rudy.GetIntracellularCalciumConcentration(), 1e-5);
    }

    void TestExceptionsToFastSlow(void) throw(Exception)
    {
        //Set up normal cell model
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus);        
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        LuoRudyIModel1991OdeSystem luo_rudy(p_solver, p_stimulus);

        //Set up fast-slow cell model
        FastSlowLuoRudyIModel1991 slow_luo_rudy(p_solver, p_stimulus);
        FastSlowLuoRudyIModel1991 fast_luo_rudy(p_solver, p_stimulus);

        slow_luo_rudy.SetState(ALL_VARS);
        fast_luo_rudy.SetState(FAST_VARS_ONLY);
        TS_ASSERT_THROWS_ANYTHING(luo_rudy.SetState(ALL_VARS));

        TS_ASSERT_EQUALS(slow_luo_rudy.IsFastOnly(), false);
        TS_ASSERT_THROWS_ANYTHING(luo_rudy.IsFastOnly());
        TS_ASSERT_EQUALS(fast_luo_rudy.IsFastOnly(), true);

        std::vector<double> slow_values;
        slow_luo_rudy.GetSlowValues(slow_values);
        TS_ASSERT_EQUALS(slow_luo_rudy.GetNumSlowValues(), 4U);
        fast_luo_rudy.SetSlowValues(slow_values);

        TS_ASSERT_THROWS_ANYTHING(luo_rudy.GetSlowValues(slow_values));
        TS_ASSERT_THROWS_ANYTHING(luo_rudy.GetNumSlowValues());
        TS_ASSERT_THROWS_ANYTHING(luo_rudy.SetSlowValues(slow_values));
   }
 };

#endif /*TESTFASTSLOWLUORUDYIMODEL1991_HPP_*/
