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
#ifndef TESTFASTSLOWBACKWARDEULERNOBLE98_HPP_
#define TESTFASTSLOWBACKWARDEULERNOBLE98_HPP_

#include <cxxtest/TestSuite.h>
#include "FastSlowBackwardEulerNoble98.hpp"
#include "BackwardEulerNobleVargheseKohlNoble1998.hpp"
#include "ZeroStimulus.hpp"

class TestFastSlowBackwardEulerNoble98 : public CxxTest::TestSuite
{
public:
    void TestFastSlowBackwardEulerNoble98FastMode(void) throw(Exception)
    {
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus);

        BackwardEulerNobleVargheseKohlNoble1998 noble(p_stimulus);

        TS_ASSERT_EQUALS(noble.GetNumberOfStateVariables(), 22u);

        FastSlowBackwardEulerNoble98 fast_noble(p_stimulus);
        fast_noble.SetState(FAST_VARS_ONLY);

        TS_ASSERT_EQUALS(fast_noble.IsFastOnly(), true);
        TS_ASSERT_EQUALS(fast_noble.GetNumberOfStateVariables(), 3u);
        TS_ASSERT_EQUALS(fast_noble.GetNumSlowValues(), 19u);

        std::vector<double> slow_values(19);
        slow_values[0] = noble.rGetStateVariables()[1];
        slow_values[1] = noble.rGetStateVariables()[2];
        slow_values[2] = noble.rGetStateVariables()[3];
        slow_values[3] = noble.rGetStateVariables()[6];
        slow_values[4] = noble.rGetStateVariables()[7];
        slow_values[5] = noble.rGetStateVariables()[8];
        slow_values[6] = noble.rGetStateVariables()[9];
        slow_values[7] = noble.rGetStateVariables()[10];
        slow_values[8] = noble.rGetStateVariables()[11];
        slow_values[9] = noble.rGetStateVariables()[12];
        slow_values[10] = noble.rGetStateVariables()[13];
        slow_values[11] = noble.rGetStateVariables()[14];
        slow_values[12] = noble.rGetStateVariables()[15];
        slow_values[13] = noble.rGetStateVariables()[16];
        slow_values[14] = noble.rGetStateVariables()[17];
        slow_values[15] = noble.rGetStateVariables()[18];
        slow_values[16] = noble.rGetStateVariables()[19];
        slow_values[17] = noble.rGetStateVariables()[20];
        slow_values[18] = noble.rGetStateVariables()[21];

        fast_noble.SetSlowValues(slow_values);

        //Test that the fast slow calculates the ionic current correctly
        double noble_i = noble.GetIIonic();
        double fast_noble_i = fast_noble.GetIIonic();

        TS_ASSERT_DELTA( noble_i, fast_noble_i, 1e-6);

        //Test that the fast model updates the m & h gates correctly
        noble.ComputeOneStepExceptVoltage(0.0);
        fast_noble.ComputeOneStepExceptVoltage(0.0);

        TS_ASSERT_DELTA(noble.rGetStateVariables()[4], fast_noble.rGetStateVariables()[1], 1e-6);
        TS_ASSERT_DELTA(noble.rGetStateVariables()[5], fast_noble.rGetStateVariables()[2], 1e-6);


        // set the slow values as less than zero
        std::vector<double> slow_value_out_of_range(19, -0.01);
        fast_noble.AdjustOutOfRangeSlowValues(slow_value_out_of_range);
        for(unsigned i=0; i<19; i++)
        {
            TS_ASSERT_DELTA(slow_value_out_of_range[i], 0.0, 1e-9);
        }

        // set the slow values of the first eleven (gating vars or fractions)
        // to slighter greater than 1.0
        for(unsigned i=0; i<11; i++)
        {
            slow_value_out_of_range[i] = 1.01;
        }
        fast_noble.AdjustOutOfRangeSlowValues(slow_value_out_of_range);
        for(unsigned i=0; i<11; i++)
        {
            TS_ASSERT_DELTA(slow_value_out_of_range[i], 1.0, 1e-9);
        }

    }

    void TestSlowMode() throw(Exception)
    {
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus);
        
        BackwardEulerNobleVargheseKohlNoble1998 noble(p_stimulus);

        TS_ASSERT_EQUALS(noble.GetNumberOfStateVariables(), 22u);

        FastSlowBackwardEulerNoble98 slow_noble(p_stimulus);
        slow_noble.SetState(ALL_VARS);

        TS_ASSERT_EQUALS(slow_noble.IsFastOnly(), false);
        TS_ASSERT_EQUALS(slow_noble.GetNumberOfStateVariables(), 22u);

        //Test that the fast slow calculates the residual correctly
        double current_guess[12];
        for (unsigned i = 0; i < 12; ++i)
        {
            //For testing purposes use a perturbation of the initial values
            //as the current guess
            current_guess[i] = noble.rGetStateVariables()[i] + 0.001;
        }

        double residual[12];
        double slow_residual[12];

        noble.ComputeResidual(current_guess, residual);
        slow_noble.ComputeResidual(current_guess, slow_residual);

        for (unsigned i = 0; i < 12; ++i)
        {
            TS_ASSERT_DELTA(residual[i], slow_residual[i], 1e-6);
        }

        //Run two steps of solvers
        slow_noble.Compute(0.0, 0.02);
        noble.Compute(0.0, 0.02);
        for (unsigned i = 0; i < noble.GetNumberOfStateVariables() ; ++i)
        {
            TS_ASSERT_DELTA(noble.rGetStateVariables()[i],
                              slow_noble.rGetStateVariables()[i],
                              1e-14);
        }
    }
 };

#endif /*TESTFASTSLOWBACKWARDEULERNOBLE98_HPP_*/
