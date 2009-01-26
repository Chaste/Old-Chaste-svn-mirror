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


#ifndef _TESTCELLPROPERTIES_HPP_
#define _TESTCELLPROPERTIES_HPP_

#include <cxxtest/TestSuite.h>
//#include <iostream>

#include "OdeSolution.hpp"
#include "CellProperties.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

class TestCellProperties : public CxxTest::TestSuite
{
public:

    void TestExceptionalBehaviour(void)
    {
        // Check throws an exception if no data given
        std::vector<double> empty;
        TS_ASSERT_THROWS_ANYTHING(CellProperties cell_props(empty, empty));
    }

    void TestCellPhysiologicalPropertiesForRegularLr91(void)
    {
        /*
         * Set stimulus
         */
        double magnitude_of_stimulus = -80.0;
        double duration_of_stimulus  = 0.5;  // ms
        double period = 1000.0; // 1s
        double when = 100.0;
        RegularStimulus stimulus(magnitude_of_stimulus,
                                 duration_of_stimulus,
                                 period,
                                 when);

        EulerIvpOdeSolver solver;

        /*
         * Solve
         */
        double start_time = 0.0;   // ms
        double end_time = 3450.0;  // ms

        LuoRudyIModel1991OdeSystem lr91_ode_system(&solver, &stimulus);

        OdeSolution solution = lr91_ode_system.Compute(start_time, end_time);

solution.WriteToFile("", __FUNCTION__, &lr91_ode_system, "ms");


        // Now calculate the properties
        std::vector<double> voltage=solution.GetVariableAtIndex(4);
        CellProperties  cell_props(voltage, solution.rGetTimes()); // Use default threshold

//        std::cout << "Max upstroke vel: " << cell_props.GetMaxUpstrokeVelocity() << std::endl;
//        std::cout << "Cycle length: " << cell_props.GetCycleLength() << std::endl;
//        std::cout << "Max potential: " << cell_props.GetMaxPotential() << std::endl;
//        std::cout << "Min potential: " << cell_props.GetMinPotential() << std::endl;
//        std::cout << "AP amplitude: " << cell_props.GetActionPotentialAmplitude() << std::endl;
//        std::cout << "APD20: " << cell_props.GetActionPotentialDuration(20) << std::endl;
//        std::cout << "APD50: " << cell_props.GetActionPotentialDuration(50) << std::endl;
//        std::cout << "APD90: " << cell_props.GetActionPotentialDuration(90) << std::endl;

        TS_ASSERT_DELTA(cell_props.GetMaxUpstrokeVelocity(), 418.4795, 0.001);
        TS_ASSERT_DELTA(cell_props.GetCycleLength(), 1000.00, 0.01);
        TS_ASSERT_DELTA(cell_props.GetMaxPotential(), 43.1665, 0.0001);
        TS_ASSERT_DELTA(cell_props.GetMinPotential(), -84.4395, 0.0001);
        TS_ASSERT_DELTA(cell_props.GetActionPotentialAmplitude(), 127.606, 0.001);
        TS_ASSERT_DELTA(cell_props.GetActionPotentialDuration(20), 6.66416, 0.00001);
        TS_ASSERT_DELTA(cell_props.GetActionPotentialDuration(50), 271.184, 0.001);
        TS_ASSERT_DELTA(cell_props.GetActionPotentialDuration(90), 361.544, 0.001); // Should use penultimate AP
        TS_ASSERT_DELTA(cell_props.GetTimeAtMaxUpstrokeVelocity(), 3100.7300, 0.001);
    }
};

#endif //_TESTCELLPROPERTIES_HPP_
