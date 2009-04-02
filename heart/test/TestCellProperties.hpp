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
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
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

        //Creating an artificial flat potential profile
        std::vector<double> times;
        std::vector<double> flat_v;
        for (unsigned i = 0; i<100; i++)
        {
            times.push_back(i);
            flat_v.push_back(-85.0);
        }
        CellProperties  cell_properties(flat_v, times);

        //Should throw exceptions because the cached vector of onset times (mOnsets) is empty
        TS_ASSERT_THROWS_ANYTHING(cell_properties.GetLastActionPotentialDuration(90));
        TS_ASSERT_THROWS_ANYTHING(cell_properties.GetAllActionPotentialDurations(90)[0]);

        //Should throw exceptions because upstroke was never crossed
        TS_ASSERT_THROWS_ANYTHING(cell_properties.GetTimeAtLastMaxUpstrokeVelocity());
        TS_ASSERT_THROWS_ANYTHING(cell_properties.GetLastMaxUpstrokeVelocity());

        //Now make it cross the threshold so the onset vector isn't empty any longer
        times.push_back(100);
        flat_v.push_back(20.0);

        CellProperties  new_cell_properties(flat_v, times);

        //Now this should throw an exception because the vectors of APs is empty...
        TS_ASSERT_THROWS_ANYTHING(new_cell_properties.GetLastActionPotentialDuration(90));

        //...but we can calculate peak properties for the last AP (though incomplete)
        TS_ASSERT_EQUALS(new_cell_properties.GetTimeAtLastMaxUpstrokeVelocity(), 100);
        TS_ASSERT_EQUALS(new_cell_properties.GetLastMaxUpstrokeVelocity(),105);
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
        double end_time = 3450;  // ms

        LuoRudyIModel1991OdeSystem lr91_ode_system(&solver, &stimulus);

        OdeSolution solution = lr91_ode_system.Compute(start_time, end_time);

        solution.WriteToFile("", __FUNCTION__, &lr91_ode_system, "ms");


        // Now calculate the properties
        std::vector<double> voltage=solution.GetVariableAtIndex(4);
        CellProperties  cell_props(voltage, solution.rGetTimes()); // Use default threshold
        unsigned size = cell_props.GetMaxUpstrokeVelocities().size();

        TS_ASSERT_DELTA(cell_props.GetMaxUpstrokeVelocities()[size-1], 418.4795, 0.001);
        TS_ASSERT_DELTA(cell_props.GetCycleLengths()[size-2], 1000.00, 0.01);//last apd is not finished, get cycle lengths from before
        TS_ASSERT_DELTA(cell_props.GetPeakPotentials()[size-1], 43.1665, 0.0001);
        TS_ASSERT_DELTA(cell_props.GetRestingPotentials()[size-1], -84.4395, 0.0001);
        TS_ASSERT_DELTA(cell_props.GetActionPotentialAmplitudes()[size-1], 127.606, 0.001);
        TS_ASSERT_DELTA(cell_props.GetLastActionPotentialDuration(20), 6.66416, 0.00001);
        TS_ASSERT_DELTA(cell_props.GetLastActionPotentialDuration(50), 271.184, 0.001);
        TS_ASSERT_DELTA(cell_props.GetLastActionPotentialDuration(90), 361.544, 0.001); // Should use penultimate AP
        TS_ASSERT_DELTA(cell_props.GetTimesAtMaxUpstrokeVelocity()[size-1], 3100.7300, 0.001);
    }

    /**
     * Further tests in a more tricky case.
     * The APs tested here are more bumpy and less regular.
     */
    void TestTrickyActionPotential()
    {
        std::ifstream apd_file("heart/test/data/TrickyAPD.dat");
        TS_ASSERT(apd_file.is_open());

        // Create the vectors to be passed to the CellProperties object
        std::vector<double> voltages(15001);
        std::vector<double> times(15001);
        for (unsigned i=0; i <15001; i++)
        {
            apd_file >> voltages[i];
            times[i] = i;
        }
        apd_file.close();

        CellProperties  cell_properties(voltages, times); // Use default threshold
        std::vector<double> apds = cell_properties.GetAllActionPotentialDurations(50);
        unsigned size = apds.size();

        //First check that the "GetPropertyAtLastAP" actually returns a value equal to
        // the last element of the "GetAllOfThem" vector
        TS_ASSERT_EQUALS(apds[size-1],cell_properties.GetLastActionPotentialDuration(50));
        TS_ASSERT_EQUALS(cell_properties.GetTimesAtMaxUpstrokeVelocity()[size-1], cell_properties.GetTimeAtLastMaxUpstrokeVelocity());
        TS_ASSERT_EQUALS(cell_properties.GetMaxUpstrokeVelocities()[size-1], cell_properties.GetLastMaxUpstrokeVelocity() );

        // Then check against hardcoded values (checked manually from the file)
        TS_ASSERT_DELTA(apds[0], 186.3, 0.1);
        TS_ASSERT_DELTA(apds[1], 186.3, 0.1);
        TS_ASSERT_DELTA(apds[2], 185.3, 0.1);
        TS_ASSERT_DELTA(apds[3], 184.3, 0.1);
        TS_ASSERT_DELTA(apds[4], 183.3, 0.1);
        TS_ASSERT_DELTA(apds[5], 183.3, 0.1);
        TS_ASSERT_DELTA(apds[6], 180.2, 0.1);
        TS_ASSERT_DELTA(apds[7], 177.3, 0.1);
        TS_ASSERT_DELTA(apds[8], 175.3, 0.1);
        TS_ASSERT_DELTA(apds[size-1], 175.3, 0.1);

        // Check against hardcoded resting values (checked manually from the file)
        std::vector<double> resting_values = cell_properties.GetRestingPotentials();
        size = resting_values.size();

        TS_ASSERT_DELTA(resting_values[0], -85.1, 0.1);
        TS_ASSERT_DELTA(resting_values[1], -84.6, 0.1);
        TS_ASSERT_DELTA(resting_values[2], -84.6, 0.1);
        TS_ASSERT_DELTA(resting_values[3], -84.6, 0.1);
        TS_ASSERT_DELTA(resting_values[4], -84.6, 0.1);
        TS_ASSERT_DELTA(resting_values[5], -84.7, 0.1);
        TS_ASSERT_DELTA(resting_values[6], -85.0, 0.1);
        TS_ASSERT_DELTA(resting_values[7], -85.1, 0.1);
        TS_ASSERT_DELTA(resting_values[8], -85.1, 0.1);
        TS_ASSERT_DELTA(resting_values[size-1], -85.2, 0.1);

        //This file comes from a frequency drop tissue simulation.
        //First five beats at 500 bcl and last five at 2500 bcl
        std::vector<double> cycle_lengths = cell_properties.GetCycleLengths();
        size = cycle_lengths.size();

        TS_ASSERT_DELTA(cycle_lengths[0], 500, 1);
        TS_ASSERT_DELTA(cycle_lengths[1], 500, 1);
        TS_ASSERT_DELTA(cycle_lengths[2], 500, 1);
        TS_ASSERT_DELTA(cycle_lengths[3], 500, 1);
        TS_ASSERT_DELTA(cycle_lengths[4], 500, 1);
        TS_ASSERT_DELTA(cycle_lengths[5], 2500, 1);
        TS_ASSERT_DELTA(cycle_lengths[6], 2500, 1);
        TS_ASSERT_DELTA(cycle_lengths[7], 2500, 1);
        TS_ASSERT_DELTA(cycle_lengths[8], 2500, 1);
        TS_ASSERT_DELTA(cycle_lengths[size-1], 2500, 1);
     }
};

#endif //_TESTCELLPROPERTIES_HPP_
