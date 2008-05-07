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


#ifndef _TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_
#define _TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <iostream>

#include "PropagationPropertiesCalculator.hpp"


class TestPropagationPropertiesCalculator : public CxxTest::TestSuite
{
public:

    void TestConductionVelocity1D(void) throw (Exception)
    {
        Hdf5DataReader simulation_data("heart/test/data/Monodomain1d",
                                       "MonodomainLR91_1d", false);
        PropagationPropertiesCalculator ppc(&simulation_data);
        double velocity=ppc.CalculateConductionVelocity(5,15,0.1);
        TS_ASSERT_DELTA(velocity, 0.0499, 0.001);
        
        // Should throw because AP does not propagate far enough in simulation time
        TS_ASSERT_THROWS_ANYTHING(ppc.CalculateConductionVelocity(5,95,0.9));
        
        // Should throw because AP does not propagate far enough in simulation time
        TS_ASSERT_THROWS_ANYTHING(ppc.CalculateConductionVelocity(90,100,0.1));
        
        TS_ASSERT_DELTA(ppc.CalculateMaximumUpstrokeVelocity(1),343.9429,0.001);
        
        TS_ASSERT_DELTA(ppc.CalculatePeakMembranePotential(5),23.4467,0.001);
        
        TS_ASSERT_DELTA(ppc.CalculateActionPotentialDuration(50,5),0,0.001);
    }
};

#endif //_TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_
