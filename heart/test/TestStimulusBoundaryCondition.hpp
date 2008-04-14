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

#ifndef TESTSTIMULUSBOUNDARYCONDITION_HPP_
#define TESTSTIMULUSBOUNDARYCONDITION_HPP_

#include <cxxtest/TestSuite.h>
#include "StimulusBoundaryCondition.hpp"
#include "InitialStimulus.hpp"

class TestStimulusBoundaryCondition : public CxxTest::TestSuite
{
public:
    void TestStimulusBoundaryConditionMethod()
    {
        ChastePoint<1> zero(0);
        InitialStimulus sq_wave(23.0, 2.0, 1.0); // magnitude, duration, start time
        StimulusBoundaryCondition<1> stim_bc(&sq_wave);
        
        PdeSimulationTime::SetTime(0.5);
        TS_ASSERT_EQUALS(0.0, stim_bc.GetValue(zero));
       
        PdeSimulationTime::SetTime(1.5);
        TS_ASSERT_EQUALS(23.0, stim_bc.GetValue(zero));
        
        PdeSimulationTime::SetTime(5.0);
        TS_ASSERT_EQUALS(0.0, stim_bc.GetValue(zero));
    }
};

#endif /*TESTSTIMULUSBOUNDARYCONDITION_HPP_*/
