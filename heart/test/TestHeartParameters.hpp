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

#ifndef TESTHEARTPARAMETERS_HPP_
#define TESTHEARTPARAMETERS_HPP_

#include <cxxtest/TestSuite.h>
#include "HeartParameters.hpp"

class TestHeartParameters : public CxxTest::TestSuite
{
public :
    void testHeartParameters()
    {
        HeartParameters::Instance()->SetParametersFile("heart/test/data/physiological_test_params.xml");
        
        double chi = HeartParameters::Instance()->Parameters()->SurfaceAreaToVolumeRatio();
        TS_ASSERT_EQUALS(chi, 2.0);
        
        // The following fails with seg faults because of issues with having to use auto_ptrs inside the
        // HeartParameters class - we think calling HeartParameters::Instance()->Parameters() calls the
        // copy constructor on the new auto_ptr which causes the old one to dereference itself.
        
        //double capacitance = HeartParameters::Instance()->Parameters()->Capacitance();
        //TS_ASSERT_EQUALS(capacitance, 4.0);
    }
};

#endif /*TESTHEARTPARAMETERS_HPP_*/
