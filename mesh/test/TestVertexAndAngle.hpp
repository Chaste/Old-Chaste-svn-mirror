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

#ifndef TESTVERTEXANDANGLE_HPP_
#define TESTVERTEXANDANGLE_HPP_


#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "VertexAndAngle.hpp"
#include "Exception.hpp"

#include <cmath>

class TestVertexAndAngle : public CxxTest::TestSuite
{
public:

    void TestMethods() throw (Exception)
    {
        VertexAndAngle<3> va;

        // First check the computation for points inside each of the four quadrants

        // x>0, y>0
        va.ComputeAndSetAngle(1.0, sqrt(3.0));
        double computed_angle = va.GetAngle();
        TS_ASSERT_DELTA(computed_angle, M_PI/3.0, 1e-7);

        // x>0, y<0
        va.ComputeAndSetAngle(1.0, -sqrt(3.0));
        computed_angle = va.GetAngle();
        TS_ASSERT_DELTA(computed_angle, -M_PI/3.0, 1e-7);

        // x<0, y>0
        va.ComputeAndSetAngle(-1.0, sqrt(3.0));
        computed_angle = va.GetAngle();
        TS_ASSERT_DELTA(computed_angle, 2.0*M_PI/3.0, 1e-7);

        // x<0, y<0
        va.ComputeAndSetAngle(-1.0, -sqrt(3.0));
        computed_angle = va.GetAngle();
        TS_ASSERT_DELTA(computed_angle, -2.0*M_PI/3.0, 1e-7);

        // Now check some boundary cases

        va.ComputeAndSetAngle(1.0, 0.0);
        computed_angle = va.GetAngle();
        TS_ASSERT_DELTA(computed_angle, 0.0, 1e-7);

        va.ComputeAndSetAngle(0.0, 1.0);
        computed_angle = va.GetAngle();
        TS_ASSERT_DELTA(computed_angle, M_PI/2.0, 1e-7);

        va.ComputeAndSetAngle(-1.0, 0.0);
        computed_angle = va.GetAngle();
        TS_ASSERT_DELTA(computed_angle, M_PI, 1e-7);

        va.ComputeAndSetAngle(0.0, -1.0);
        computed_angle = va.GetAngle();
        TS_ASSERT_DELTA(computed_angle, -M_PI/2.0, 1e-7);

        TS_ASSERT_THROWS_THIS(va.ComputeAndSetAngle(0.0, 0.0),"Tried to compute polar angle of (0,0)");

        // Coverage of templated dimensions
        VertexAndAngle<1> va_1d;
        va_1d.ComputeAndSetAngle(1.0, sqrt(3.0));
        computed_angle = va_1d.GetAngle();
        TS_ASSERT_DELTA(computed_angle, M_PI/3.0, 1e-7);

        VertexAndAngle<2> va_2d;
        va_2d.ComputeAndSetAngle(1.0, sqrt(3.0));
        computed_angle = va_2d.GetAngle();
        TS_ASSERT_DELTA(computed_angle, M_PI/3.0, 1e-7);
    }
};

#endif /*TESTVERTEXANDANGLE_HPP_*/
