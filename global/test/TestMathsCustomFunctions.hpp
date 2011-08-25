/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTMATHSCUSTOMFUNCTIONS_HPP_
#define TESTMATHSCUSTOMFUNCTIONS_HPP_

#include <cxxtest/TestSuite.h>

#include "MathsCustomFunctions.hpp"

class TestMathsCustomFunctions : public CxxTest::TestSuite
{
public:

    void TestSmallPow() throw(Exception)
    {
        for (unsigned i=0; i<10; i++)
        {
            TS_ASSERT_DELTA(SmallPow(0.0, i), pow(0.0, double(i)), 1e-16);
            TS_ASSERT_DELTA(SmallPow(-1.67e3, i), pow(-1.67e3, double(i)), 1e14);//(1e3)^10/DBL_EPSILON
            TS_ASSERT_DELTA(SmallPow(75.0, i), pow(75.0, double(i)), 1e-16);
        }
    }

    void TestDivides() throw(Exception)
    {
        TS_ASSERT_EQUALS(Divides(0.7, 0.1),  false);
        TS_ASSERT_EQUALS(Divides(0.07, 0.1),  false);
        TS_ASSERT_EQUALS(Divides(0.1, 0.1),  true);
        TS_ASSERT_EQUALS(Divides(1e10, 1e10),  true);
        TS_ASSERT_EQUALS(Divides(5.7e10, 5.7e20),  true);
        TS_ASSERT_EQUALS(Divides(0.01, 0.1),  true);
        TS_ASSERT_EQUALS(Divides(0.01, 1.0),  true);
        TS_ASSERT_EQUALS(Divides(0.01, 10.0),  true);
        TS_ASSERT_EQUALS(Divides(0.01, 100.01),  true);

        // Note that Divides() returns false if you attempt to divide zero by a non-zero number
        TS_ASSERT_EQUALS(Divides(0.01, 0.00),  false);
    }
};

#endif /*TESTMATHSCUSTOMFUNCTIONS_HPP_*/
