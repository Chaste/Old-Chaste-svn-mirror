/*

Copyright (C) University of Oxford, 2005-2010

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

#ifndef _TESTSENSITIVITYMODIFIERS_HPP_
#define _TESTSENSITIVITYMODIFIERS_HPP_

#include <cxxtest/TestSuite.h>

#include "Exception.hpp"
#include "SensitivityModifiers.hpp"

class TestSensitivityModifiers : public CxxTest::TestSuite
{
public:
    void TestDummyModifiers(void) throw(Exception)
    {
        DummyModifier dummymod;

        double parameter = 2;
        double returned = dummymod.calc(parameter, 0.0);

        TS_ASSERT_DELTA(parameter, returned, 1e-9);
    }

    void TestFactorModifiers(void) throw(Exception)
    {
        double factor = 2;
        FactorModifier mod(factor);

        double parameter = 2;
        double returned = mod.calc(parameter, 0.0);

        TS_ASSERT_DELTA(parameter*factor, returned, 1e-9);
    }

    void TestFixedModifiers(void) throw(Exception)
    {
        double fixed = 32;
        FixedModifier mod(fixed);

        double parameter = 2;
        double returned = mod.calc(parameter, 1.0);

        TS_ASSERT_DELTA(fixed, returned, 1e-9);
    }

    void TestTimeModifiers(void) throw(Exception)
    {
        // This class just provides an example of how to make a time modifier you might want.
        TimeModifier mod;
        double parameter = 2;

        for (unsigned i=0; i<7; i++)
        {
            double time = (double)i;
            double returned = mod.calc(parameter, time);
            TS_ASSERT_DELTA(parameter*sin(time), returned, 1e-9);
        }
    }

};


#endif //_TESTSENSITIVITYMODIFIERS_HPP_
