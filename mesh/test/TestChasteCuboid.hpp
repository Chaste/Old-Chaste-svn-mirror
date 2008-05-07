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


#ifndef TESTCUBOID_HPP_
#define TESTCUBOID_HPP_

#include <cxxtest/TestSuite.h>
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"

class TestCuboid : public CxxTest::TestSuite
{
public:

    void TestCreationAndContained(void) throw(Exception)
    {
        ChastePoint<3> point_a(-3, -3, -3);
        ChastePoint<3> point_b(3, 3, 3);
        ChastePoint<3> point_inside(0, 0, 0);
        ChastePoint<3> point_outside(-4, -4, -4);
                        
                        
        ChasteCuboid cuboid_a_b(point_a, point_b);
        TS_ASSERT_THROWS_ANYTHING(ChasteCuboid cuboid_b_a(point_b, point_a));
        
        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(point_inside), true);
        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(point_a), true);
        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(point_b), true);

        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(point_outside), false);    
        
        // just outside point counts as inside to deal with rounding errors
        double just=3.00000000000000008882; // taken from error in cuboid mesh generation
        ChastePoint<3> just_outside(just, just, just);
        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(just_outside), true);               
    }
};
    
#endif /*TESTCUBOID_HPP_*/
