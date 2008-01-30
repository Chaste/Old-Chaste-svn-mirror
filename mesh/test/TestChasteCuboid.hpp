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
