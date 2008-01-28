#ifndef TESTCUBOID_HPP_
#define TESTCUBOID_HPP_

#include <cxxtest/TestSuite.h>
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"

class TestCuboid : public CxxTest::TestSuite
{
public:

    void TestCreationAndContained(void)
    {
        ChastePoint<3> point_a(-3, -3, -3);
        ChastePoint<3> point_b(3, 3, 3);
        ChastePoint<3> point_inside(0, 0, 0);
        ChastePoint<3> point_outside(-4, -4, -4);                
                        
        ChasteCuboid cuboid_a_b(point_a, point_b);
        ChasteCuboid cuboid_b_a(point_b, point_a);
        
        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(point_inside), true);
        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(point_a), true);
        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(point_b), true);

        TS_ASSERT_EQUALS(cuboid_b_a.DoesContain(point_inside), true);
        TS_ASSERT_EQUALS(cuboid_b_a.DoesContain(point_a), true);
        TS_ASSERT_EQUALS(cuboid_b_a.DoesContain(point_b), true);

        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(point_outside), false);        
        TS_ASSERT_EQUALS(cuboid_b_a.DoesContain(point_outside), false);        
    }
};
    
#endif /*TESTCUBOID_HPP_*/
