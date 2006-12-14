#ifndef _TESTPOINT_HPP_
#define _TESTPOINT_HPP_

#include <cxxtest/TestSuite.h>
#include "Point.hpp"

class TestPoint : public CxxTest::TestSuite
{
public:
    /**
     * Test that values set at a coordinate are the same when accessed.
     * Also check constructors.
     * We only test dimensions from 1 to 3 inclusive.
     */
    void TestSetAndGetCoordinate(void)
    {
    
        Point<1> point1;
        
        double value = 12.0;
        int index = 0;
        point1.SetCoordinate(index, value);
        TS_ASSERT_DELTA(value, point1[index], 1e-12);
        
        Point<2> point2;
        point2.SetCoordinate(index, value);
        TS_ASSERT_DELTA(value, point2[index], 1e-12);
        index = 1;
        value = -13.56;
        point2.SetCoordinate(index, value);
        TS_ASSERT_DELTA(value, point2[index], 1e-12);
        
        Point<3> point3;
        index = 0;
        point3.SetCoordinate(index, value);
        TS_ASSERT_DELTA(value, point3[index], 1e-12);
        index = 1;
        value = 1e5;
        point3.SetCoordinate(index, value);
        TS_ASSERT_DELTA(value, point3[index], 1e-12);
        index = 2;
        value = 1e-5;
        point3.SetCoordinate(index, value);
        TS_ASSERT_DELTA(value, point3[index], 1e-12);
        
        Point<1> point4(1);
        TS_ASSERT_DELTA(point4[0], 1, 1e-12);
        
        Point<2> point5(2,3);
        TS_ASSERT_DELTA(point5[0], 2, 1e-12);
        TS_ASSERT_DELTA(point5[1], 3, 1e-12);
        
        Point<3> point6(4,5,6);
        TS_ASSERT_DELTA(point6[0], 4, 1e-12);
        TS_ASSERT_DELTA(point6[1], 5, 1e-12);
        TS_ASSERT_DELTA(point6[2], 6, 1e-12);
        
        Point<1> point7;
        TS_ASSERT_DELTA(point7[0], 0, 1e-12);
    }
    
    void TestGetLocation(void)
    
    {
        Point<3> point1(1.0,2.0,3.0);
        
        c_vector<double, 3> &point_location = point1.rGetLocation();
        
        TS_ASSERT_EQUALS(point_location(1), 2.0);
        
        point_location(0) = 0;
    }
    
    void TestZeroDimPoint(void)
    {
        Point<0> zero_dim_point;
        TS_ASSERT_THROWS_ANYTHING(zero_dim_point[0]);
    }
    
    void TestCreateFromCvector(void)
    {
        c_vector<double, 1> location;
        location[0]=34.0;
        Point<1> point(location);
        TS_ASSERT_EQUALS(point[0], 34.0);
    }
};

#endif //_TESTPOINT_HPP_
