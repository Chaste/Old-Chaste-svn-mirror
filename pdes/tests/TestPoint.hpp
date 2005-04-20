#ifndef _TESTPOINT_HPP_
#define _TESTPOINT_HPP_

#include "Point.hpp"
#include <cxxtest/TestSuite.h>

class TestPoint : public CxxTest::TestSuite 
{
	public:
	/**
	 * Test that values set at a coordinate are the same when accessed.
	 * Also check constructors.
	 * We only test dimensions from 1 to 3 inclusive.
	 */
	void testSetAndGetCoordinate(void)
	{
		Point<1> point1;
		double value = 12.0;
		int index = 0;
		point1.SetCoordinate(index, value);
		TS_ASSERT_DELTA(value, point1[index], 1e-12);
		
		Point<2> point2;
		point2.SetCoordinate(index, value);
		TS_ASSERT_DELTA(value, point2[index], 1e-12);
		index = 1; value = -13.56;
		point2.SetCoordinate(index, value);
		TS_ASSERT_DELTA(value, point2[index], 1e-12);
		
		Point<3> point3;
		index = 0;
		point3.SetCoordinate(index, value);
		TS_ASSERT_DELTA(value, point3[index], 1e-12);
		index = 1; value = 1e5;
		point3.SetCoordinate(index, value);
		TS_ASSERT_DELTA(value, point3[index], 1e-12);
		index = 2; value = 1e-5;
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
};


#endif //_TESTPOINT_HPP_
