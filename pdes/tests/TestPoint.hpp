#ifndef _TESTPOINT_HPP_
#define _TESTPOINT_HPP_

#include "Point.hpp"

class TestPoint : public CxxTest::TestSuite 
{
	public:
	/**
	 * Test that values set at a coordinate are the same when accessed.
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
				
	}
};


#endif //_TESTPOINT_HPP_
