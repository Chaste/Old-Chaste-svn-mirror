#ifndef _TESTVECTORDOUBLE_HPP_
#define _TESTVECTORDOUBLE_HPP_

#include "VectorDouble.hpp"

class TestVectorDouble : public CxxTest::TestSuite
{
	public:
	
	void testConstructor()
	{
		VectorDouble A(3);
		TS_ASSERT_DELTA(A(1), 0.0, 0.0000000001);
	}
	
	
	
	void testCopyConstructor()
	{
		VectorDouble A(3);
		double value = 5.0;
		A(2) = value;
		VectorDouble B(A);
		TS_ASSERT_DELTA(B(2), value, 0.0000000001);
		
		double othervalue = 3.0;
		B(2) = othervalue;
		TS_ASSERT_DELTA(A(2), value, 0.0000000001);
		TS_ASSERT_DELTA(B(2), othervalue, 0.0000000001);
	}
	
	
	
	void testOverloadedEqualsOperator()
	{
		VectorDouble A(2);
		double value = 5.0;
		A(1) = value;
		VectorDouble B(2);
		B = A;
		TS_ASSERT_DELTA(A(1), value, 0.0000000001);
		TS_ASSERT_DELTA(B(1), value, 0.0000000001);
	}
	
	void testSize( void )
	{
		VectorDouble A(12);
		TS_ASSERT_EQUALS( A.Size(), 12);
		
	}
}; 

#endif //_TESTVECTORDOUBLE_HPP_
