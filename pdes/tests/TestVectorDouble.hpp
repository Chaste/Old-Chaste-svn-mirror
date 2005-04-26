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
	
	void testDotProduct( void )
	{
		VectorDouble A(5), B(5);
		for (int i=0; i<5; i++)
		{
			A(i) = i;
			B(i) = i;
		}
		TS_ASSERT_DELTA(A.dot(B), 30.0, 1e-12);
	}
	
	void testResetToZero()
	{
		VectorDouble A(3);
		A(0)=10; A(1)=10; A(2)=10;
		A.ResetToZero();
		TS_ASSERT_DELTA(A(0), 0.0, 0.0000000001);
		TS_ASSERT_DELTA(A(1), 0.0, 0.0000000001);
		TS_ASSERT_DELTA(A(2), 0.0, 0.0000000001);
	}
	
	void testVectorProduct( void )
	{
		VectorDouble A(3);
		VectorDouble B(3);
		VectorDouble C(3);
		
		A(0)=1; A(1)=2; A(2)=3;
		B(0)=2; B(1)=3; B(2)=4;
		
		C=A.VectorProduct(B);
		
		TS_ASSERT_DELTA(C(0), -1, 0.0000001);
		TS_ASSERT_DELTA(C(1), 2, 0.0000001);
		TS_ASSERT_DELTA(C(2), -1, 0.0000001);
		
	}
	
}; 

#endif //_TESTVECTORDOUBLE_HPP_
