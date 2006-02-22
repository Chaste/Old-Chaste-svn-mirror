#ifndef _TESTVECTORDOUBLE_HPP_
#define _TESTVECTORDOUBLE_HPP_

#include "VectorDouble.hpp"

class TestVectorDouble : public CxxTest::TestSuite
{
	public:
	
	void TestConstructor()
	{
		VectorDouble A(3);
		TS_ASSERT_DELTA(A(1), 0.0, 0.0000000001);
	}
	
	
	
	void TestCopyConstructor()
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
	
	
	
	void TestOverloadedEqualsOperator()
	{
		VectorDouble A(2);
		double value = 5.0;
		A(1) = value;
		VectorDouble B(2);
		B = A;
		TS_ASSERT_DELTA(A(1), value, 0.0000000001);
		TS_ASSERT_DELTA(B(1), value, 0.0000000001);
	}
	void TestAddition( void )
	{
		VectorDouble A(4), B(4), C(4);
		for (int i=0; i<4; i++)
		{
			A(i) = (double)(i);
			B(i) = 2.0*(double)(i);
		}
		C=A+B;
		for (int i=0; i<4; i++)
		{
			TS_ASSERT_DELTA(C(i), 3.0*(double)(i), 1e-12);
		}
		
	}
	
	void TestSubtraction( void )
	{
		VectorDouble A(4), B(4), C(4);
		for (int i=0; i<4; i++)
		{
			A(i) = (double)(i);
			B(i) = 2.0*(double)(i);
		}
		C=A-B;
		for (int i=0; i<4; i++)
		{
			TS_ASSERT_DELTA(C(i), -1.0*(double)(i), 1e-12);
		}
		
	}
	
	void TestSize( void )
	{
		VectorDouble A(4);
		TS_ASSERT_EQUALS( A.Size(), 4);
		
	}
	
	void TestDotProduct( void )
	{
		VectorDouble A(4), B(4);
		for (int i=0; i<4; i++)
		{
			A(i) = i;
			B(i) = i;
		}
		TS_ASSERT_DELTA(A.dot(B), 14.0, 1e-12);
	}
	
	void TestResetToZero()
	{
		VectorDouble A(3);
		A(0)=10; A(1)=10; A(2)=10;
		A.ResetToZero();
		TS_ASSERT_DELTA(A(0), 0.0, 0.0000000001);
		TS_ASSERT_DELTA(A(1), 0.0, 0.0000000001);
		TS_ASSERT_DELTA(A(2), 0.0, 0.0000000001);
	}
	
	void TestVectorProduct( void )
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
	
	void TestScalarVectorAndVectorScalarMultiply()
	{
		VectorDouble A(3);
		VectorDouble B(3);
		A(0)=1; A(1)=2; A(2)=3;
		B = 3.0 * A;
		
		TS_ASSERT_DELTA(B(0), 3.0, 0.0000000001);
		TS_ASSERT_DELTA(B(1), 6.0, 0.0000000001);
		TS_ASSERT_DELTA(B(2), 9.0, 0.0000000001);
		
		A = B * 3.0;
		
		TS_ASSERT_DELTA(A(0), 9.0, 0.0000000001);
		TS_ASSERT_DELTA(A(1), 18.0, 0.0000000001);
		TS_ASSERT_DELTA(A(2), 27.0, 0.0000000001);
	}
	
	void TestL2Norm()
	{
		VectorDouble A(2);
		A(0)=3; A(1)=4;
		TS_ASSERT_DELTA(A.L2Norm(), 5, 0.000001);
	}
}; 

#endif //_TESTVECTORDOUBLE_HPP_
