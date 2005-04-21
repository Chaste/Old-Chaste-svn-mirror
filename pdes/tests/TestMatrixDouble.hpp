#ifndef _TESTMATRIXDOUBLE_HPP_
#define _TESTMATRIXDOUBLE_HPP_

#include "MatrixDouble.hpp"

class TestMatrixDouble : public CxxTest::TestSuite
{
	public:
	
	void testConstructor()
	{
		MatrixDouble A(3,2);
		TS_ASSERT_DELTA(A(1,1), 0.0, 0.0000000001);
	}
	
	
	
	void testCopyConstructor()
	{
		MatrixDouble A(3,4);
		double value = 5.0;
		A(2,2) = value;
		MatrixDouble B(A);
		TS_ASSERT_DELTA(B(2,2), value, 0.0000000001);
		
		double othervalue = 3.0;
		B(2,2) = othervalue;
		TS_ASSERT_DELTA(A(2,2), value, 0.0000000001);
		TS_ASSERT_DELTA(B(2,2), othervalue, 0.0000000001);
	}
	
	
	
	void testOverloadedEqualsOperator()
	{
		MatrixDouble A(2,2);
		double value = 5.0;
		A(0,1) = value;
		MatrixDouble B(2,2);
		B = A;
		TS_ASSERT_DELTA(A(0,1), value, 0.0000000001);
		TS_ASSERT_DELTA(B(0,1), value, 0.0000000001);
	}
	
	
	
	void testIdentity()
	{
		MatrixDouble A=MatrixDouble::Identity(3);
		
		for (int i=0; i<3; i++)
		{
			for (int j=0; j<3; j++)
			{
				if (i == j)
				{
					TS_ASSERT_DELTA(A(i,j), 1.0, 0.0000000001);
				}
				else
				{
					TS_ASSERT_DELTA(A(i,j), 0.0, 0.0000000001);
				}
			}
		}
	}
	void testRows( void )
	{
		MatrixDouble A(33,3);
		TS_ASSERT_EQUALS( A.Rows(), 33);
		
	}
	void testColumns( void )
	{
		MatrixDouble A(3,7);
		TS_ASSERT_EQUALS( A.Columns(), 7);
		
	}

}; 



#endif //_TESTMATRIXDOUBLE_HPP_
