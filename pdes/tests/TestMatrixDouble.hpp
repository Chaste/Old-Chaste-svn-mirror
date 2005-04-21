#ifndef _TESTMATRIXDOUBLE_HPP_
#define _TESTMATRIXDOUBLE_HPP_

#include "MatrixDouble.hpp"
#include "VectorDouble.hpp"

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
		// TODO: If you update B now, does A also change?
	}
	
	
	void testScalarMultiplication()
	{
		for (int i=2; i<4; i++)
		{
			for (int j=1; j<3; j++)
			{
				MatrixDouble A(i,j);
				A(0,0) = 1.0;
				A(1,0) = 2.0;
				A = A * 3.0;
				TS_ASSERT_DELTA(A(0,0), 3.0, 0.0000000001);
				TS_ASSERT_DELTA(A(1,0), 6.0, 0.0000000001);
				if (j>1) TS_ASSERT_DELTA(A(1,1), 0.0, 0.0000000001);
			}
		}
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
	
	void testDeterminant()
	{
		MatrixDouble C(1,1);
		double OneOneDeterminant;
		C(0,0)=5.6;
		OneOneDeterminant = C.Determinant();
		TS_ASSERT_DELTA( OneOneDeterminant, 5.6, 0.0000000001);
		MatrixDouble A(3,3);
		A(0,0) = 2.4;
		A(0,1) = 5;
		A(0,2) = 5;
		A(1,0) = 5;
		A(1,1) = 6;
		A(1,2) = 7;
		A(2,0) = 6;
		A(2,1) = 8;
		A(2,2) = 9;
		double ThreeThreeDeterminant = A.Determinant();
		TS_ASSERT_DELTA( ThreeThreeDeterminant, 0.2, 0.0000000001);
		MatrixDouble B(2,2);
		B(0,0) = 2.4;
		B(0,1) = 5;
		B(1,0) = 5;
		B(1,1) = 6;
		double TwoTwoDeterminant = B.Determinant();
		TS_ASSERT_DELTA( TwoTwoDeterminant, -10.6, 0.0000000001);
	}
	
	void testInverse( void )
	{
		MatrixDouble C(1,1);
		C(0,0) = 8;
		MatrixDouble invC(1,1);
		invC = C.Inverse();
		TS_ASSERT_DELTA(invC(0,0), 0.125, 0.0000000001);
		MatrixDouble A(3,3);
		A(0,0) = 2.4;
		A(0,1) = 5;
		A(0,2) = 5;
		A(1,0) = 5;
		A(1,1) = 6;
		A(1,2) = 7;
		A(2,0) = 6;
		A(2,1) = 8;
		A(2,2) = 9;
		MatrixDouble invA(3,3);
		MatrixDouble invAMatlab(3,3);		
		invAMatlab(0,0) = -10.00;
		invAMatlab(0,1) = -25.00;
		invAMatlab(0,2) = 25.00;
		invAMatlab(1,0) = -15.00;
		invAMatlab(1,1) = -42.00;
		invAMatlab(1,2) = 41.00;
		invAMatlab(2,0) = 20.00;
		invAMatlab(2,1) = 54.00;
		invAMatlab(2,2) = -53.00;
		invA = A.Inverse();
		for( int i = 0; i < 3; i++)
		{
			for( int j = 0; j < 3; j++)
			{
				TS_ASSERT_DELTA( invA(i,j), invAMatlab(i,j), 0.0000000001);
			}
		}
		MatrixDouble B(2,2);
		B(0,0) = 2.4;
		B(0,1) = 5;
		B(1,0) = 5;
		B(1,1) = 6;
		MatrixDouble invB(2,2);
		invB = B.Inverse();
		MatrixDouble invBMatlab(2,2);
		invBMatlab(0,0) = -0.5660;
		invBMatlab(0,1) = 0.4717;
		invBMatlab(1,0) = 0.4717;
		invBMatlab(1,1) = -0.2264;
		for( int i = 0; i < 2; i++)
		{
			for( int j = 0; j < 2; j++)
			{
				TS_ASSERT_DELTA( invB(i,j), invBMatlab(i,j), 0.0001);
			}
		}
			
	}
	
	void TestMatrixVectorMultiplication( void )
	{
		MatrixDouble A(3,2);
		VectorDouble b(2);
		VectorDouble c(3);
		VectorDouble matlab_calc_c(3);
		A(0,0) = 2.4;
		A(0,1) = 5;
		A(1,0) = 5;
		A(1,1) = 6;
		A(2,0) = 6;
		A(2,1) = 8;
		b(0) = 1;
		b(1) = 2;
		matlab_calc_c(0) = 12.4;
		matlab_calc_c(1) = 17.0;
		matlab_calc_c(2) = 22.0;
		c = A * b;
		for( int i = 0; i < 3; i++)
		{
			TS_ASSERT_DELTA( c(i), matlab_calc_c(i), 0.000001);
		}
		
	}

}; 



#endif //_TESTMATRIXDOUBLE_HPP_
