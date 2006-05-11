#ifndef _TESTMATRIXDOUBLE_HPP_
#define _TESTMATRIXDOUBLE_HPP_

#include "MatrixDouble.hpp"
#include "VectorDouble.hpp"

class TestMatrixDouble : public CxxTest::TestSuite
{
	public:
	
	void TestConstructor()
	{
		MatrixDouble A(3,3);
		TS_ASSERT_DELTA(A(1,1), 0.0, 0.0000000001);
	}
	
	
	
	void TestCopyConstructor()
	{
		MatrixDouble A(3,3);
		double value = 5.0;
		A(2,2) = value;
		MatrixDouble B(A);
		TS_ASSERT_DELTA(B(2,2), value, 0.0000000001);
		
		double othervalue = 3.0;
		B(2,2) = othervalue;
		TS_ASSERT_DELTA(A(2,2), value, 0.0000000001);
		TS_ASSERT_DELTA(B(2,2), othervalue, 0.0000000001);
	}
	
	
	
	void TestOverloadedEqualsOperator()
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
	
	
	void TestScalarMultiplication()
	{
		MatrixDouble A(3,3);
		A(0,0) = 1.0;
		A(1,0) = 2.0;
		A = A * 3.0;
		
		for (int j=0; j<3; j++)
		{

			TS_ASSERT_DELTA(A(0,0), 3.0, 0.0000000001);
			TS_ASSERT_DELTA(A(1,0), 6.0, 0.0000000001);
			if (j>0) TS_ASSERT_DELTA(A(1,1), 0.0, 0.0000000001);
		}
		
	}
	
	void TestScalarMultiplication2()
	{
		MatrixDouble A(3,3);
		A(0,0) = 1.0;
		A(1,0) = 2.0;
		A =  3.0 * A;
		
		for (int j=0; j<3; j++)
		{

			TS_ASSERT_DELTA(A(0,0), 3.0, 0.0000000001);
			TS_ASSERT_DELTA(A(1,0), 6.0, 0.0000000001);
			if (j>0) TS_ASSERT_DELTA(A(1,1), 0.0, 0.0000000001);
		}
		
	}	
	
	void TestIsSquare( void )
	{
		MatrixDouble A(3,3);
		A(1,1)=1;
		MatrixDouble B(3,3);
		B(1,1)=1;
		
		TS_ASSERT(A.IsSquare());
		
		TS_ASSERT((B.IsSquare()));
	}
	
	void TestIdentity()
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
	void TestRows( void )
	{
		MatrixDouble A(3,3);
		TS_ASSERT_EQUALS( A.Rows(), 3);
		
	}
	void TestColumns( void )
	{
		MatrixDouble A(3,3);
		TS_ASSERT_EQUALS( A.Columns(), 3);
		
	}
	
	void TestDeterminant()
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
	
	void TestInverse( void )
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
		MatrixDouble a(3,3);
		VectorDouble b(3);
		VectorDouble c(3);
		VectorDouble matlab_calc_c(3);
		a(0,0) = 2.4;
		a(0,1) = 5;
		a(1,0) = 5;
		a(1,1) = 6;
		a(2,0) = 6;
		a(2,1) = 8;
		b(0) = 1;
		b(1) = 2;
		matlab_calc_c(0) = 12.4;
		matlab_calc_c(1) = 17.0;
		matlab_calc_c(2) = 22.0;
		c = a * b;
		for( int i = 0; i < 3; i++)
		{
			TS_ASSERT_DELTA( c(i), matlab_calc_c(i), 0.000001);
		}
        
	}
    
    void dontTestMatrixVectorMultiplicationSpeed (void)
    {
        MatrixDouble m(3,3);
        VectorDouble v(3);
        VectorDouble w(3);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++)
                m (i, j) = 3 * i + j;
            v (i) = i;
        }
	    
        
        for (int k = 0; k < 5E7; k++) {
            w = m * v;
            w = v * m;
        }
        
        for (int j = 0; j < 3; ++ j)  std::cout << w(j) << " ";
        
        std::cout << std::endl;
        
    }
    
	void TestTranspose()
	{
		MatrixDouble A(3,3), B(3,3);
		A(0,0) = 2.4;
		A(0,1) = 5;
		A(1,0) = 5;
		A(1,1) = 6;
		A(2,0) = 6;
		A(2,1) = 8;
		B=A.Transpose();
		for (int i=0; i<3; i++)
		{
			for (int j=0; j<2; j++)
			{
				TS_ASSERT_DELTA(A(i,j), B(j,i), 0.00000001);
			}
		}
	}

	void TestVectorMatrixMultiplication( void )
	{
		MatrixDouble A(3,3);
		VectorDouble b(3);
		VectorDouble c(3);
		VectorDouble matlab_calc_b(3);
		A(0,0) = 2.4;
		A(0,1) = 5;
		A(1,0) = 5;
		A(1,1) = 6;
		A(2,0) = 6;
		A(2,1) = 8;
		c(0) = 1;
		c(1) = 2;
		c(2) = 3;
		matlab_calc_b(0) = 30.4;
        matlab_calc_b(1) = 41.0;
        matlab_calc_b(2) = 0.0;
		b = c * A;
		for( int i = 0; i < 3; i++)
		{
			TS_ASSERT_DELTA( b(i), matlab_calc_b(i), 0.000001);
		}
		
	}
	
	void TestResetToZero( void )
	{
		MatrixDouble A(3,3);
		A(0,0) = 2.4;
		A(0,1) = 5;
		A(1,0) = 5;
		A(1,1) = 6;
		A(2,0) = 6;
		A(2,1) = 8;
		A.ResetToZero();
		for( int i = 0; i < 3; i++)
		{
			for( int j = 0; j<3; j++)
			{
				TS_ASSERT_DELTA( A(i,j), 0.0, 0.000001);
			}
		}
	}

	void TestTraceAndInvariants()
	{
		MatrixDouble A(1,1);
		A(0,0) = 10;
		TS_ASSERT_DELTA( A.GetTrace(), 10, 1e-12);
		TS_ASSERT_DELTA( A.GetFirstInvariant(), 10, 1e-12);

		MatrixDouble B(2,2);
		B(0,0) = 1;
		B(0,1) = 2;
		B(1,0) = 3;
		B(1,1) = 4;

		TS_ASSERT_DELTA( B.GetTrace(), 5, 1e-12);
		TS_ASSERT_DELTA( B.GetFirstInvariant(),   5, 1e-12);
		TS_ASSERT_DELTA( B.GetSecondInvariant(), -2, 1e-12);
		
		MatrixDouble C(3,3);
		C(0,0) = 1;
		C(0,1) = 2;
		C(0,2) = 3;
		C(1,0) = 4;
		C(1,1) = 5;
		C(1,2) = 6;
		C(2,0) = 7;
		C(2,1) = 8;
		C(2,2) = 9;

		TS_ASSERT_DELTA( C.GetTrace(),           15, 1e-12);
		TS_ASSERT_DELTA( C.GetFirstInvariant(),  15, 1e-12);
		TS_ASSERT_DELTA( C.GetSecondInvariant(), 18, 1e-12);
		TS_ASSERT_DELTA( C.GetThirdInvariant(),   0, 1e-12);
	}
	
	void TestMatrixMatrixMultiplication()
	{
		MatrixDouble A(4,4);
		MatrixDouble B(4,4);
	
		for(int i=0; i<3; i++)
		{	
			for(int j=0; j<4; j++)
			{
				A(i,j) = i+j;
			}
		}

		for(int i=0; i<4; i++)
		{	
			for(int j=0; j<2; j++)
			{
				B(i,j) = (i+1)*(j+1);
			}
		}

		MatrixDouble C = A*B;

		
		TS_ASSERT_EQUALS(C.Rows(), 4);
		TS_ASSERT_EQUALS(C.Columns(), 4);
		
		TS_ASSERT_DELTA( C(0,0), 20, 1e-12);
		TS_ASSERT_DELTA( C(0,1), 40, 1e-12);
		TS_ASSERT_DELTA( C(1,0), 30, 1e-12);
		TS_ASSERT_DELTA( C(1,1), 60, 1e-12);
		TS_ASSERT_DELTA( C(2,0), 40, 1e-12);
		TS_ASSERT_DELTA( C(2,1), 80, 1e-12);
	}
	
	void TestMatrixAdditionAndSubtraction()
	{
		
		MatrixDouble A(4,4);
		MatrixDouble B(4,4);
		
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<4; j++)
			{
				A(i,j) = i*j + i+j + i*(10-j) + i/(j+1);
				B(i,j) = 23 + (i+4)*j + 12*j;
			}
		}

		MatrixDouble C = A+B;
		MatrixDouble D = A-B;


		for(int i=0; i<3; i++)
		{
			for(int j=0; j<4; j++)
			{
				TS_ASSERT_DELTA( C(i,j), i*j + i+j + i*(10-j) + i/(j+1) +  23 + (i+4)*j + 12*j, 1e-12);
				TS_ASSERT_DELTA( D(i,j), i*j + i+j + i*(10-j) + i/(j+1) - (23 + (i+4)*j + 12*j), 1e-12);
			}
		}
	}
    
    void TestConvertMatrixDoubleToUblas()
    {
        MatrixDouble z(1,1);
        z(0,0)=3;
        c_matrix<double, 1,1>& z_ublas_matrix = z.rGetUblasHandle1();
        TS_ASSERT_EQUALS(z(0,0), z_ublas_matrix(0,0));
        
        MatrixDouble a(2,2);
        a(0,0)=3; a(0,1)=4; a(1,0)=2; a(1,1)=5; 
        c_matrix<double, 2,2>& a_ublas_matrix = a.rGetUblasHandle2();
        TS_ASSERT_EQUALS(a(0,0), a_ublas_matrix(0,0));
        TS_ASSERT_EQUALS(a(0,1), a_ublas_matrix(0,1));
        TS_ASSERT_EQUALS(a(1,0), a_ublas_matrix(1,0));
        TS_ASSERT_EQUALS(a(1,1), a_ublas_matrix(1,1));
        
        MatrixDouble b(3,3);
        b(0,0)=3; b(2,1)=4; b(1,2)=2; 
        c_matrix<double, 3,3>& b_ublas_matrix = b.rGetUblasHandle3();
        TS_ASSERT_EQUALS(b(0,0), b_ublas_matrix(0,0));
        TS_ASSERT_EQUALS(b(2,1), b_ublas_matrix(2,1));
        TS_ASSERT_EQUALS(b(1,2), b_ublas_matrix(1,2));
        
        MatrixDouble c(4,4);
        c(0,0)=3; c(3,1)=4; c(2,3)=2; 
        c_matrix<double, 4,4>& c_ublas_matrix = c.rGetUblasHandle4();
        TS_ASSERT_EQUALS(c(0,0), c_ublas_matrix(0,0));
        TS_ASSERT_EQUALS(c(3,1), c_ublas_matrix(3,1));
        TS_ASSERT_EQUALS(c(2,3), c_ublas_matrix(2,3));
        
        c_ublas_matrix(1,1) = 10;
        TS_ASSERT_EQUALS(c(1,1), 10);
     }
    
    
}; 

#endif //_TESTMATRIXDOUBLE_HPP_
