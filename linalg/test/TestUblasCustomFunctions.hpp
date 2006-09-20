#ifndef TESTUBLASCUSTOMFUNCTIONS_HPP_
#define TESTUBLASCUSTOMFUNCTIONS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "UblasCustomFunctions.hpp"

class TestUblasCustomFunctions : public CxxTest::TestSuite
{

public:

    void TestDeterminant()
    {
        using namespace boost::numeric::ublas;
        c_matrix<double, 1, 1> C;
        double OneOneDeterminant;
        C(0,0)=5.6;
        OneOneDeterminant = Determinant(C);
        TS_ASSERT_DELTA( OneOneDeterminant, 5.6, 0.0000000001);
        
        c_matrix<double, 3, 3> A;
        A(0,0) = 2.4;
        A(0,1) = 5;
        A(0,2) = 5;
        A(1,0) = 5;
        A(1,1) = 6;
        A(1,2) = 7;
        A(2,0) = 6;
        A(2,1) = 8;
        A(2,2) = 9;
        double ThreeThreeDeterminant = Determinant(A);
        TS_ASSERT_DELTA( ThreeThreeDeterminant, 0.2, 0.0000000001);
        
        c_matrix<double, 2, 2> B;
        B(0,0) = 2.4;
        B(0,1) = 5;
        B(1,0) = 5;
        B(1,1) = 6;
        double TwoTwoDeterminant = Determinant(B);
        TS_ASSERT_DELTA( TwoTwoDeterminant, -10.6, 0.0000000001);
    }
    
    void TestInverse( void )
    {
        using namespace boost::numeric::ublas;
        c_matrix<double, 1, 1> C;
        C(0,0) = 8;
        c_matrix<double, 1, 1> invC;
        invC = Inverse(C);
        TS_ASSERT_DELTA(invC(0,0), 0.125, 0.0000000001);
        c_matrix<double, 3, 3> A;
        A(0,0) = 2.4;
        A(0,1) = 5;
        A(0,2) = 5;
        A(1,0) = 5;
        A(1,1) = 6;
        A(1,2) = 7;
        A(2,0) = 6;
        A(2,1) = 8;
        A(2,2) = 9;
        c_matrix<double, 3, 3> invA;
        c_matrix<double, 3, 3> invAMatlab;
        invAMatlab(0,0) = -10.00;
        invAMatlab(0,1) = -25.00;
        invAMatlab(0,2) = 25.00;
        invAMatlab(1,0) = -15.00;
        invAMatlab(1,1) = -42.00;
        invAMatlab(1,2) = 41.00;
        invAMatlab(2,0) = 20.00;
        invAMatlab(2,1) = 54.00;
        invAMatlab(2,2) = -53.00;
        invA = Inverse(A);
        for ( int i = 0; i < 3; i++)
        {
            for ( int j = 0; j < 3; j++)
            {
                TS_ASSERT_DELTA( invA(i,j), invAMatlab(i,j), 0.0000000001);
            }
        }
        c_matrix<double, 2, 2> B;
        B(0,0) = 2.4;
        B(0,1) = 5;
        B(1,0) = 5;
        B(1,1) = 6;
        c_matrix<double, 2, 2> invB;
        invB = Inverse(B);
        c_matrix<double, 2, 2> invBMatlab;
        invBMatlab(0,0) = -0.5660;
        invBMatlab(0,1) = 0.4717;
        invBMatlab(1,0) = 0.4717;
        invBMatlab(1,1) = -0.2264;
        
        for ( int i = 0; i < 2; i++)
        {
            for ( int j = 0; j < 2; j++)
            {
                TS_ASSERT_DELTA( invB(i,j), invBMatlab(i,j), 0.0001);
            }
        }
        
    }
};

#endif /*TESTUBLASCUSTOMFUNCTIONS_HPP_*/
