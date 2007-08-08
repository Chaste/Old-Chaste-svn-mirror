#ifndef TESTFOURTHORDERTENSOR_HPP_
#define TESTFOURTHORDERTENSOR_HPP_

#include <cxxtest/TestSuite.h>
#include "FourthOrderTensor.hpp"

class TestFourthOrderTensor : public CxxTest::TestSuite
{
public :
    void testFourthOrderTensor() throw(Exception)
    {
        FourthOrderTensor<2> x;
        
        for(unsigned M=0; M<2; M++)
        {
            for(unsigned N=0; N<2; N++)
            {
                for(unsigned P=0; P<2; P++)
                {
                    for(unsigned Q=0; Q<2; Q++)
                    {
                        TS_ASSERT_DELTA(x(M,N,P,Q), 0.0, 1e-9);
                        x(M,N,P,Q) = M+N+P+Q;
                    }
                }
            }
        }

        for(unsigned M=0; M<2; M++)
        {
            for(unsigned N=0; N<2; N++)
            {
                for(unsigned P=0; P<2; P++)
                {
                    for(unsigned Q=0; Q<2; Q++)
                    {
                        TS_ASSERT_DELTA(x(M,N,P,Q), M+N+P+Q, 1e-9);
                    }
                }
            }
        }


        FourthOrderTensor<3> xx;
        
        for(unsigned M=0; M<3; M++)
        {
            for(unsigned N=0; N<3; N++)
            {
                for(unsigned P=0; P<3; P++)
                {
                    for(unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA(xx(M,N,P,Q), 0.0, 1e-9);
                        xx(M,N,P,Q) = M+N+P+Q;
                    }
                }
            }
        }

        for(unsigned M=0; M<3; M++)
        {
            for(unsigned N=0; N<3; N++)
            {
                for(unsigned P=0; P<3; P++)
                {
                    for(unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA(xx(M,N,P,Q), M+N+P+Q, 1e-9);
                    }
                }
            }
        }
        
        xx.Zero();

        for(unsigned M=0; M<3; M++)
        {
            for(unsigned N=0; N<3; N++)
            {
                for(unsigned P=0; P<3; P++)
                {
                    for(unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA(xx(M,N,P,Q), 0.0, 1e-9);
                    }
                }
            }
        }
    }

    
    // Test the first of the four possibilities for SetAsProduct
    void TestSetAsProduct0() throw(Exception)
    {
        FourthOrderTensor<3> X;
        Tensor<2,3> A;
        
        // check throws if bad component passed in..
        FourthOrderTensor<3> Z;
        TS_ASSERT_THROWS_ANYTHING(Z.SetAsProduct(X,A,5));
        
        
        for(unsigned M=0; M<3; M++)
        {
            for(unsigned N=0; N<3; N++)
            {
                A[M][N] = M+N;
                for(unsigned P=0; P<3; P++)
                {
                    for(unsigned Q=0; Q<3; Q++)
                    {
                        X(M,N,P,Q) = 3*M+N+P+Q;
                    }
                }
            }
        }
        
        FourthOrderTensor<3> Y;
        Y.SetAsProduct(X,A,0);
        
        for(unsigned M=0; M<3; M++)
        {
            for(unsigned N=0; N<3; N++)
            {
                for(unsigned P=0; P<3; P++)
                {
                    for(unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA( Y(M,N,P,Q), 15+3*(Q+3*M+N+P)+(Q+N+P)*3*M, 1e-9);
                    }
                }
            }
        }
    }
    

    // Test the second of the four possibilities for SetAsProduct
    void TestSetAsProduct1() throw(Exception)
    {
        FourthOrderTensor<3> X;
        Tensor<2,3> A;
        
        for(unsigned M=0; M<3; M++)
        {
            for(unsigned N=0; N<3; N++)
            {
                A[M][N] = M+N;
                for(unsigned P=0; P<3; P++)
                {
                    for(unsigned Q=0; Q<3; Q++)
                    {
                        X(M,N,P,Q) = M+3*N+P+Q;
                    }
                }
            }
        }
        
        FourthOrderTensor<3> Y;
        Y.SetAsProduct(X,A,1);
        
        for(unsigned M=0; M<3; M++)
        {
            for(unsigned N=0; N<3; N++)
            {
                for(unsigned P=0; P<3; P++)
                {
                    for(unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA( Y(M,N,P,Q), 15+3*(Q+M+3*N+P)+(Q+M+P)*3*N, 1e-9);
                    }
                }
            }
        }
    }

    // Test the third of the four possibilities for SetAsProduct
    void TestSetAsProduct2() throw(Exception)
    {
        FourthOrderTensor<3> X;
        Tensor<2,3> A;
        
        for(unsigned M=0; M<3; M++)
        {
            for(unsigned N=0; N<3; N++)
            {
                A[M][N] = M+N;
                for(unsigned P=0; P<3; P++)
                {
                    for(unsigned Q=0; Q<3; Q++)
                    {
                        X(M,N,P,Q) = M+N+3*P+Q;
                    }
                }
            }
        }
        
        FourthOrderTensor<3> Y;
        Y.SetAsProduct(X,A,2);
        
        for(unsigned M=0; M<3; M++)
        {
            for(unsigned N=0; N<3; N++)
            {
                for(unsigned P=0; P<3; P++)
                {
                    for(unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA( Y(M,N,P,Q), 15+3*(Q+M+N+3*P)+(Q+M+N)*3*P, 1e-9);
                    }
                }
            }
        }
    }
    
    // Test the last of the four possibilities for SetAsProduct
    void TestSetAsProduct3() throw(Exception)
    {
        FourthOrderTensor<3> X;
        Tensor<2,3> A;
        
        for(unsigned M=0; M<3; M++)
        {
            for(unsigned N=0; N<3; N++)
            {
                A[M][N] = M+N;
                for(unsigned P=0; P<3; P++)
                {
                    for(unsigned Q=0; Q<3; Q++)
                    {
                        X(M,N,P,Q) = M+N+P+3*Q;
                    }
                }
            }
        }
        
        FourthOrderTensor<3> Y;
        Y.SetAsProduct(X,A,3);
        
        for(unsigned M=0; M<3; M++)
        {
            for(unsigned N=0; N<3; N++)
            {
                for(unsigned P=0; P<3; P++)
                {
                    for(unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA( Y(M,N,P,Q), 15+3*(3*Q+M+N+P)+(M+N+P)*3*Q, 1e-9);
                    }
                }
            }
        }
    }
};
#endif /*TESTFOURTHORDERTENSOR_HPP_*/
