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
};
#endif /*TESTFOURTHORDERTENSOR_HPP_*/
