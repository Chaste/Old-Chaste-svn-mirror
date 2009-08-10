/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


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

        for (unsigned M=0; M<2; M++)
        {
            for (unsigned N=0; N<2; N++)
            {
                for (unsigned P=0; P<2; P++)
                {
                    for (unsigned Q=0; Q<2; Q++)
                    {
                        TS_ASSERT_DELTA(x(M,N,P,Q), 0.0, 1e-9);
                        x(M,N,P,Q) = M+N+P+Q;
                    }
                }
            }
        }

        for (unsigned M=0; M<2; M++)
        {
            for (unsigned N=0; N<2; N++)
            {
                for (unsigned P=0; P<2; P++)
                {
                    for (unsigned Q=0; Q<2; Q++)
                    {
                        TS_ASSERT_DELTA(x(M,N,P,Q), M+N+P+Q, 1e-9);
                    }
                }
            }
        }


        FourthOrderTensor<3> y;

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA(y(M,N,P,Q), 0.0, 1e-9);
                        y(M,N,P,Q) = M+N+P+Q;
                    }
                }
            }
        }

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA(y(M,N,P,Q), M+N+P+Q, 1e-9);
                    }
                }
            }
        }

        y.Zero();

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA(y(M,N,P,Q), 0.0, 1e-9);
                    }
                }
            }
        }
    }


    // Test the first of the four possibilities for SetAsProduct
    void TestSetAsProduct0() throw(Exception)
    {
        FourthOrderTensor<3> X;
        c_matrix<double,3,3> A;

        // check throws if bad component passed in..
        FourthOrderTensor<3> Z;
        TS_ASSERT_THROWS_THIS(Z.SetAsProduct(X,A,5), "Component not 0, 1, 2, or 3");


        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                A(M,N) = M+N;
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        X(M,N,P,Q) = 3*M+N+P+Q;
                    }
                }
            }
        }

        FourthOrderTensor<3> Y;
        Y.SetAsProduct(X,A,0);

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
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
        c_matrix<double,3,3> A;

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                A(M,N) = M+N;
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        X(M,N,P,Q) = M+3*N+P+Q;
                    }
                }
            }
        }

        FourthOrderTensor<3> Y;
        Y.SetAsProduct(X,A,1);

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
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
        c_matrix<double,3,3> A;

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                A(M,N) = M+N;
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        X(M,N,P,Q) = M+N+3*P+Q;
                    }
                }
            }
        }

        FourthOrderTensor<3> Y;
        Y.SetAsProduct(X,A,2);

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
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
        c_matrix<double,3,3> A;

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                A(M,N) = M+N;
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        X(M,N,P,Q) = M+N+P+3*Q;
                    }
                }
            }
        }

        FourthOrderTensor<3> Y;
        Y.SetAsProduct(X,A,3);

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA( Y(M,N,P,Q), 15+3*(3*Q+M+N+P)+(M+N+P)*3*Q, 1e-9);
                    }
                }
            }
        }
    }
};
#endif /*TESTFOURTHORDERTENSOR_HPP_*/
