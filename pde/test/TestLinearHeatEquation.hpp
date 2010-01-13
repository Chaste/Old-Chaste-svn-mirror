/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef _TESTLINEARHEATEQUATION_HPP_
#define _TESTLINEARHEATEQUATION_HPP_

#include <cxxtest/TestSuite.h>
#include "SimplePoissonEquation.hpp"
#include "Node.hpp"

class TestSimplePoissonEquation : public CxxTest::TestSuite
{
public:
    void TestComputeConstantInUSourceTermAtNode()
    {
        Node<2> zero(0);
        SimplePoissonEquation<2,2> heat_equation;

        TS_ASSERT_DELTA(heat_equation.ComputeConstantInUSourceTermAtNode(zero), 1.0, 1e-12);
    }

    void TestSimplePoissonEquationMethod()
    {
        ChastePoint<1> zero1(0);
        ChastePoint<2> zero2(0,0);
        ChastePoint<3> zero3(0,0,0);

        SimplePoissonEquation<1,1> heat_equation1;
        SimplePoissonEquation<2,2> heat_equation2;
        SimplePoissonEquation<3,3> heat_equation3;

        // diffusion matrices should be equal to identity
        c_matrix<double,1,1> diff1 = heat_equation1.ComputeDiffusionTerm(zero1);
        c_matrix<double,2,2> diff2 = heat_equation2.ComputeDiffusionTerm(zero2);
        c_matrix<double,3,3> diff3 = heat_equation3.ComputeDiffusionTerm(zero3);

        TS_ASSERT_DELTA(diff1(0,0),1,1e-12);

        TS_ASSERT_DELTA(diff2(0,0),1,1e-12);
        TS_ASSERT_DELTA(diff2(1,1),1,1e-12);
        TS_ASSERT_DELTA(diff2(0,1),0,1e-12);

        TS_ASSERT_DELTA(diff3(0,0),1,1e-12);
        TS_ASSERT_DELTA(diff3(1,1),1,1e-12);
        TS_ASSERT_DELTA(diff3(2,2),1,1e-12);
        TS_ASSERT_DELTA(diff3(0,1),0,1e-12);
        TS_ASSERT_DELTA(diff3(0,2),0,1e-12);
        TS_ASSERT_DELTA(diff3(1,2),0,1e-12);
    }
};


#endif //_TESTLINEARHEATEQUATION_HPP_
