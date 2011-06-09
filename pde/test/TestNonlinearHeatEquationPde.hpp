/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef _TESTNONLINEARHEATEQUATIONPDE_HPP_
#define _TESTNONLINEARHEATEQUATIONPDE_HPP_

#include "NonlinearEquationPde.hpp"
#include <cxxtest/TestSuite.h>

class TestNonlinearEquationPde : public CxxTest::TestSuite
{
public:

    void TestNonlinearEquationPdeMethod()
    {
        ChastePoint<1> zero1(0);
        ChastePoint<2> zero2(0,0);
        ChastePoint<3> zero3(0,0,0);
        double u = 2.0;

        NonlinearEquationPde<1> heat_equation1;
        NonlinearEquationPde<2> heat_equation2;
        NonlinearEquationPde<3> heat_equation3;

        TS_ASSERT_DELTA(heat_equation1.ComputeNonlinearSourceTerm(zero1,u),0.0,1e-12);
        TS_ASSERT_DELTA(heat_equation2.ComputeNonlinearSourceTerm(zero2,u),0.0,1e-12);
        TS_ASSERT_DELTA(heat_equation3.ComputeNonlinearSourceTerm(zero3,u),0.0,1e-12);

        // Diffusion matrices should be equal to identity * u;
        c_matrix<double, 1, 1> diff1 = heat_equation1.ComputeDiffusionTerm(zero1,u);
        c_matrix<double, 2, 2> diff2 = heat_equation2.ComputeDiffusionTerm(zero2,u);
        c_matrix<double, 3, 3> diff3 = heat_equation3.ComputeDiffusionTerm(zero3,u);

        TS_ASSERT_DELTA(diff1(0,0),u,1e-12);

        TS_ASSERT_DELTA(diff2(0,0),u,1e-12);
        TS_ASSERT_DELTA(diff2(1,1),u,1e-12);
        TS_ASSERT_DELTA(diff2(0,1),0,1e-12);

        TS_ASSERT_DELTA(diff3(0,0),u,1e-12);
        TS_ASSERT_DELTA(diff3(1,1),u,1e-12);
        TS_ASSERT_DELTA(diff3(2,2),u,1e-12);
        TS_ASSERT_DELTA(diff3(0,1),0,1e-12);
        TS_ASSERT_DELTA(diff3(0,2),0,1e-12);
        TS_ASSERT_DELTA(diff3(1,2),0,1e-12);
    }
};

#endif //_TESTNONLINEARHEATEQUATIONPDE_HPP_
