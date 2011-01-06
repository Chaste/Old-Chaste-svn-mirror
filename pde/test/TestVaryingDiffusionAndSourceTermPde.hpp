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
#ifndef _TESTVARYINGDIFFUSIONANDSOURCETERMPDE_HPP_
#define _TESTVARYINGDIFFUSIONANDSOURCETERMPDE_HPP_

#include <cxxtest/TestSuite.h>
#include "VaryingDiffusionAndSourceTermPde.hpp"

class TestVaryingDiffusionAndSourceTermPde: public CxxTest::TestSuite
{
public:
    void TestVaryingPde1D()
    {
        VaryingDiffusionAndSourceTermPde<1> pde;
        ChastePoint<1> evaluation_point(2);
        TS_ASSERT_EQUALS(pde.ComputeConstantInUSourceTerm(evaluation_point),8.0);
        c_matrix<double, 1, 1> diffusion_term=pde.ComputeDiffusionTerm(evaluation_point);
        TS_ASSERT_EQUALS(diffusion_term(0,0),4);
    }

    void TestVaryingPde2D()
    {
        VaryingDiffusionAndSourceTermPde<2> pde;
        ChastePoint<2> evaluation_point(3,4);
        TS_ASSERT_EQUALS(pde.ComputeConstantInUSourceTerm(evaluation_point),125.0);
        c_matrix<double, 2, 2> diffusion_term=pde.ComputeDiffusionTerm(evaluation_point);
        TS_ASSERT_EQUALS(diffusion_term(0,0),25.0);
        TS_ASSERT_EQUALS(diffusion_term(0,1),0.0);
        TS_ASSERT_EQUALS(diffusion_term(1,0),0.0);
        TS_ASSERT_EQUALS(diffusion_term(1,1),25.0);
    }
};

#endif //_TESTVARYINGDIFFUSIONANDSOURCETERMPDE_HPP_
