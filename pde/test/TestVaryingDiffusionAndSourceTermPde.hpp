#ifndef _TESTVARYINGDIFFUSIONANDSOURCETERMPDE_HPP_
#define _TESTVARYINGDIFFUSIONANDSOURCETERMPDE_HPP_

#include <cxxtest/TestSuite.h>
#include "VaryingDiffusionAndSourceTermPde.hpp"

class TestVaryingDiffusionAndSourceTermPde: public CxxTest::TestSuite
{
public:
    void TestVaryingPde1D ( void )
    {
        VaryingDiffusionAndSourceTermPde<1> pde;
        ChastePoint<1> evaluation_point(2);
        TS_ASSERT_EQUALS(pde.ComputeNonlinearSourceTerm(evaluation_point,1.0),0.0);
        TS_ASSERT_EQUALS(pde.ComputeLinearSourceTerm(evaluation_point),8.0);
        c_matrix<double, 1, 1> diffusion_term=pde.ComputeDiffusionTerm(evaluation_point);
        TS_ASSERT_EQUALS(diffusion_term(0,0),4);
    }
    
    void TestVaryingPde2D ( void )
    {
        VaryingDiffusionAndSourceTermPde<2> pde;
        ChastePoint<2> evaluation_point(3,4);
        TS_ASSERT_EQUALS(pde.ComputeNonlinearSourceTerm(evaluation_point,1.0),0.0);
        TS_ASSERT_EQUALS(pde.ComputeLinearSourceTerm(evaluation_point),125.0);
        c_matrix<double, 2, 2> diffusion_term=pde.ComputeDiffusionTerm(evaluation_point);
        TS_ASSERT_EQUALS(diffusion_term(0,0),25.0);
        TS_ASSERT_EQUALS(diffusion_term(0,1),0.0);
        TS_ASSERT_EQUALS(diffusion_term(1,0),0.0);
        TS_ASSERT_EQUALS(diffusion_term(1,1),25.0);
    }
    
};

#endif //_TESTVARYINGDIFFUSIONANDSOURCETERMPDE_HPP_
