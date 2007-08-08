#ifndef _TESTLINEARHEATEQUATION_HPP_
#define _TESTLINEARHEATEQUATION_HPP_

#include <cxxtest/TestSuite.h>
#include "LinearHeatEquationPde.hpp"
#include "Node.hpp"

class TestLinearHeatEquationPde : public CxxTest::TestSuite
{
public:
    void TestComputeLinearSourceTermAtNode()
    {
        Node<2> zero(0);
        LinearHeatEquationPde<2> heat_equation;
        
        TS_ASSERT_DELTA(heat_equation.ComputeLinearSourceTermAtNode(zero), 1.0, 1e-12);
    }
    
    void TestLinearHeatEquationPdeMethod()
    {
        ChastePoint<1> zero1(0);
        ChastePoint<2> zero2(0,0);
        ChastePoint<3> zero3(0,0,0);
        double u = 2.0;
        
        LinearHeatEquationPde<1> heat_equation1;
        LinearHeatEquationPde<2> heat_equation2;
        LinearHeatEquationPde<3> heat_equation3;
        
        TS_ASSERT_DELTA(heat_equation1.ComputeNonlinearSourceTerm(zero1,u),0.0,1e-12);
        TS_ASSERT_DELTA(heat_equation2.ComputeNonlinearSourceTerm(zero2,u),0.0,1e-12);
        TS_ASSERT_DELTA(heat_equation3.ComputeNonlinearSourceTerm(zero3,u),0.0,1e-12);
        
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
