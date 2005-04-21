#ifndef _TESTNONLINEARHEATEQUATIONPDE_HPP_
#define _TESTNONLINEARHEATEQUATIONPDE_HPP_

#include "NonlinearHeatEquationPde.hpp"

class TestNonlinearHeatEquationPde : public CxxTest::TestSuite
{
public:
	void testNonlinearHeatEquationPde()
	{
		Point<1> zero1(0);
		Point<2> zero2(0,0);
		Point<3> zero3(0,0,0);
		double u = 2.0;
		
		NonlinearHeatEquationPde<1> heatEquation1;
		NonlinearHeatEquationPde<2> heatEquation2;
		NonlinearHeatEquationPde<3> heatEquation3;
			
		TS_ASSERT_DELTA(heatEquation1.ComputeNonlinearSourceTerm(zero1,u),0.0,1e-12);
		TS_ASSERT_DELTA(heatEquation2.ComputeNonlinearSourceTerm(zero2,u),0.0,1e-12);
		TS_ASSERT_DELTA(heatEquation3.ComputeNonlinearSourceTerm(zero3,u),0.0,1e-12);

		// diffusion matrices should be equal to identity * u;	
		MatrixDouble diff1 = heatEquation1.ComputeDiffusionTerm(zero1,u);
		MatrixDouble diff2 = heatEquation2.ComputeDiffusionTerm(zero2,u);
		MatrixDouble diff3 = heatEquation3.ComputeDiffusionTerm(zero3,u);

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
