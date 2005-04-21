#ifndef _TESTLINEARHEATEQUATION_HPP_
#define _TESTLINEARHEATEQUATION_HPP_

#include <cxxtest/TestSuite.h>
#include "LinearHeatEquationPde.hpp"

class TestLinearHeatEquationPde : public CxxTest::TestSuite
{
public:
	void testLinearHeatEquationPde()
	{
		Point<1> zero1(0);
		Point<2> zero2(0,0);
		Point<3> zero3(0,0,0);
		double u = 2.0;
		
		LinearHeatEquationPde<1> heatEquation1;
		LinearHeatEquationPde<2> heatEquation2;
		LinearHeatEquationPde<3> heatEquation3;
			
		TS_ASSERT_DELTA(heatEquation1.ComputeNonlinearSourceTerm(zero1,u),0.0,1e-12);
		TS_ASSERT_DELTA(heatEquation2.ComputeNonlinearSourceTerm(zero2,u),0.0,1e-12);
		TS_ASSERT_DELTA(heatEquation3.ComputeNonlinearSourceTerm(zero3,u),0.0,1e-12);

		// diffusion matrices should be equal to identity
		MatrixDouble diff1 = heatEquation1.ComputeDiffusionTerm(zero1);
		MatrixDouble diff2 = heatEquation2.ComputeDiffusionTerm(zero2);
		MatrixDouble diff3 = heatEquation3.ComputeDiffusionTerm(zero3);

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
