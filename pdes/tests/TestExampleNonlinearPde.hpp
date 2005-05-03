#ifndef _TESTEXAMPLENONLINEARPDE_HPP_
#define _TESTEXAMPLENONLINEARPDE_HPP_

#include <cxxtest/TestSuite.h>
#include "ExampleNonlinearPde.hpp"

class TestExampleNonlinearPde : public CxxTest::TestSuite
{
public:
	void TestExampleNonlinearPdeMethod()
	{
		Point<1> zero(0);
		double u = 2.0;
		ExampleNonlinearPde example_equation;
					
		TS_ASSERT_DELTA(example_equation.ComputeNonlinearSourceTerm(zero,u), 0.0, 1e-12);

		// Diffusion matrix should be equal to 2
		MatrixDouble diff = example_equation.ComputeDiffusionTerm(zero,u);

		TS_ASSERT_DELTA(diff(0,0), 2, 1e-12);
 
	}
};

#endif //_TESTEXAMPLENONLINEARPDE_HPP_
