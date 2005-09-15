#ifndef _TESTABSTRACTODESYSTEM_HPP_
#define _TESTABSTRACTODESYSTEM_HPP_

// TestAbstractOdeSystem.hpp

#include <cmath>
#include <iostream>
#include <vector>
#include "AbstractOdeSystem.hpp"
#include "TestOde1.hpp"
#include "TestOde2.hpp"
#include "TestOde3.hpp"

// Tolerance for tests
double tol=0.01;

class TestAbstractOdeSystem : public CxxTest::TestSuite
{
	public:
		
	void TestOdeSystemOne(void)
	{
		// pointer to TestOde1 class		
		TestOde1 ode1;
		// Yprime
		std::vector<double> YPrime;
        YPrime = ode1.EvaluateYDerivatives(1.0, ode1.mInitialConditions);
		TS_ASSERT_DELTA(YPrime[0],1.0,tol);
	}
	
	
	void TestOdeSystemTwo(void)
	{
		TestOde2 ode2;
		std::vector<double> YPrime;
		YPrime = ode2.EvaluateYDerivatives(2.0, ode2.mInitialConditions);
		TS_ASSERT_DELTA(YPrime[0],8.0,tol);
	}
	
	void TestOdeSystemThree(void)
	{
		TestOde3 ode3;
		std::vector<double> YPrime;
		YPrime = ode3.EvaluateYDerivatives(2.0, ode3.mInitialConditions);
		TS_ASSERT_DELTA(YPrime[0],8.0,tol);
		TS_ASSERT_DELTA(YPrime[1],16.0,tol);
	}
	
};



#endif //_TESTABSTRACTODESYSTEM_HPP_
