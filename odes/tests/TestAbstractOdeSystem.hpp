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
		
	void testAddition(void)
	{
		TS_ASSERT( 1 + 1 > 1 );
	}
		
	void TestOdeSystemOne(void)
	{
		std::vector<double> yInit(1);
		yInit[0]=0.0;
		
		// pointer to TestOde1 class		
		TestOde1* pOde1 = new TestOde1();
		// Yprime
		std::vector<double> YPrime;
		YPrime = pOde1->EvaluateYDerivatives(1.0, yInit);
		TS_ASSERT_DELTA(YPrime[0],1.0,tol);
	}
	
	
	void TestOdeSystemTwo(void)
	{
		std::vector<double> yInit(1);
		yInit[0] = 4.0;
		TestOde2* pOde2 = new TestOde2();
		std::vector<double> YPrime;
		YPrime = pOde2->EvaluateYDerivatives(2.0, yInit);
		TS_ASSERT_DELTA(YPrime[0],8.0,tol);
	}
	
	void TestOdeSystemThree(void)
	{
		std::vector<double> yInit(2);
		yInit[0] = 4.0;
		yInit[1] = 8.0;

		TestOde3* pOde3 = new TestOde3();
		std::vector<double> YPrime;
		YPrime = pOde3->EvaluateYDerivatives(2.0, yInit);
		TS_ASSERT_DELTA(YPrime[0],8.0,tol);
		TS_ASSERT_DELTA(YPrime[1],16.0,tol);
	}
	
};



#endif //_TESTABSTRACTODESYSTEM_HPP_
