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
		std::cout << "\n Start of OdeSystem tests\n";
		
		TS_ASSERT(1+1 >1)
	}
	
		
	void TestOdeOne(void)
	{
		std::vector<double> yInit(1);
		yInit[0]=0.0;
		
		// pointer to TestOde1 class		
		TestOde1 * pOde1 = new TestOde1();
		// Yprime
		std::vector<double> rYPrime2;
		rYPrime2=pOde1->EvaluateYDerivatives(1.0, yInit);
		std::cout << rYPrime2[0];
		TS_ASSERT_DELTA(rYPrime2[0],1.0,tol);
		std::cout << "\n End of third set of tests\n";
	}
	
	
	void TestOdeTwo(void)
	{
		std::cout << "Start of fourth test\n";
		std::vector<double> yInit(1);
		
		yInit[0]=4.0;
		
		TestOde2 * pOde2 = new TestOde2();
		std::vector<double> rYPrime2;
		rYPrime2=pOde2->EvaluateYDerivatives(2.0, yInit);
		std::cout << rYPrime2[0];
		TS_ASSERT_DELTA(rYPrime2[0],8.0,tol);
		std::cout << "\n End of fourth set of tests\n";
	}
	
	void TestOdeThree(void)
	{
		std::cout << "Start of fifth test\n";
		std::vector<double> yInit(2);
		
		yInit[0]=4.0;
		yInit[1]=8.0;
		
		TestOde3 * pOde3 = new TestOde3();
		std::vector<double> rYPrime3;
		rYPrime3=pOde3->EvaluateYDerivatives(2.0, yInit);
		std::cout << rYPrime3[0] << std::endl;
		std::cout << rYPrime3[1] << std::endl;
		TS_ASSERT_DELTA(rYPrime3[0],8.0,tol);
		TS_ASSERT_DELTA(rYPrime3[1],16.0,tol);
		std::cout << "\n End of fifth set of tests\n";
		
		
	std::cout << "\n End of OdeSystem tests\n";
	}
	
};



#endif //_TESTABSTRACTODESYSTEM_HPP_
