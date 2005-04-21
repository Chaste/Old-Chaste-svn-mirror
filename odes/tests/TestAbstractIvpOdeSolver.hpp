#ifndef _TESTABSTRACTIVPODESOLVER_HPP_
#define _TESTABSTRACTIVPODESOLVER_HPP_
#include <cxxtest/TestSuite.h>

#include <vector>

// The superb classes for our assignment
#include "AbstractIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "TestOde1.hpp"

class TestAbstractIvpOdeSolver: public CxxTest::TestSuite
{
    public:
    
    void testAddition( void )
    {
    	std::cout<<"Start of OdeSolver Tests"<<std::endl;
        //TS_TRACE("Running our test ****************");
        TS_ASSERT( 1 + 1 > 1 );
    }
    
    

	// Test we can construct an ODESolver
	void testconstructOdeSolver()
	{
		TS_ASSERT_THROWS_NOTHING(EulerIvpOdeSolver anOdeSolver);		
	}
	
	
	void testIvpOdeSolver()
	{
		TestOde1* pMyOdeSystem = new TestOde1();
		
		int SystemSize=1;
		
		std::vector<double> yInit(SystemSize);
	
	// Initialising the instance of our solver class	
		EulerIvpOdeSolver* myEulerSolver = new EulerIvpOdeSolver;
	// Initialising the instance of our solution class
		OdeSolution solutions;
		
	// Solving the ode problem and writing to solution		
	    solutions = myEulerSolver->Solve(pMyOdeSystem, 0.0, 2.0, 0.001, yInit);
		
		int last = solutions.mNumberOfTimeSteps;
		
		
	// Test to see if this worked		
		double testvalue = solutions.mSolutions[last-1][0]	;
		
		TS_ASSERT_DELTA(testvalue,2.0,0.01);
	}
	
	
};

#endif //_TESTABSTRACTIVPODESOLVER_HPP_
