#ifndef _TESTABSTRACTIVPODESOLVER_HPP_
#define _TESTABSTRACTIVPODESOLVER_HPP_
#include <cxxtest/TestSuite.h>

#include <vector>

#include "AbstractIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RK2IvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "TestOde1.hpp"

class TestAbstractIvpOdeSolver: public CxxTest::TestSuite
{
    public:
    
    void testAddition( void )
    {
        TS_ASSERT( 1 + 1 > 1 );
    }	
	
	void testEulerSolver()
	{
		TestOde1* pMyOdeSystem = new TestOde1();
		int SystemSize = 1;		
		std::vector<double> yInit(SystemSize);
	
     	// Initialising the instance of our solver class	
		EulerIvpOdeSolver* myEulerSolver = new EulerIvpOdeSolver;
	    // Initialising the instance of our solution class
		OdeSolution solutions;
		
	    // Solving the ode problem and writing to solution		
	    solutions = myEulerSolver->Solve(pMyOdeSystem, 0.0, 2.0, 0.001, yInit);
		
		int last = solutions.mNumberOfTimeSteps;		
	    // Test to see if this worked		
		double testvalue = solutions.mSolutions[last-1][0];
		
		TS_ASSERT_DELTA(testvalue,2.0,0.01);
	}
	
	void testRK2Solver()
	{
		TestOde1* pMyOdeSystem = new TestOde1();
		int SystemSize = 1;		
		std::vector<double> yInit(SystemSize);
		yInit[0]=2.0;
	
     	// Initialising the instance of our solver class	
		RK2IvpOdeSolver* myRK2Solver = new RK2IvpOdeSolver;
	    // Initialising the instance of our solution class
		OdeSolution solutions;
		
	    // Solving the ode problem and writing to solution		
	    solutions = myRK2Solver->Solve(pMyOdeSystem, 0.0, 2.0, 0.001, yInit);
		
		int last = solutions.mNumberOfTimeSteps;		
	    // Test to see if this worked		
		double testvalue = solutions.mSolutions[last-1][0];
		
		TS_ASSERT_DELTA(testvalue,4.0,0.01);
		std::cout << "RK2 test result = " << testvalue << std::endl;
	}
	
	
	void testLastTimeStep()
	{
		TestOde1* pMyOdeSystem = new TestOde1();
		
		int SystemSize=1;
		
		std::vector<double> yInit(SystemSize);
	
	// Initialising the instance of our solver class	
		EulerIvpOdeSolver* myEulerSolver = new EulerIvpOdeSolver;
	// Initialising the instance of our solution class
		OdeSolution solutions;
		
	// Solving the ode problem and writing to solution		
	    solutions = myEulerSolver->Solve(pMyOdeSystem, 0.0, 2.0, 0.000037, yInit);
		
		int last = solutions.mNumberOfTimeSteps;
	// Test to see if this worked		
		double testvalue = solutions.mSolutions[last-1][0]	;
		
		TS_ASSERT_DELTA(testvalue,2.0,0.001);
	}
	
	
};

#endif //_TESTABSTRACTIVPODESOLVER_HPP_
