#ifndef _TESTABSTRACTIVPODESOLVER_HPP_
#define _TESTABSTRACTIVPODESOLVER_HPP_
#include <cxxtest/TestSuite.h>

#include <vector>
#include <iostream>

#include "AbstractIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "AdamsBashforthIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
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
	
	void testAdamsBashforthSolver()
	{
		TestOde1* pMyOdeSystem = new TestOde1();
		int SystemSize = 1;		
		std::vector<double> yInit(SystemSize);
		yInit[0] = 0.0;
	
     	// Initialising the instance of our solver class	
		AdamsBashforthIvpOdeSolver* myABSolver = new AdamsBashforthIvpOdeSolver;
	    // Initialising the instance of our solution class
		OdeSolution solutions;
		
	    // Solving the ode problem and writing to solution		
	    solutions = myABSolver->Solve(pMyOdeSystem, 0.0, 2.0, 0.01, yInit);
		
		int last = solutions.mNumberOfTimeSteps;		
	    // Test to see if this worked		
		double testvalue = solutions.mSolutions[last-1][0];
		
		TS_ASSERT_DELTA(testvalue,2.0,0.01);
	}
	
	
	void testRungeKutta2Solver()
	{
		TestOde1* pMyOdeSystem = new TestOde1();
		int SystemSize = 1;		
		std::vector<double> yInit(SystemSize);
		yInit[0]=2.0;
	
     	// Initialising the instance of our solver class	
		RungeKutta2IvpOdeSolver* myRK2Solver = new RungeKutta2IvpOdeSolver;
	    // Initialising the instance of our solution class
		OdeSolution solutions;
		
	    // Solving the ode problem and writing to solution		
	    solutions = myRK2Solver->Solve(pMyOdeSystem, 0.0, 2.0, 0.001, yInit);
		
		int last = solutions.mNumberOfTimeSteps;		
	    // Test to see if this worked		
		double testvalue = solutions.mSolutions[last][0];
		
		TS_ASSERT_DELTA(testvalue,4.0,0.0001);
	}
	
	void testRungeKutta4Solver()
	{
		TestOde1* pMyOdeSystem = new TestOde1();
		int SystemSize = 1;		
		std::vector<double> yInit(SystemSize);
		yInit[0]=2.0;
	
     	// Initialising the instance of our solver class	
		RungeKutta4IvpOdeSolver* myRungeKutta4Solver = new RungeKutta4IvpOdeSolver;
	    // Initialising the instance of our solution class
		OdeSolution solutions;
		
	    // Solving the ode problem and writing to solution		
	    solutions = myRungeKutta4Solver->Solve(pMyOdeSystem, 0.0, 2.0, 0.001, yInit);
		
		int last = solutions.mNumberOfTimeSteps;		
	    // Test to see if this worked		
		double testvalue = solutions.mSolutions[last][0];
		
		TS_ASSERT_DELTA(testvalue,4.0,0.000001);
	}
	
	void testLastTimeStep()
	{
		TestOde1* pMyOdeSystem = new TestOde1();
		
		int SystemSize=1;
		
		std::vector<double> yInit(SystemSize);
		yInit[0] = 0.0;
	
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
