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
#include "TestOdeOrder.hpp"

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
	
	void testGlobalError()
	{
		TestOdeOrder* pMyOdeSystem = new TestOdeOrder();
		
		int SystemSize=1;
		
		std::vector<double> yInit(SystemSize);
		yInit[0] = 1.0;
		double hValue=0.01;
	
		//Euler solver solution worked out	
		EulerIvpOdeSolver* myEulerSolver = new EulerIvpOdeSolver;
		OdeSolution solutionsEuler;
		
		solutionsEuler = myEulerSolver->Solve(pMyOdeSystem, 0.0, 2.0, hValue, yInit);
		int last = solutionsEuler.mNumberOfTimeSteps;
		double testvalueEuler = solutionsEuler.mSolutions[last][0]	;
		
		//Runge Kutta 2 solver solution worked out	
		RungeKutta2IvpOdeSolver* myRungeKutta2Solver = new RungeKutta2IvpOdeSolver;
		OdeSolution solutionsRungeKutta2;
		
		solutionsRungeKutta2 = myRungeKutta2Solver->Solve(pMyOdeSystem, 0.0, 2.0, hValue, yInit);
		int last2 = solutionsRungeKutta2.mNumberOfTimeSteps;
		double testvalueRungeKutta2 = solutionsRungeKutta2.mSolutions[last2][0]	;
		
		//Runge Kutta 4 solver solution worked out	
		RungeKutta4IvpOdeSolver* myRungeKutta4Solver = new RungeKutta4IvpOdeSolver;
		OdeSolution solutionsRungeKutta4;
		
		solutionsRungeKutta4 = myRungeKutta4Solver->Solve(pMyOdeSystem, 0.0, 2.0, hValue, yInit);
		int last3 = solutionsRungeKutta4.mNumberOfTimeSteps;
		double testvalueRungeKutta4 = solutionsRungeKutta4.mSolutions[last2][0]	;
		
		// The tests
		double GlobalErrorEuler, workedOutSolnEuler=exp(2);
		GlobalErrorEuler=.5*exp(2)*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueEuler,workedOutSolnEuler,GlobalErrorEuler);
		
		double GlobalErrorRungeKutta2, workedOutSolnRungeKutta2=exp(2);
		GlobalErrorRungeKutta2=(1.0/6.0)*hValue*exp(2)*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueRungeKutta2,workedOutSolnRungeKutta2,GlobalErrorRungeKutta2);
		
		double GlobalErrorRungeKutta4, workedOutSolnRungeKutta4=exp(2);
		GlobalErrorRungeKutta4=(1.0/24.0)*pow(hValue,3)*exp(2)*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueRungeKutta4,workedOutSolnRungeKutta4,GlobalErrorRungeKutta4);
	}
	
	
	
	
};

#endif //_TESTABSTRACTIVPODESOLVER_HPP_
