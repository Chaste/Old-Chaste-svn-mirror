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
#include "TestOdeOrderSystem.hpp"
#include "TestOdeOrderSystemOf3.hpp"
#include "TestOde4.hpp"

class TestAbstractIvpOdeSolver: public CxxTest::TestSuite
{
    public:
	
	void testEulerSolver()
	{
		TestOde1 ode_system;
		int SystemSize = 1;		
		std::vector<double> yInit(SystemSize);
	    yInit[0] = 0.0;
        
     	// Initialising the instance of our solver class	
		EulerIvpOdeSolver EulerSolver;
	    // Initialising the instance of our solution class
		OdeSolution solutions;
		
	    // Solving the ode problem and writing to solution		
	    solutions = EulerSolver.Solve(&ode_system, 0.0, 2.0, 0.001, yInit);
		
		int last = solutions.GetNumberOfTimeSteps();		
	    // Test to see if this worked		
		double testvalue = solutions.mSolutions[last][0];
		
		TS_ASSERT_DELTA(testvalue,2.0,0.01);
	}
	
	void testAdamsBashforthSolver()
	{
		TestOde1 ode_system;
		int SystemSize = 1;		
		std::vector<double> yInit(SystemSize);
		yInit[0] = 0.0;
	
     	// Initialising the instance of our solver class	
		AdamsBashforthIvpOdeSolver ABSolver;
	    // Initialising the instance of our solution class
		OdeSolution solutions;
		
	    // Solving the ode problem and writing to solution		
	    solutions = ABSolver.Solve(&ode_system, 0.0, 2.0, 0.01, yInit);
		
		int last = solutions.GetNumberOfTimeSteps();		
	    // Test to see if this worked		
		double testvalue = solutions.mSolutions[last][0];
		
		TS_ASSERT_DELTA(testvalue,2.0,0.01);
	}
	
	
	void testRungeKutta2Solver()
	{
		TestOde1 ode_system;
		int SystemSize = 1;		
		std::vector<double> yInit(SystemSize);
		yInit[0]=2.0;
	
     	// Initialising the instance of our solver class	
		RungeKutta2IvpOdeSolver RK2Solver;
	    // Initialising the instance of our solution class
		OdeSolution solutions;
		
	    // Solving the ode problem and writing to solution		
	    solutions = RK2Solver.Solve(&ode_system, 0.0, 2.0, 0.001, yInit);
		
		int last = solutions.GetNumberOfTimeSteps();		
	    // Test to see if this worked		
		double testvalue = solutions.mSolutions[last][0];
		
		TS_ASSERT_DELTA(testvalue,4.0,0.0001);
	}
	
	void testRungeKutta4Solver()
	{
		TestOde1 ode_system;
		int SystemSize = 1;		
		std::vector<double> yInit(SystemSize);
		yInit[0]=2.0;
	
     	// Initialising the instance of our solver class	
		RungeKutta4IvpOdeSolver RungeKutta4Solver;
	    // Initialising the instance of our solution class
		OdeSolution solutions;
		
	    // Solving the ode problem and writing to solution		
	    solutions = RungeKutta4Solver.Solve(&ode_system, 0.0, 2.0, 0.001, yInit);
		
		int last = solutions.GetNumberOfTimeSteps();		
	    // Test to see if this worked		
		double testvalue = solutions.mSolutions[last][0];
		
		TS_ASSERT_DELTA(testvalue,4.0,0.000001);
	}
	
	void testLastTimeStep()
	{
		TestOde1 ode_system;
		
		int SystemSize=1;
		
		std::vector<double> yInit(SystemSize);
		yInit[0] = 0.0;
	
		// Initialising the instance of our solver class
		EulerIvpOdeSolver EulerSolver;
		// Initialising the instance of our solution class
		OdeSolution solutions;
		
		// Solving the ode problem and writing to solution
	    solutions = EulerSolver.Solve(&ode_system, 0.0, 2.0, 0.000037, yInit);
		
		int last = solutions.GetNumberOfTimeSteps();
		// Test to see if this worked
		double testvalue = solutions.mSolutions[last-1][0]	;
		
		TS_ASSERT_DELTA(testvalue,2.0,0.001);
	}
	
	void testGlobalError()
	{
		TestOdeOrder ode_system;
		
		int SystemSize=1;
		
		std::vector<double> yInit(SystemSize);
		yInit[0] = 1.0;
		double hValue=0.01;
	
		//Euler solver solution worked out	
		EulerIvpOdeSolver euler_solver;
		OdeSolution solutions_euler;
		
		solutions_euler = euler_solver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last = solutions_euler.GetNumberOfTimeSteps();
		double testvalueEuler = solutions_euler.mSolutions[last][0];
		
		//Runge Kutta 2 solver solution worked out	
		RungeKutta2IvpOdeSolver RungeKutta2Solver;
		OdeSolution solutionsRungeKutta2;
		
		solutionsRungeKutta2 = RungeKutta2Solver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last2 = solutionsRungeKutta2.GetNumberOfTimeSteps();
		double testvalueRungeKutta2 = solutionsRungeKutta2.mSolutions[last2][0];
		
		//Runge Kutta 4 solver solution worked out	
		RungeKutta4IvpOdeSolver RungeKutta4Solver;
		OdeSolution solutionsRungeKutta4;
		
		solutionsRungeKutta4 = RungeKutta4Solver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last3 = solutionsRungeKutta4.GetNumberOfTimeSteps();
		double testvalueRungeKutta4 = solutionsRungeKutta4.mSolutions[last3][0];
		
		//Adams-Bashforth solver solution worked out	
		AdamsBashforthIvpOdeSolver AdamsBashforthSolver;
		OdeSolution solutionsAdamsBashforth;
		
		solutionsAdamsBashforth = AdamsBashforthSolver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last4 = solutionsAdamsBashforth.GetNumberOfTimeSteps();
		double testvalueAdamsBashforth = solutionsAdamsBashforth.mSolutions[last4][0];
		
		// The tests
		double exactSolution=exp(2);
				
		double GlobalErrorEuler;
		GlobalErrorEuler = 0.5*exp(2)*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueEuler,exactSolution,GlobalErrorEuler);
		
		double GlobalErrorRungeKutta2;
		GlobalErrorRungeKutta2 = (1.0/6.0)*hValue*exp(2)*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueRungeKutta2,exactSolution,GlobalErrorRungeKutta2);
		
		double GlobalErrorRungeKutta4;
		GlobalErrorRungeKutta4 = (1.0/24.0)*pow(hValue,3)*exp(2)*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueRungeKutta4,exactSolution,GlobalErrorRungeKutta4);
		
		double GlobalErrorAdamsBashforth;
		GlobalErrorAdamsBashforth = (1.0/6.0)*pow(hValue,3)*exp(2)*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueAdamsBashforth,exactSolution,GlobalErrorAdamsBashforth);
		
	}
		
	void testGlobalErrorSystemOf2Equations()
	{
		TestOdeOrderSystem ode_system;
		
		int SystemSize=2;
		
		std::vector<double> yInit(SystemSize);
		yInit[0] = 0.0;
		yInit[1] = 1.0;
		double hValue=0.01;
	
		//Euler solver solution worked out	
		EulerIvpOdeSolver EulerSolver;
		OdeSolution solutionsEuler;
		
		solutionsEuler = EulerSolver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last = solutionsEuler.GetNumberOfTimeSteps();
		
		double testvalueEuler[2];
		testvalueEuler[0] = solutionsEuler.mSolutions[last][0];
		testvalueEuler[1] = solutionsEuler.mSolutions[last][1];
		
		//Runge Kutta 2 solver solution worked out	
		RungeKutta2IvpOdeSolver RungeKutta2Solver;
		OdeSolution solutionsRungeKutta2;
		
		solutionsRungeKutta2 = RungeKutta2Solver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last2 = solutionsRungeKutta2.GetNumberOfTimeSteps();
		
		double testvalueRungeKutta2[2];
		testvalueRungeKutta2[0] = solutionsRungeKutta2.mSolutions[last2][0];
		testvalueRungeKutta2[1] = solutionsRungeKutta2.mSolutions[last2][1];
		
		//Runge Kutta 4 solver solution worked out	
		RungeKutta4IvpOdeSolver RungeKutta4Solver;
		OdeSolution solutionsRungeKutta4;
		
		solutionsRungeKutta4 = RungeKutta4Solver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last3 = solutionsRungeKutta4.GetNumberOfTimeSteps();
		double testvalueRungeKutta4[2];
		testvalueRungeKutta4[0] = solutionsRungeKutta4.mSolutions[last3][0];
		testvalueRungeKutta4[1] = solutionsRungeKutta4.mSolutions[last3][1];
		
		//solutionsRungeKutta4.SaveToFile("result.dat");
		
		//Adams-Bashforth solver solution worked out	
		AdamsBashforthIvpOdeSolver AdamsBashforthSolver;
		OdeSolution solutionsAdamsBashforth;
		
		solutionsAdamsBashforth = AdamsBashforthSolver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last4 = solutionsAdamsBashforth.GetNumberOfTimeSteps();
		double testvalueAdamsBashforth[2];
		testvalueAdamsBashforth[0] = solutionsAdamsBashforth.mSolutions[last4][0];
		testvalueAdamsBashforth[1] = solutionsAdamsBashforth.mSolutions[last4][1];
		
		// The tests
		double exactSolution[2];
		
		exactSolution[0] = sin(2);
		exactSolution[1] = cos(2);
				
		double GlobalErrorEuler;
		GlobalErrorEuler = 0.5*1*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueEuler[0],exactSolution[0],GlobalErrorEuler);
		TS_ASSERT_DELTA(testvalueEuler[1],exactSolution[1],GlobalErrorEuler);
		
		double GlobalErrorRungeKutta2;
		GlobalErrorRungeKutta2 = (1.0/6.0)*hValue*1*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueRungeKutta2[0],exactSolution[0],GlobalErrorRungeKutta2);
		TS_ASSERT_DELTA(testvalueRungeKutta2[1],exactSolution[1],GlobalErrorRungeKutta2);
		
		double GlobalErrorRungeKutta4;
		GlobalErrorRungeKutta4 = (1.0/24.0)*pow(hValue,3)*1*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueRungeKutta4[0],exactSolution[0],GlobalErrorRungeKutta4);
		TS_ASSERT_DELTA(testvalueRungeKutta4[1],exactSolution[1],GlobalErrorRungeKutta4);
		
		double GlobalErrorAdamsBashforth;
		GlobalErrorAdamsBashforth = (1.0/6.0)*pow(hValue,3)*1*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueAdamsBashforth[0],exactSolution[0],GlobalErrorAdamsBashforth);	
		TS_ASSERT_DELTA(testvalueAdamsBashforth[1],exactSolution[1],GlobalErrorAdamsBashforth);	
		
		
	}
	
	void testGlobalErrorSystemOf3Equations()
	{
		TestOdeOrderSystemOf3 ode_system;
		
		int SystemSize=3;
		
		std::vector<double> yInit(SystemSize);
		yInit[0] = 0.0;
		yInit[1] = 1.0;
		yInit[2] = 0.0;
		
		double hValue=0.01;
	
		//Euler solver solution worked out	
		EulerIvpOdeSolver EulerSolver;
		OdeSolution solutionsEuler;
		
		solutionsEuler = EulerSolver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last = solutionsEuler.GetNumberOfTimeSteps();
		
		double testvalueEuler[3];
		testvalueEuler[0] = solutionsEuler.mSolutions[last][0];
		testvalueEuler[1] = solutionsEuler.mSolutions[last][1];
		testvalueEuler[2] = solutionsEuler.mSolutions[last][2];
		
		//Runge Kutta 2 solver solution worked out	
		RungeKutta2IvpOdeSolver RungeKutta2Solver;
		OdeSolution solutionsRungeKutta2;
		
		solutionsRungeKutta2 = RungeKutta2Solver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last2 = solutionsRungeKutta2.GetNumberOfTimeSteps();
		
		double testvalueRungeKutta2[3];
		testvalueRungeKutta2[0] = solutionsRungeKutta2.mSolutions[last2][0];
		testvalueRungeKutta2[1] = solutionsRungeKutta2.mSolutions[last2][1];
		testvalueRungeKutta2[2] = solutionsRungeKutta2.mSolutions[last2][2];
		
		//Runge Kutta 4 solver solution worked out	
		RungeKutta4IvpOdeSolver RungeKutta4Solver;
		OdeSolution solutionsRungeKutta4;
		
		solutionsRungeKutta4 = RungeKutta4Solver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last3 = solutionsRungeKutta4.GetNumberOfTimeSteps();
		
		double testvalueRungeKutta4[3];
		testvalueRungeKutta4[0] = solutionsRungeKutta4.mSolutions[last3][0];
		testvalueRungeKutta4[1] = solutionsRungeKutta4.mSolutions[last3][1];
		testvalueRungeKutta4[2] = solutionsRungeKutta4.mSolutions[last3][2];
		
		//solutionsRungeKutta4.SaveToFile("result.dat");
		
		//Adams-Bashforth solver solution worked out	
		AdamsBashforthIvpOdeSolver AdamsBashforthSolver;
		OdeSolution solutionsAdamsBashforth;
		
		solutionsAdamsBashforth = AdamsBashforthSolver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last4 = solutionsAdamsBashforth.GetNumberOfTimeSteps();
		
		double testvalueAdamsBashforth[3];
		testvalueAdamsBashforth[0] = solutionsAdamsBashforth.mSolutions[last4][0];
		testvalueAdamsBashforth[1] = solutionsAdamsBashforth.mSolutions[last4][1];
		testvalueAdamsBashforth[2] = solutionsAdamsBashforth.mSolutions[last4][2];
		
		// The tests
		double exactSolution[3];
		
		exactSolution[0] = -sin(2);
		exactSolution[1] = sin(2)+cos(2);
		exactSolution[2] = 2*sin(2);
				
		double GlobalErrorEuler;
		GlobalErrorEuler = 0.5*2*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueEuler[0],exactSolution[0],GlobalErrorEuler);
		TS_ASSERT_DELTA(testvalueEuler[1],exactSolution[1],GlobalErrorEuler);
		TS_ASSERT_DELTA(testvalueEuler[2],exactSolution[2],GlobalErrorEuler);
		
		double GlobalErrorRungeKutta2;
		GlobalErrorRungeKutta2 = (1.0/6.0)*hValue*2*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueRungeKutta2[0],exactSolution[0],GlobalErrorRungeKutta2);
		TS_ASSERT_DELTA(testvalueRungeKutta2[1],exactSolution[1],GlobalErrorRungeKutta2);
		TS_ASSERT_DELTA(testvalueRungeKutta2[2],exactSolution[2],GlobalErrorRungeKutta2);
		
		double GlobalErrorRungeKutta4;
		GlobalErrorRungeKutta4 = (1.0/24.0)*pow(hValue,3)*2*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueRungeKutta4[0],exactSolution[0],GlobalErrorRungeKutta4);
		TS_ASSERT_DELTA(testvalueRungeKutta4[1],exactSolution[1],GlobalErrorRungeKutta4);
		TS_ASSERT_DELTA(testvalueRungeKutta4[2],exactSolution[2],GlobalErrorRungeKutta4);
		
		double GlobalErrorAdamsBashforth;
		GlobalErrorAdamsBashforth = (1.0/6.0)*pow(hValue,3)*2*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueAdamsBashforth[0],exactSolution[0],GlobalErrorAdamsBashforth);	
		TS_ASSERT_DELTA(testvalueAdamsBashforth[1],exactSolution[1],GlobalErrorAdamsBashforth);	
		TS_ASSERT_DELTA(testvalueAdamsBashforth[2],exactSolution[2],GlobalErrorAdamsBashforth);
		
	}
	
	void testGlobalError2()
	{
		TestOde4 ode_system;
		
		int SystemSize=1;
		
		std::vector<double> yInit(SystemSize);
		yInit[0] = 0.5;
		double hValue=0.001;
	
		//Euler solver solution worked out	
		EulerIvpOdeSolver EulerSolver;
		OdeSolution solutionsEuler;
		
		solutionsEuler = EulerSolver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last = solutionsEuler.GetNumberOfTimeSteps();
		double testvalueEuler = solutionsEuler.mSolutions[last][0];
		
		//Runge Kutta 2 solver solution worked out	
		RungeKutta2IvpOdeSolver RungeKutta2Solver;
		OdeSolution solutionsRungeKutta2;
		
		solutionsRungeKutta2 = RungeKutta2Solver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last2 = solutionsRungeKutta2.GetNumberOfTimeSteps();
		double testvalueRungeKutta2 = solutionsRungeKutta2.mSolutions[last2][0];
		
		//Runge Kutta 4 solver solution worked out	
		RungeKutta4IvpOdeSolver RungeKutta4Solver;
		OdeSolution solutionsRungeKutta4;
		
		solutionsRungeKutta4 = RungeKutta4Solver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last3 = solutionsRungeKutta4.GetNumberOfTimeSteps();
		double testvalueRungeKutta4 = solutionsRungeKutta4.mSolutions[last3][0];
		
		//Adams-Bashforth solver solution worked out	
		AdamsBashforthIvpOdeSolver AdamsBashforthSolver;
		OdeSolution solutionsAdamsBashforth;
		
		solutionsAdamsBashforth = AdamsBashforthSolver.Solve(&ode_system, 0.0, 2.0, hValue, yInit);
		int last4 = solutionsAdamsBashforth.GetNumberOfTimeSteps();
		double testvalueAdamsBashforth = solutionsAdamsBashforth.mSolutions[last4][0];
		
		// The tests
		double alpha = 100;
		double exactSolution=1/(1+exp(-alpha*2));
				
		double GlobalErrorEuler;
		GlobalErrorEuler = 0.5*1/(1+exp(-alpha*2))*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueEuler,exactSolution,GlobalErrorEuler);
		
		double GlobalErrorRungeKutta2;
		GlobalErrorRungeKutta2 = (1.0/6.0)*hValue*1/(1+exp(-alpha*2))*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueRungeKutta2,exactSolution,GlobalErrorRungeKutta2);
		
		double GlobalErrorRungeKutta4;
		GlobalErrorRungeKutta4 = (1.0/24.0)*pow(hValue,3)*1/(1+exp(-alpha*2))*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueRungeKutta4,exactSolution,GlobalErrorRungeKutta4);
		
		double GlobalErrorAdamsBashforth;
		GlobalErrorAdamsBashforth = (1.0/6.0)*pow(hValue,3)*1/(1+exp(-alpha*2))*(exp(2)-1)*hValue;
		TS_ASSERT_DELTA(testvalueAdamsBashforth,exactSolution,GlobalErrorAdamsBashforth);
	}
};

#endif //_TESTABSTRACTIVPODESOLVER_HPP_
