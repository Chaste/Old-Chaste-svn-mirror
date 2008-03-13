#ifndef TESTWNTCELLCYCLEODESYSTEM_HPP_
#define TESTWNTCELLCYCLEODESYSTEM_HPP_

#include <stdio.h>
#include <time.h>
#include <cxxtest/TestSuite.h>
#include "Lee2003WntSignallingOdeSystem.hpp"
#include <vector>
#include <iostream>
#include "RungeKutta4IvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestLee2003WntSignallingOdeSystem : public CxxTest::TestSuite
{
public:

    void TestLee2003WntSignallingEquations()
    {
        double WntLevel = 0.0;
        Lee2003WntSignallingOdeSystem lee_system(WntLevel);
        
        double time = 0.0;
        std::vector<double> initial_conditions = lee_system.GetInitialConditions();
        
        std::vector<double> derivs(initial_conditions.size());
        lee_system.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        double factor = 60.0; // to convert to hours
        TS_ASSERT_DELTA(derivs[0],factor*0, 1e-5);
        TS_ASSERT_DELTA(derivs[1],factor*-5.744861643948267e-06, 1e-5);
        TS_ASSERT_DELTA(derivs[2],factor*1.132999999999829e-04, 1e-5);
        TS_ASSERT_DELTA(derivs[3],factor*-8.799999999999364e-04, 1e-5);
        TS_ASSERT_DELTA(derivs[4],factor*-8.972455467983818e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[5],factor*2.725552844639836e-04, 1e-5);
        TS_ASSERT_DELTA(derivs[6],factor*2.649781064090187e-07, 1e-5);
        TS_ASSERT_DELTA(derivs[7],factor*0.0, 1e-5);
    }
    
    void TestLee2003WntSignallingOdeSolver() throw(Exception)
    {
        double WntLevel = 1.0;
        Lee2003WntSignallingOdeSystem lee_system(WntLevel);
        
        // Solve system using rk4 solver
        
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.
        
        double h_value = 0.0001;
        
        RungeKutta4IvpOdeSolver rk4_solver;
        
        OdeSolution solutions;
        
        std::vector<double> initial_conditions = lee_system.GetInitialConditions();
                
        double start_time, end_time, elapsed_time = 0.0;
        start_time = std::clock();
        solutions = rk4_solver.Solve(&lee_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "1. Runge-Kutta Elapsed time = " << elapsed_time << "\n";
        

        
        // Test solutions are correct for a new steady state
        int end = solutions.rGetSolutions().size() - 1;
        
        // Test the simulation is ending at the right time
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 100.0, 1e-2);
        
        // Proper values calculated using the MatLab stiff ODE solver ode15s. Note that 
        // large tolerances are required for the tests to pass (see #238 and #316).
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],9.090909090909091e+01, 1e-5);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],7.275154952501657e-04, 1e-5);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],1.862484031071281e-03, 1e-5);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],9.200760441262037e-01, 1e-5);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],1.460501031816546e-03, 1e-5);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],1.530283641433451e+02, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6],4.922155688613715e-04, 1e-5);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7],1.00, 1e-5);
    }
};

#endif /*TESTWNTCELLCYCLEODESYSTEM_HPP_*/
