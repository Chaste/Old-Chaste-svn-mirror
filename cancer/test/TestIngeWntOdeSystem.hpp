#ifndef TESTINGEWNTODESYSTEM_HPP_
#define TESTINGEWNTODESYSTEM_HPP_

#include <stdio.h>
#include <time.h>
#include <cxxtest/TestSuite.h>
#include "IngeWntOdeSystem.hpp"
#include <vector>
#include <iostream>
#include "RungeKutta4IvpOdeSolver.hpp"
#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "ColumnDataWriter.hpp"
#include "CryptCellMutationStates.hpp"


class TestIngeWntOdeSystem : public CxxTest::TestSuite
{
public:

    void TestIngeWntOdeEquations()
    {
        double WntLevel = 0.0;
        IngeWntOdeSystem wnt_ode_system(WntLevel);
        
        TS_ASSERT_EQUALS(wnt_ode_system.rGetMutationState(), HEALTHY);
        
        double time = 0.0;
        std::vector<double> initial_conditions = wnt_ode_system.GetInitialConditions();
        
        TS_ASSERT_EQUALS(initial_conditions.size(), 17u);
        TS_ASSERT_DELTA(initial_conditions[0],2.0/3.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[1],2.0/30.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[2],0.4492, 1e-4);
        TS_ASSERT_DELTA(initial_conditions[3],2.5403, 1e-4);
        TS_ASSERT_DELTA(initial_conditions[4],0.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[5],0.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[6],0.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[7],10.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[8],18.1449, 1e-4);
        TS_ASSERT_DELTA(initial_conditions[9],0.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[10],25.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[11],2.5403, 1e-4);
        TS_ASSERT_DELTA(initial_conditions[12],0.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[13],0.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[14],0.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[15],0.4835, 1e-4);
        TS_ASSERT_DELTA(initial_conditions[16],WntLevel, 1e-5);
        
        wnt_ode_system.SetUseHypothesisTwo(false);
                
        std::vector<double> derivs(initial_conditions.size());
        wnt_ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from Mathematica code)
        for (unsigned i=0 ; i<17 ; i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-5);
        }

        
        wnt_ode_system.SetUseHypothesisTwo(true);
        wnt_ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        // Test derivatives are correct at t=0 for these initial conditions
        for (unsigned i=0 ; i<17 ; i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-5);
        }
        
        /**
         * And the same for a high Wnt level
         */
        WntLevel = 0.0;
        IngeWntOdeSystem wnt_ode_system2(WntLevel,LABELLED);
        initial_conditions = wnt_ode_system2.GetInitialConditions();
        initial_conditions[16] = 1.0;
        wnt_ode_system2.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from MatLab code)
        TS_ASSERT_DELTA(derivs[0],-6.666666667, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-10.0, 1e-5);
        for (unsigned i=2 ; i<17 ; i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-5);
        }
        
        wnt_ode_system2.SetUseHypothesisTwo(true);
        wnt_ode_system2.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        // Test derivatives are correct at t=0 for these initial conditions

        TS_ASSERT_DELTA(derivs[0],-6.666666667, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-10.0, 1e-5);
        TS_ASSERT_DELTA(derivs[2],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[3],-62.7107, 1e-4);
        TS_ASSERT_DELTA(derivs[4],62.7107, 1e-4);
        for (unsigned i=5 ; i<17 ; i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-5);
        }
        
        
        /**
         * A test for the APC +/- mutation
         */
        CryptCellMutationState mutation = HEALTHY;
        WntLevel = 0.0;
        IngeWntOdeSystem wnt_ode_system3(WntLevel,mutation);
        initial_conditions = wnt_ode_system3.GetInitialConditions();
        
        wnt_ode_system3.SetMutationState(APC_ONE_HIT);
        
        wnt_ode_system3.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from MatLab code)
        TS_ASSERT_DELTA(derivs[0],-3.333333334, 1e-5);
        TS_ASSERT_DELTA(derivs[1],3.333333334, 1e-5);
        for (unsigned i=2 ; i<17 ; i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-5);
        }
        
        /**
        * A test for the beta-cat delta45 mutation
        */
        mutation = BETA_CATENIN_ONE_HIT;
        WntLevel = 0.0;
        IngeWntOdeSystem wnt_ode_system4(WntLevel,mutation);
        initial_conditions = wnt_ode_system4.GetInitialConditions();
        
        wnt_ode_system4.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from MatLab code)
        TS_ASSERT_DELTA(derivs[0],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[1],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[2],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[3],-12.5, 1e-5);
        TS_ASSERT_DELTA(derivs[4],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[5],12.5, 1e-5);
        for (unsigned i=6 ; i<17 ; i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-5);
        }
        
        /**
        * A test for APC -/- mutation
        */
        mutation = HEALTHY;
        WntLevel = 0.0;
        IngeWntOdeSystem wnt_ode_system5(WntLevel,mutation);
        mutation = APC_TWO_HIT;
        wnt_ode_system5.SetMutationState(APC_TWO_HIT);
        initial_conditions = wnt_ode_system5.GetInitialConditions();
        //std::cout << "mutation " << mutation << " beta-cat = " << initial_conditions[6] << "\n";
                
        wnt_ode_system5.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from MatLab code)
        TS_ASSERT_DELTA(derivs[0],-6.666666666667, 1e-5);
        TS_ASSERT_DELTA(derivs[1],6.66666666667, 1e-5);
        for (unsigned i=2 ; i<17 ; i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-5);
        }
    }
    
    void TestWntCellCycleSolver() throw(Exception)
    {
        double WntLevel = 0.0;
        IngeWntOdeSystem wnt_system(WntLevel,LABELLED);
        
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.
        
        double h_value=0.0001;
        
        RungeKutta4IvpOdeSolver rk4_solver;
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        
        OdeSolution solutions;
        //OdeSolution solutions2;
        
        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();
        // Set Wnt to one at time zero...
        initial_conditions[16] = 1.0;
        
        //double start_time, end_time, elapsed_time = 0.0;
        //start_time = std::clock();
        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 8.0, h_value, h_value);
        //end_time = std::clock();
        //elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        //std::cout <<  "1. Runge-Kutta Elapsed time = " << elapsed_time << "\n";
        
        h_value = 0.1;
        initial_conditions = wnt_system.GetInitialConditions();
        // Set Wnt to one at time zero...
        initial_conditions[16] = 1.0;
        
        //start_time = std::clock();
        solutions = rkf_solver.Solve(&wnt_system, initial_conditions, 0.0, 8.0, h_value, 1e-4);
        //end_time = std::clock();
        //elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        //std::cout <<  "2. Runge-Kutta-Fehlberg Elapsed time = " << elapsed_time << "\n";
        
        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Proper values from MatLab ode15s - shocking tolerances to pass though.
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],0.1428, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],0.0286, 1.01e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],0.1959, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],10.907, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6],0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7],8.1325, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8],63.3473, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][9],0.0, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][10],23.5696, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][11],10.2824, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][12],0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][13],0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][14],0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][15],1.6411, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][16],1.00, 1e-3);
    }
    
    void TestWntCellCycleSolverHypothesisTwo() throw(Exception)
    {
        double WntLevel = 0.0;
        IngeWntOdeSystem wnt_system(WntLevel,LABELLED);
        wnt_system.SetUseHypothesisTwo(true);
        
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.
        
        double h_value=0.0001;
        
        RungeKutta4IvpOdeSolver rk4_solver;
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        
        OdeSolution solutions;
        //OdeSolution solutions2;
        
        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();
        // Set Wnt to one at time zero...
        initial_conditions[16] = 1.0;
        
        //double start_time, end_time, elapsed_time = 0.0;
        //start_time = std::clock();
        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 8.0, h_value, h_value);
        //end_time = std::clock();
        //elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        //std::cout <<  "1. Runge-Kutta Elapsed time = " << elapsed_time << "\n";
        
//        h_value = 0.0001;
//        initial_conditions = wnt_system.GetInitialConditions();
//        // Set Wnt to one at time zero...
//        initial_conditions[16] = 1.0;
//        
//        start_time = std::clock();
//        solutions = rkf_solver.Solve(&wnt_system, initial_conditions, 0.0, 8.0, h_value, 1e-3);
//        end_time = std::clock();
//        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
//        std::cout <<  "2. Runge-Kutta-Fehlberg Elapsed time = " << elapsed_time << "\n";
        
        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Proper values from MatLab ode15s - shocking tolerances to pass though.
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],0.1428, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],0.0286, 1.01e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],0.2113, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],0.9385, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],13.2706, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6],0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7],10.0016, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8],6.705, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][9],0.0, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][10],23.8358, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][11],0.8948, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][12],12.6524, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][13],0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][14],0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][15],2.1008, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][16],1.00, 1e-3);
    }
    
    
    
};

#endif /*TESTINGEWNTODESYSTEM_HPP_*/
