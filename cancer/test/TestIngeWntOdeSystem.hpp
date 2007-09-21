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
        TS_ASSERT_DELTA(initial_conditions[2],0.4492, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[3],2.5403, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[4],0.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[5],0.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[6],0.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[7],10.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[8],18.1449, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[9],18.1449, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[10],25.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[11],2.5403, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[12],0.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[13],0.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[14],0.0, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[15],0.4835, 1e-5);
        TS_ASSERT_DELTA(initial_conditions[16],WntLevel, 1e-5);
        
        wnt_ode_system.SetUseHypothesisTwo(false);
                
        std::vector<double> derivs(initial_conditions.size());
        wnt_ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from MatLab code)
        for (unsigned i=0 ; i<17 ; i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-5);
        }
        TS_ASSERT_DELTA(derivs[0],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[1],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[2],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[3],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[4],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[5],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[9],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[10],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[11],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[12],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[13],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[14],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[15],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[16],0.0, 1e-5);
        
        wnt_ode_system.SetUseHypothesisTwo(true);
        wnt_ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        // Test derivatives are correct at t=0 for these initial conditions
        for (unsigned i=0 ; i<17 ; i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-5);
        }
        TS_ASSERT_DELTA(derivs[0],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[1],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[2],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[3],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[4],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[5],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[9],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[10],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[11],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[12],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[13],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[14],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[15],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[16],0.0, 1e-5);
        
        /**
         * And the same for a high Wnt level
         */
        WntLevel = 1.0;
        IngeWntOdeSystem wnt_ode_system2(WntLevel,LABELLED);
        initial_conditions = wnt_ode_system2.GetInitialConditions();
        //std::cout << "mutation 0 beta-cat = " << initial_conditions[6] << "\n";
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
        TS_ASSERT_DELTA(derivs[3],-62.7108, 1e-5);
        TS_ASSERT_DELTA(derivs[4],62.7108, 1e-5);
        for (unsigned i=5 ; i<17 ; i++)
        {
            TS_ASSERT_DELTA(derivs[i],0.0, 1e-5);
        }
        
        /**
         * A test for the APC +/- mutation
         */
        CryptCellMutationState mutation = APC_ONE_HIT;
        WntLevel = 0.0;
        IngeWntOdeSystem wnt_ode_system3(WntLevel,mutation);
        initial_conditions = wnt_ode_system3.GetInitialConditions();
        //std::cout << "mutation " << mutation << " beta-cat = " << initial_conditions[6] << "\n";
        TS_ASSERT_DELTA(initial_conditions[6]+initial_conditions[7],0.750207,1e-6);
        
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
        mutation = APC_TWO_HIT;
        WntLevel = 0.0;
        IngeWntOdeSystem wnt_ode_system5(WntLevel,mutation);
        initial_conditions = wnt_ode_system5.GetInitialConditions();
        //std::cout << "mutation " << mutation << " beta-cat = " << initial_conditions[6] << "\n";
        
        TS_ASSERT_DELTA(initial_conditions[6]+initial_conditions[7],1.0,1e-6);
        
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
    
    void xTestWntCellCycleSolver() throw(Exception)
    {
        double WntLevel = 1.0;
        IngeWntOdeSystem wnt_system(WntLevel,LABELLED);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.
        
        double h_value=0.001;
        
        RungeKutta4IvpOdeSolver rk4_solver;
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        BackwardEulerIvpOdeSolver back_solver(9);
        
        OdeSolution solutions;
        //OdeSolution solutions2;
        
        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();
                
        double start_time, end_time, elapsed_time = 0.0;
        start_time = std::clock();
        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "1. Runge-Kutta Elapsed time = " << elapsed_time << "\n";
        
        h_value = 0.1;
        
        initial_conditions = wnt_system.GetInitialConditions();
        start_time = std::clock();
        solutions = rkf_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, 1e-4);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "2. Runge-Kutta-Fehlberg Elapsed time = " << elapsed_time << "\n";
        
//        IngeWntOdeSystem wnt_system_2(WntLevel);
//        initial_conditions = wnt_system.GetInitialConditions();
//
//        h_value = 0.001;
//
//        start_time = std::clock();
//        solutions = back_solver.Solve(&wnt_system_2, initial_conditions, 0.0, 100.0, h_value, h_value);
//        end_time = std::clock();
//        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
//        std::cout <<  "1. BackwardEuler Elapsed time = " << elapsed_time << "\n";

        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 5.971 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end] , 5.971 , 1e-2);
        // Proper values from MatLab ode15s - shocking tolerances to pass though.
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],2.880603485931000e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],1.000220438771564e+00, 1.01e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],2.453870380958196e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],1.446185835615586e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],1.383272155041549e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],4.975124378109454e-03, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6]+solutions.rGetSolutions()[end][7],6.002649406788524e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8],1.00, 1e-3);
    }
    
    void xTestWntCellCycleSolverWithAPCSingleHit() throw(Exception)
    {
        double WntLevel = 1.0;
        IngeWntOdeSystem wnt_system(WntLevel,APC_ONE_HIT);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.
        
        double h_value=0.001;
        
        RungeKutta4IvpOdeSolver rk4_solver;
        BackwardEulerIvpOdeSolver back_solver(9);
        
        OdeSolution solutions;
        //OdeSolution solutions2;
        
        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();
                
        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        
        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 3.94 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end] , 4.804 , 1e-2);
        // Proper values from MatLab ode15s - shocking tolerances to pass though.
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],2.493601889546602e-01, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],1.0, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],2.922278616458120e+00, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],1.913484075138688e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],1.583557222931072e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],4.975124378109454e-03, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6],0.375, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7],0.375, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8],1.00, 1e-3);
    }
    
    void xTestWntCellCycleSolverWithBetaCateninHit() throw(Exception)
    {
        double WntLevel = 0.0;
        IngeWntOdeSystem wnt_system(WntLevel,BETA_CATENIN_ONE_HIT);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.
        
        double h_value=0.001;
        
        RungeKutta4IvpOdeSolver rk4_solver;
        BackwardEulerIvpOdeSolver back_solver(9);
        
        OdeSolution solutions;
        //OdeSolution solutions2;
        
        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();
        
        //double start_time, end_time, elapsed_time = 0.0;
        //start_time = std::clock();
        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        //end_time = std::clock();
        //elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        //std::cout <<  "1. Runge-Kutta Elapsed time = " << elapsed_time << "\n";
        

        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 7.81 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end] , 7.82 , 1e-2);
        // Proper values from MatLab ode15s - shocking tolerances to pass though.
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],3.242663439545868e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],0.999, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],2.153277726022381e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],1.145736207245425e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],1.234257806221668e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],1.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6],3.707764665147012e-03, 1e-5);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7],0.5, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8],0.00, 1e-3);
    }
    
    void xTestWntCellCycleSolverWithAPCDoubleHit() throw(Exception)
    {
        double WntLevel = 0.0;
        IngeWntOdeSystem wnt_system(WntLevel,APC_TWO_HIT);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.
        
        double h_value=0.001;
        
        RungeKutta4IvpOdeSolver rk4_solver;
        BackwardEulerIvpOdeSolver back_solver(9);
        
        OdeSolution solutions;
        //OdeSolution solutions2;
        
        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();
        
        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);

        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 3.94 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end] , 3.9435 , 1e-2);
        // Proper values from MatLab ode15s - shocking tolerances to pass though.
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],2.058373151310055e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],0.999, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],3.699024514542648e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],2.687523235896298e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],1.834144555072084e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],0.0, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6],0.5, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7],0.5, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][8],0.00, 1e-3);
    }
    
};

#endif /*TESTINGEWNTODESYSTEM_HPP_*/
