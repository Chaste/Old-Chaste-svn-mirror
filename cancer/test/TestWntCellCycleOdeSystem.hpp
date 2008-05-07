/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef TESTWNTCELLCYCLEODESYSTEM_HPP_
#define TESTWNTCELLCYCLEODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include <stdio.h>
#include <time.h>
#include <vector>
#include <iostream>

#include "WntCellCycleOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CellMutationStates.hpp"


class TestWntCellCycleOdeSystem : public CxxTest::TestSuite
{
public:

    void TestWntCellCycleEquations()
    {
        double wnt_level = 0.0;
        WntCellCycleOdeSystem wnt_cell_cycle_system(wnt_level);
        
        double time = 0.0;
        std::vector<double> initial_conditions = wnt_cell_cycle_system.GetInitialConditions();
        
        std::vector<double> derivs(initial_conditions.size());
        wnt_cell_cycle_system.EvaluateYDerivatives(time, initial_conditions, derivs);
        TS_ASSERT_DELTA(initial_conditions[6]+initial_conditions[7],0.0074,1e-4);
        
        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from MatLab code)
        TS_ASSERT_DELTA(derivs[0],-1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2],-9.370533804903016e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[3],-7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4],1.549680000000000e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[5],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8],0.0, 1e-5);
        
        /**
         * And the same for a high Wnt level
         */
        wnt_level = 1.0;
        WntCellCycleOdeSystem wnt_cell_cycle_system2(wnt_level,LABELLED);
        initial_conditions = wnt_cell_cycle_system2.GetInitialConditions();
        //std::cout << "mutation 0 beta-cat = " << initial_conditions[6] << "\n";
        
        TS_ASSERT_DELTA(initial_conditions[6]+initial_conditions[7],0.6002,1e-4);
        
        wnt_cell_cycle_system2.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from MatLab code)
        TS_ASSERT_DELTA(derivs[0],-1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2],5.845305754146019e+00, 1e-5);
        TS_ASSERT_DELTA(derivs[3],-7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4],1.549680000000000e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[5],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8],0.0, 1e-5);
        
        /**
         * A test for the case mutation = 1
         * (An APC +/- mutation)
         */
        CellMutationState mutation = APC_ONE_HIT;
        wnt_level = 1.0;
        WntCellCycleOdeSystem wnt_cell_cycle_system3(wnt_level,mutation);
        initial_conditions = wnt_cell_cycle_system3.GetInitialConditions();
        //std::cout << "mutation " << mutation << " beta-cat = " << initial_conditions[6] << "\n";
        TS_ASSERT_DELTA(initial_conditions[6]+initial_conditions[7],0.750207,1e-6);
        
        wnt_cell_cycle_system3.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from MatLab code)
        TS_ASSERT_DELTA(derivs[0],-1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2],7.3260e+00, 1e-4);
        TS_ASSERT_DELTA(derivs[3],-7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4],1.549680000000000e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[5],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8],0.0, 1e-5);
        
        /**
        * A test for the case mutation = 2
        * (A beta-cat delta45 mutation)
        */
        mutation = BETA_CATENIN_ONE_HIT;
        wnt_level = 1.0;
        WntCellCycleOdeSystem wnt_cell_cycle_system4(wnt_level,mutation);
        initial_conditions = wnt_cell_cycle_system4.GetInitialConditions();
        //std::cout << "mutation " << mutation << " beta-cat = " << initial_conditions[6] << "\n";
        TS_ASSERT_DELTA(initial_conditions[6]+initial_conditions[7],0.8001,1e-4);
        
        wnt_cell_cycle_system4.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from MatLab code)
        TS_ASSERT_DELTA(derivs[0],-1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2],7.8190e+00, 1e-4);
        TS_ASSERT_DELTA(derivs[3],-7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4],1.549680000000000e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[5],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8],0.0, 1e-5);
        
        /**
        * A test for the case mutation = 3
        * (An APC -/- mutation)
        */
        mutation = APC_TWO_HIT;
        wnt_level = 1.0;
        WntCellCycleOdeSystem wnt_cell_cycle_system5(wnt_level,mutation);
        initial_conditions = wnt_cell_cycle_system5.GetInitialConditions();
        //std::cout << "mutation " << mutation << " beta-cat = " << initial_conditions[6] << "\n";
        
        TS_ASSERT_DELTA(initial_conditions[6]+initial_conditions[7],1.0,1e-6);
        
        wnt_cell_cycle_system5.EvaluateYDerivatives(time, initial_conditions, derivs);
        
        // Test derivatives are correct at t=0 for these initial conditions
        // (figures from MatLab code)
        TS_ASSERT_DELTA(derivs[0],-1.586627673253325e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-5.532201118824132e-05, 1e-5);
        TS_ASSERT_DELTA(derivs[2],9.7928e+00, 1e-4);
        TS_ASSERT_DELTA(derivs[3],-7.449833887043188e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[4],1.549680000000000e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[5],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[7],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[8],0.0, 1e-5);
    }
    
    void TestWntCellCycleSolver() throw(Exception)
    {
        double wnt_level = 1.0;
        WntCellCycleOdeSystem wnt_system(wnt_level,LABELLED);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.
        
        double h_value_rk4=1e-4;
        
        RungeKutta4IvpOdeSolver rk4_solver;
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        BackwardEulerIvpOdeSolver back_solver(9);
        
        OdeSolution solutions_rk4;
        OdeSolution solutions_rkf;
        
        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();
                
        double start_time, end_time, elapsed_time = 0.0;
        start_time = std::clock();
        solutions_rk4 = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value_rk4, h_value_rk4);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "1. Runge-Kutta Elapsed time = " << elapsed_time << "\n";
        
        double h_value_rkf = 0.1;
        
        initial_conditions = wnt_system.GetInitialConditions();
        start_time = std::clock();
        solutions_rkf = rkf_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value_rkf, 1e-4);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "2. Runge-Kutta-Fehlberg Elapsed time = " << elapsed_time << "\n";
        
//        WntCellCycleOdeSystem wnt_system_2(wnt_level);
//        initial_conditions = wnt_system.GetInitialConditions();
//
//        h_value = 0.001;
//
//        start_time = std::clock();
//        solutions = back_solver.Solve(&wnt_system_2, initial_conditions, 0.0, 100.0, h_value, h_value);
//        end_time = std::clock();
//        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
//        std::cout <<  "1. BackwardEuler Elapsed time = " << elapsed_time << "\n";



        
        // Testing RK4 solution
        // Test solutions are OK for a small time increase...
        int end = solutions_rk4.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 5.971 hours)
        TS_ASSERT_DELTA(solutions_rk4.rGetTimes()[end] , 5.971 , 1e-2);
        // Proper values from MatLab ode15s - shocking tolerances to pass though.
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][0],2.880603485931000e-01, 1e-3);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][1],1.000220438771564e+00, 1.02e-2);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][2],2.453870380958196e+00, 1e-3);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][3],1.446185835615586e+00, 1e-3);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][4],1.383272155041549e-01, 1e-3);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][5],4.975124378109454e-03, 1e-3);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][6]+solutions_rk4.rGetSolutions()[end][7],6.002649406788524e-01, 1e-3);
        TS_ASSERT_DELTA(solutions_rk4.rGetSolutions()[end][8],1.00, 1e-3);
        
        
        // Testing RKF solution
        // Test solutions are OK for a small time increase...
        end = solutions_rkf.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 5.971 hours)
        TS_ASSERT_DELTA(solutions_rkf.rGetTimes()[end] , 5.971 , 1e-2);
        // Proper values from MatLab ode15s - shocking tolerances to pass though.
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][0],2.880603485931000e-01, 1e-3);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][1],1.000220438771564e+00, 1.02e-2);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][2],2.453870380958196e+00, 1e-3);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][3],1.446185835615586e+00, 1e-3);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][4],1.383272155041549e-01, 1e-3);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][5],4.975124378109454e-03, 1e-3);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][6]+solutions_rkf.rGetSolutions()[end][7],6.002649406788524e-01, 1e-3);
        TS_ASSERT_DELTA(solutions_rkf.rGetSolutions()[end][8],1.00, 1e-3);               
    }
    
    void TestWntCellCycleSolverWithAPCSingleHit() throw(Exception)
    {
        double wnt_level = 1.0;
        WntCellCycleOdeSystem wnt_system(wnt_level,APC_ONE_HIT);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.
        
        double h_value=0.0001;
        
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
    
    void TestWntCellCycleSolverWithBetaCateninHit() throw(Exception)
    {
        double wnt_level = 0.0;
        WntCellCycleOdeSystem wnt_system(wnt_level,BETA_CATENIN_ONE_HIT);
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
        
//
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
    
    void TestWntCellCycleSolverWithAPCDoubleHit() throw(Exception)
    {
        double wnt_level = 0.0;
        WntCellCycleOdeSystem wnt_system(wnt_level,APC_TWO_HIT);
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

#endif /*TESTWNTCELLCYCLEODESYSTEM_HPP_*/
