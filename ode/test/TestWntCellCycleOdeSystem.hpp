#ifndef TESTWNTCELLCYCLEODESYSTEM_HPP_
#define TESTWNTCELLCYCLEODESYSTEM_HPP_

#include <stdio.h>
#include <time.h>
#include <cxxtest/TestSuite.h>
#include "WntCellCycleOdeSystem.hpp"
#include <vector>
#include <iostream>
#include "RungeKutta4IvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "ColumnDataWriter.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestWntCellCycleOdeSystem : public CxxTest::TestSuite
{
public:

    void testWntCellCycleEquations()
    {
    	double WntLevel = 0.0;
		WntCellCycleOdeSystem wnt_cell_cycle_system(WntLevel);
        
        double time = 0.0;
        std::vector<double> initial_conditions = wnt_cell_cycle_system.GetInitialConditions();
        
        std::vector<double> derivs(initial_conditions.size());
        wnt_cell_cycle_system.EvaluateYDerivatives(time, initial_conditions, derivs);
        
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
        
        /** 
         * And the same for a high Wnt level
         */
        WntLevel = 1.0;
        WntCellCycleOdeSystem wnt_cell_cycle_system2(WntLevel);
		initial_conditions = wnt_cell_cycle_system2.GetInitialConditions();
		
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
    }
    
    void testWntCellCycleSolver() throw(Exception)
    {
    	double WntLevel = 1.0;
        WntCellCycleOdeSystem wnt_system(WntLevel);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.

        double h_value=0.001;
        
        RungeKutta4IvpOdeSolver rk4_solver;
        BackwardEulerIvpOdeSolver back_solver(8);

        OdeSolution solutions;
        //OdeSolution solutions2;
                
        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();
              
              
        double start_time, end_time, elapsed_time = 0.0;
        start_time = std::clock();                
        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "1. Runge-Kutta Elapsed time = " << elapsed_time << "\n";
        
//        WntCellCycleOdeSystem wnt_system_2(WntLevel);
//        initial_conditions = wnt_system.GetInitialConditions();
//        
//        h_value = 0.001;
//        
//        start_time = std::clock();                
//        solutions = back_solver.Solve(&wnt_system_2, initial_conditions, 0.0, 100.0, h_value, h_value);
//        end_time = std::clock();
//        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
//        std::cout <<  "1. BackwardEuler Elapsed time = " << elapsed_time << "\n";
        
        
        int my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        if (my_rank==0) // if master process
        {

            int step_per_row = 1;
            ColumnDataWriter writer("WntCellCycle","WntCellCycle");
            int time_var_id = writer.DefineUnlimitedDimension("Time","s");

            std::vector<int> var_ids;
            for (unsigned i=0; i<wnt_system.rGetVariableNames().size(); i++)
            {
                var_ids.push_back(writer.DefineVariable(wnt_system.rGetVariableNames()[i],
                                                    wnt_system.rGetVariableUnits()[i]));
            }
            writer.EndDefineMode();

            for (unsigned i = 0; i < solutions.rGetSolutions().size(); i+=step_per_row)
            {
                writer.PutVariable(time_var_id, solutions.rGetTimes()[i]);
                for (unsigned j=0; j<var_ids.size(); j++)
                {
                    writer.PutVariable(var_ids[j], solutions.rGetSolutions()[i][j]);
                }
                writer.AdvanceAlongUnlimitedDimension();
            }
            writer.Close();
        }
        MPI_Barrier(PETSC_COMM_WORLD);

        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
    	// Tests the simulation is ending at the right time...(going into S phase at 5.971 hours)
    	TS_ASSERT_DELTA(solutions.rGetTimes()[end] , 5.971 , 1e-2);
        // Proper values from MatLab ode15s - shocking tolerances to pass though.
		TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],2.880603485931000e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],1.000220438771564e+00, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],2.453870380958196e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],1.446185835615586e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],1.383272155041549e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],4.975124378109454e-03, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6],6.002649406788524e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7],1.00, 1e-3);
    }
    
};

#endif /*TESTWNTCELLCYCLEODESYSTEM_HPP_*/
