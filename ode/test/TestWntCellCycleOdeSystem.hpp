#ifndef TESTWNTCELLCYCLEODESYSTEM_HPP_
#define TESTWNTCELLCYCLEODESYSTEM_HPP_

#include <stdio.h>
#include <time.h>
#include <cxxtest/TestSuite.h>
#include "WntCellCycleOdeSystem.hpp"
#include <vector>
#include <iostream>
#include "RungeKutta4IvpOdeSolver.hpp"
#include "ColumnDataWriter.hpp"

#include "PetscSetupAndFinalize.hpp"


class TestWntCellCycleOdeSystem : public CxxTest::TestSuite
{
public:

    void testWntCellCycleEquations()
    {
		WntCellCycleOdeSystem wnt_cell_cycle_system;
        
        double time = 0.0;
        std::vector<double> initialConditions;
        
        // these at start of G1 phase values
        initialConditions.push_back(7.357000000000000e-01);
        initialConditions.push_back(1.713000000000000e-01);
        initialConditions.push_back(6.900000000000001e-02);
        initialConditions.push_back(3.333333333333334e-03);
        initialConditions.push_back(1.000000000000000e-04);
        
        // These three at equilibrium values
        double a2d = 2.57e-03;
        double a3d = 3.44e-01;
        double a4d = 1000.0;
        double a5d = 5.0;
        
        double WntLevel = 0.0;
        double destruction_level = a5d/(a4d*WntLevel+a5d);
        double beta_cat_level = a2d/(a2d+a3d*destruction_level);
        
        initialConditions.push_back(destruction_level);
        initialConditions.push_back(beta_cat_level);
        initialConditions.push_back(WntLevel);// Wnt level
        std::vector<double> derivs = wnt_cell_cycle_system.EvaluateYDerivatives(time,initialConditions);
        
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
        destruction_level = a5d/(a4d*WntLevel+a5d);
        beta_cat_level = a2d/(a2d+a3d*destruction_level);
        
        initialConditions[5] = destruction_level;
        initialConditions[6] = beta_cat_level;
        initialConditions[7] = WntLevel;// Wnt level
		derivs = wnt_cell_cycle_system.EvaluateYDerivatives(time,initialConditions);
        
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
        WntCellCycleOdeSystem wnt_system;
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.

        double h_value=0.001;
        
        RungeKutta4IvpOdeSolver rk4_solver;

        OdeSolution solutions;
                
        std::vector<double> initialConditions;
        // these at start of G1 phase values
        initialConditions.push_back(7.357000000000000e-01);
        initialConditions.push_back(1.713000000000000e-01);
        initialConditions.push_back(6.900000000000001e-02);
        initialConditions.push_back(3.333333333333334e-03);
        initialConditions.push_back(1.000000000000000e-04);
        // These three at equilibrium values
        double a2d = 2.57e-03;
        double a3d = 3.44e-01;
        double a4d = 1000.0;
        double a5d = 5.0;
        double WntLevel = 1.0;
        double destruction_level = a5d/(a4d*WntLevel+a5d);
        double beta_cat_level = a2d/(a2d+a3d*destruction_level);
        
        initialConditions.push_back(destruction_level);
        initialConditions.push_back(beta_cat_level);
        initialConditions.push_back(WntLevel);// Wnt level
                        
        solutions = rk4_solver.Solve(&wnt_system, initialConditions, 0.0, 100.0, h_value, h_value);
        
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

        // Test backward euler solutions are OK for a very small time increase...
        int end = solutions.rGetSolutions().size() - 1;
    	// Tests the simulation is ending at the right time...(going into S phase at 5.971 hours)
    	TS_ASSERT_DELTA(solutions.rGetTimes()[end] , 5.971 , 1e-3);
        // Proper values from MatLab ode15s - shocking tolerances to pass though.
		TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],2.880603485931000e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],1.000220438771564e+00, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],2.453870380958196e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],1.446185835615586e+00, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],1.383272155041549e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],4.975124378109454e-03, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6],6.002649406788524e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][7],1.00, 1e-4);
    }
    
};

#endif /*TESTWNTCELLCYCLEODESYSTEM_HPP_*/
