#ifndef TESTGARYWNTODESYSTEM_HPP_
#define TESTGARYWNTODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>
#include "RungeKutta4IvpOdeSolver.hpp"
#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
#include "GaryWntOdeSystem.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"

class TestGaryWntOdeSystem : public CxxTest::TestSuite
{
private:
    void CheckDerivativesZero(CryptCellMutationState mutation, double wnt_level, double total_beta_cat)
    {
        double time = 0.0;
        GaryWntOdeSystem wnt_ode_system(wnt_level,mutation);
        std::vector<double> initial_conditions = wnt_ode_system.GetInitialConditions();
        std::vector<double> derivs(initial_conditions.size());

        TS_ASSERT_DELTA(initial_conditions[1]+initial_conditions[2],total_beta_cat,1e-4);
        
        wnt_ode_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        TS_ASSERT_DELTA(derivs[0],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[1],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[2],0.0, 1e-5);
        TS_ASSERT_DELTA(derivs[3],0.0, 1e-5);        
    }
    
public:
    void TestDerivatesZero()
    {
        CheckDerivativesZero(HEALTHY, 0.0, 0.0074);
        CheckDerivativesZero(HEALTHY, 1.0, 0.6002);
        CheckDerivativesZero(LABELLED,1.0, 0.6002);
        CheckDerivativesZero(APC_ONE_HIT, 1.0, 0.750207);
        CheckDerivativesZero(BETA_CATENIN_ONE_HIT, 1.0, 0.8001);
        CheckDerivativesZero(APC_TWO_HIT, 1.0, 1.0);
    }
    
    void TestWntCellCycleSolver() throw(Exception)
    {
        double wnt_level = 1.0;
        GaryWntOdeSystem wnt_system(wnt_level,LABELLED);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.
        
        double h_value=0.001;
        
        RungeKutta4IvpOdeSolver rk4_solver;
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        
        OdeSolution solutions;
        //OdeSolution solutions2;
        
        std::vector<double> initial_conditions = wnt_system.GetInitialConditions();
                
        double start_time, end_time, elapsed_time = 0.0;
        start_time = std::clock();
        solutions = rk4_solver.Solve(&wnt_system, initial_conditions, 0.0, 5.971, h_value, h_value);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "1. Runge-Kutta Elapsed time = " << elapsed_time << "\n";
        
        h_value = 0.1;
        
        initial_conditions = wnt_system.GetInitialConditions();
        start_time = std::clock();
        solutions = rkf_solver.Solve(&wnt_system, initial_conditions, 0.0, 5.971, h_value, 1e-4);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "2. Runge-Kutta-Fehlberg Elapsed time = " << elapsed_time << "\n";
        


//        int my_rank;
//        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
//        if (my_rank==0) // if master process
//        {
//        
//            int step_per_row = 1;
//            ColumnDataWriter writer("WntCellCycle","WntCellCycle");
//            int time_var_id = writer.DefineUnlimitedDimension("Time","s");
//            
//            std::vector<int> var_ids;
//            for (unsigned i=0; i<wnt_system.rGetVariableNames().size(); i++)
//            {
//                var_ids.push_back(writer.DefineVariable(wnt_system.rGetVariableNames()[i],
//                                                        wnt_system.rGetVariableUnits()[i]));
//            }
//            writer.EndDefineMode();
//            
//            for (unsigned i = 0; i < solutions.rGetSolutions().size(); i+=step_per_row)
//            {
//                writer.PutVariable(time_var_id, solutions.rGetTimes()[i]);
//                for (unsigned j=0; j<var_ids.size(); j++)
//                {
//                    writer.PutVariable(var_ids[j], solutions.rGetSolutions()[i][j]);
//                }
//                writer.AdvanceAlongUnlimitedDimension();
//            }
//            writer.Close();
//        }
//        MPI_Barrier(PETSC_COMM_WORLD);
        
        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 5.971 hours)
//        TS_ASSERT_DELTA(solutions.rGetTimes()[end] , 5.971 , 1e-2);
        // Proper values from MatLab ode15s - shocking tolerances to pass though.
//        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],2.880603485931000e-01, 1e-3);
//        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],1.000220438771564e+00, 1.01e-2);
//        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],2.453870380958196e+00, 1e-3);
//        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],1.446185835615586e+00, 1e-3);
//        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],1.383272155041549e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],4.975124378109454e-03, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1]+solutions.rGetSolutions()[end][2],6.002649406788524e-01, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],1.00, 1e-3);
    }    
};

#endif /*TESTGARYWNTODESYSTEM_HPP_*/
