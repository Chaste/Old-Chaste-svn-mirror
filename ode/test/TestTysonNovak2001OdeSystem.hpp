#ifndef TESTTYSONNOVAK2001ODESYSTEM_HPP_
#define TESTTYSONNOVAK2001ODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>
#include "TysonNovak2001OdeSystem.hpp"
#include <vector>
#include <iostream>
#include "BackwardEulerIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "ColumnDataWriter.hpp"

#include "PetscSetupAndFinalize.hpp"


class TestTysonNovak2001OdeSystem : public CxxTest::TestSuite
{
public:

    void testTysonNovak() throw(Exception)
    {
        TysonNovak2001OdeSystem tyson_novak_system;
        
        double time = 0.0;
        std::vector<double> initialConditions; 
        initialConditions.push_back(0.6);
        initialConditions.push_back(0.1);
        initialConditions.push_back(1.5);
        initialConditions.push_back(0.6);
        initialConditions.push_back(0.6);
        initialConditions.push_back(0.85);
        std::vector<double> derivs = tyson_novak_system.EvaluateYDerivatives(time,initialConditions);
        
//		Test derivatives are correct
        TS_ASSERT_DELTA(derivs[0],-4.400000000000000e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[1],-6.047872340425530e+00, 1e-5);
        TS_ASSERT_DELTA(derivs[2],3.361442884485838e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[3],4.016602000735009e-02, 1e-5);
        TS_ASSERT_DELTA(derivs[4],8.400000000000001e-03, 1e-5);
        TS_ASSERT_DELTA(derivs[5],7.777500000000001e-03, 1e-5);
        
////		Solve system using backward Euler solver
//        double h_value=0.0001;
//        
//        //Euler solver solution worked out
//        BackwardEulerIvpOdeSolver backward_euler_solver;
//        RungeKutta4IvpOdeSolver rk4_solver;
//        
//        OdeSolution solutions;
//        OdeSolution solutions2;
//
//        
//        std::vector<double> state_variables = tyson_novak_system.GetInitialConditions();
//        solutions = backward_euler_solver.Solve(&tyson_novak_system, state_variables, 0.0, 0.001, h_value, h_value);
//        solutions2 = rk4_solver.Solve(&tyson_novak_system, state_variables, 0.0, 100.0, h_value, h_value);
//        TS_ASSERT_EQUALS(solutions.GetNumberOfTimeSteps(), 10);
//        
//        int step_per_row = 100;
//        ColumnDataWriter writer("TysonNovak","TysonNovak");
//        int time_var_id = writer.DefineUnlimitedDimension("Time","s");
//        
//        std::vector<int> var_ids;
//        for (unsigned i=0; i<tyson_novak_system.rGetVariableNames().size(); i++)
//        {
//            var_ids.push_back(writer.DefineVariable(tyson_novak_system.rGetVariableNames()[i],
//                                                    tyson_novak_system.rGetVariableUnits()[i]));
//        }
//        writer.EndDefineMode();
//        
//        for (unsigned i = 0; i < solutions2.rGetSolutions().size(); i+=step_per_row)
//        {
//            writer.PutVariable(time_var_id, solutions2.rGetTimes()[i]);
//            for (unsigned j=0; j<var_ids.size(); j++)
//            {
//                writer.PutVariable(var_ids[j], solutions2.rGetSolutions()[i][j]);
//            }
//            writer.AdvanceAlongUnlimitedDimension();
//        }
//        writer.Close();
//        
//
//		// Test backward euler solutions are OK for a very small time increase...
//		int end = solutions.rGetSolutions().size() - 1;
//	
//        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],0.59995781827316, 1e-5);
//        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],0.09406711653612, 1e-5);
//        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],1.50003361032032, 1e-5);
//        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],0.60004016820575, 1e-5);
//        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],0.60000839905560, 1e-5);
//        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],0.85000777753272, 1e-5);
//        
    }
};

#endif /*TESTTYSONNOVAK2001ODESYSTEM_HPP_*/
