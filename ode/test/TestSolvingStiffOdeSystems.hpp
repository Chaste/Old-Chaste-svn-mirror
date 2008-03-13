#ifndef TESTSOLVINGSTIFFODESYSTEMS_HPP_
#define TESTSOLVINGSTIFFODESYSTEMS_HPP_

#include <cxxtest/TestSuite.h>

#include <vector>
#include <iostream>

#include "BackwardEulerIvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "FieldNoyesReactionSystem.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestSolvingStiffOdeSystems : public CxxTest::TestSuite
{
public:

    /**
     *  Solve the Field-Noyes system. 
     *  This is a stiff ode so Runge-Kutta won't work - program hangs if
     *  end time > 0.01 and variables go out of bounds.
     *  Can be solved ok with backward Euler.
     */
    void TestFieldNoyesReactionSystem()
    {
        FieldNoyesReactionSystem ode_system;
        RungeKutta4IvpOdeSolver rk4_solver;
        OdeSolution ode_solution;
        std::vector<double> ode_state_variables = ode_system.GetInitialConditions();
        
        // note small end time
        ode_solution = rk4_solver.Solve(&ode_system, ode_state_variables, 0.0, 0.01, 0.001, 0.001);
        
        int last = ode_solution.GetNumberOfTimeSteps();
        for (int i=0; i<=last; i++)
        {
            double x = ode_solution.rGetSolutions()[i][0];
            TS_ASSERT_LESS_THAN(x, 9.7e3);
            TS_ASSERT_LESS_THAN(0.99, x);
        }
        
        BackwardEulerIvpOdeSolver backward_euler_solver(ode_system.GetNumberOfStateVariables());
        OdeSolution ode_solution2;
        std::vector<double> ode_state_variables2 = ode_system.GetInitialConditions();
        ode_solution2 = backward_euler_solver.Solve(&ode_system, ode_state_variables2, 0.0, 0.1, 0.001, 0.001);
        
        last = ode_solution2.GetNumberOfTimeSteps();
        for (int i=0; i<=last; i++)
        {
            double x = ode_solution2.rGetSolutions()[i][0];
            TS_ASSERT_LESS_THAN(x, 9.7e3);
            TS_ASSERT_LESS_THAN(0.99, x);
        }
    }
};

#endif /*TESTSOLVINGSTIFFODESYSTEMS_HPP_*/
