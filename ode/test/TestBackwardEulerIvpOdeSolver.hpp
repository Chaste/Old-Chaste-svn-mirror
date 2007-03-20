#ifndef TESTBACKWARDEULERIVPODESOLVER_HPP_
#define TESTBACKWARDEULERIVPODESOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include<iostream>

#include "OdeThirdOrder.hpp"
#include "OdeThirdOrderWithEvents.hpp"
#include "Ode4.hpp"
#include "Ode5.hpp"
#include "Ode5Jacobian.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestBackwardEulerIvpOdeSolver: public CxxTest::TestSuite
{
public:
    void TestBackwardEulerSystemOf3Equations()
    {
        OdeThirdOrder ode_system;
        
        double h_value=0.01;
        
        //Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(ode_system.GetNumberOfStateVariables());
        
        // cover the SetEpsilonForNumericalJacobian() method
        backward_euler_solver.SetEpsilonForNumericalJacobian(1e-6);
        
        OdeSolution solutions;
        
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = backward_euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last = solutions.GetNumberOfTimeSteps();
        
        double numerical_solution[3];
        numerical_solution[0] = solutions.rGetSolutions()[last][0];
        numerical_solution[1] = solutions.rGetSolutions()[last][1];
        numerical_solution[2] = solutions.rGetSolutions()[last][2];
        
        // The tests
        double analytical_solution[3];
        
        analytical_solution[0] = -sin(2);
        analytical_solution[1] = sin(2)+cos(2);
        analytical_solution[2] = 2*sin(2);
        
        double global_error_euler = 0.5*2*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(numerical_solution[0],analytical_solution[0],global_error_euler);
        TS_ASSERT_DELTA(numerical_solution[1],analytical_solution[1],global_error_euler);
        TS_ASSERT_DELTA(numerical_solution[2],analytical_solution[2],global_error_euler);
    }
    
    void TestBackwardEulerNonlinearEquation()
    {
        Ode4 ode_system;
        
        double h_value=0.01;
        
        //Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(ode_system.GetNumberOfStateVariables());
        OdeSolution solutions;
        
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        
        solutions = backward_euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last = solutions.GetNumberOfTimeSteps();
        
        double numerical_solution;
        numerical_solution = solutions.rGetSolutions()[last][0];
        
        // The tests
        double analytical_solution = 1.0/(1.0+exp(-12.5));
        
        TS_ASSERT_DELTA(numerical_solution,analytical_solution,1.0e-4);
    }
    
    void TestBackwardEulerAnotherNonlinearEquation()
    {
        Ode5 ode_system;
        
        double h_value=0.01;
        double end_time = 1.0;
        
        //Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(ode_system.GetNumberOfStateVariables());
        OdeSolution solutions;
        
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        
        solutions = backward_euler_solver.Solve(&ode_system, state_variables, 0.0, end_time, h_value, h_value);
        int last = solutions.GetNumberOfTimeSteps();
        
        double numerical_solution;
        numerical_solution = solutions.rGetSolutions()[last][0];
        
        // The tests
        double analytical_solution = 1.0/(1.0+4.0*exp(-100.0*end_time));
        
        TS_ASSERT_DELTA(numerical_solution,analytical_solution,1.0e-3);
    }
    
    void TestBackwardEulerSystemOf3EquationsWithEvents()
    {
        OdeThirdOrderWithEvents ode_system_with_events;
        
        double h_value=0.01;
        
        //Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(ode_system_with_events.GetNumberOfStateVariables());
        OdeSolution solutions;
        
        std::vector<double> state_variables = ode_system_with_events.GetInitialConditions();
        solutions = backward_euler_solver.Solve(&ode_system_with_events, state_variables, 0.0, 2.0, h_value, h_value);
        int last = solutions.GetNumberOfTimeSteps();
        
        // final time should be pi/6 (?)
        TS_ASSERT_DELTA( solutions.rGetTimes()[last], 0.5236, 0.01);
        
        // penultimate y0 should be greater than -0.5
        TS_ASSERT_LESS_THAN(-0.5,solutions.rGetSolutions()[last-1][0]);
        // final y0 should be less than -0.5
        TS_ASSERT_LESS_THAN( solutions.rGetSolutions()[last][0], -0.5);
        
        // solver should correctly state the stopping event occured
        TS_ASSERT_EQUALS(backward_euler_solver.StoppingEventOccured(), true);
    }
    
    void TestBackwardEulerAnotherNonlinearEquationAnalytic() throw(Exception)
    {
        Ode5Jacobian ode_system;
        
        double h_value = 0.01;
        double end_time = 1.0;
        
        //Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver(ode_system.GetNumberOfStateVariables());
        OdeSolution solutions;
        
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        
        solutions = backward_euler_solver.Solve(&ode_system, state_variables, 0.0, end_time, h_value, h_value);
        int last = solutions.GetNumberOfTimeSteps();
        
        double numerical_solution;
        numerical_solution = solutions.rGetSolutions()[last][0];
        
        // The tests
        double analytical_solution = 1.0/(1.0+4.0*exp(-100.0*end_time));
        
        TS_ASSERT_DELTA(numerical_solution,analytical_solution,1.0e-3);
    }
    
    
};

#endif /*TESTBACKWARDEULERIVPODESOLVER_HPP_*/
