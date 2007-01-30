#ifndef TESTBACKWARDEULERIVPODESOLVER_HPP_
#define TESTBACKWARDEULERIVPODESOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include<iostream>

#include "OdeThirdOrder.hpp"
#include "OdeThirdOrderWithEvents.hpp"
#include "Ode4.hpp"
#include "Ode5.hpp"
#include "Ode5Jacobian.hpp"
#include "BackwardEulerIvpOdeSolver.cpp"
#include "BetterBackwardEulerIvpOdeSolver.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestBackwardEulerIvpOdeSolver: public CxxTest::TestSuite
{
public:
    void testBackwardEulerSystemOf3Equations()
    {
        OdeThirdOrder ode_system;
        
        double h_value=0.01;
        
        //Euler solver solution worked out
        BetterBackwardEulerIvpOdeSolver<3> backward_euler_solver;
        
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
    
    void testBackwardEulerNonlinearEquation()
    {
        Ode4 ode_system;
        
        double h_value=0.01;
        
        //Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver;
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
    
    void testBackwardEulerAnotherNonlinearEquation()
    {
        Ode5 ode_system;
        
        double h_value=0.01;
        double end_time = 1.0;
        
        //Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver;
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
    
    // Superceded test
    void notestComputeResidual()
    {
        double h_value=1.0;
        OdeThirdOrder ode_system;
        BackwardEulerIvpOdeSolver backward_euler_solver;
        BackwardEulerStructure backward_euler_structure;
        backward_euler_structure.TimeStep = h_value;
        backward_euler_structure.Time = 0.0;
        backward_euler_structure.pAbstractOdeSystem = &ode_system;
        std::vector<double> current_y_value;
        current_y_value.push_back(1.0);
        current_y_value.push_back(2.0);
        current_y_value.push_back(3.0);
        backward_euler_structure.currentYValues = current_y_value;
        
        Vec solution_guess, residual;
        int indices[3] = {0,1,2};
        double values[3] = {1.0, 2.0, 3.0};
        SNES snes;
        
        VecCreate(PETSC_COMM_WORLD,&solution_guess);
        VecSetSizes(solution_guess,PETSC_DECIDE,3);
        VecSetFromOptions(solution_guess);
        VecDuplicate(solution_guess,&residual);
        
        VecSetValues(solution_guess,3,indices,values,INSERT_VALUES);
        VecAssemblyBegin(solution_guess);
        VecAssemblyEnd(solution_guess);
        
        ComputeResidual(snes,solution_guess,residual,&backward_euler_structure);
        
        PetscScalar *p_residual_array;
        VecGetArray(residual, &p_residual_array);
        
        TS_ASSERT_DELTA(p_residual_array[0],-2.0, 1e-6);
        TS_ASSERT_DELTA(p_residual_array[1], 1.0, 1e-6);
        TS_ASSERT_DELTA(p_residual_array[2],-1.0, 1e-6);
        
        VecRestoreArray(residual, &p_residual_array);
        
    }
    
    // Superceded test
    void notestComputeJacobian()
    {
        double h_value=1.0;
        OdeThirdOrder ode_system;
        BackwardEulerIvpOdeSolver backward_euler_solver;
        BackwardEulerStructure backward_euler_structure;
        backward_euler_structure.TimeStep = h_value;
        backward_euler_structure.Time = 0.0;
        backward_euler_structure.pAbstractOdeSystem = &ode_system;
        std::vector<double> current_y_value;
        current_y_value.push_back(1.0);
        current_y_value.push_back(2.0);
        current_y_value.push_back(3.0);
        backward_euler_structure.currentYValues = current_y_value;
        
        Vec solution_guess;
        int indices[3] = {0,1,2};
        double values[3] = {1.0, 2.0, 3.0};
        
        SNES snes;
        
        VecCreate(PETSC_COMM_WORLD,&solution_guess);
        VecSetSizes(solution_guess,PETSC_DECIDE,3);
        VecSetFromOptions(solution_guess);
        
        VecSetValues(solution_guess,3,indices,values,INSERT_VALUES);
        VecAssemblyBegin(solution_guess);
        VecAssemblyEnd(solution_guess);
        
        Mat jacobian;
        Mat preconditioner;
        MatStructure mat_structure;
#if (PETSC_VERSION_MINOR == 2) //Old API
        MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 3, 3, &jacobian);
#else
        MatCreate(PETSC_COMM_WORLD, &jacobian);
        MatSetSizes(jacobian, PETSC_DECIDE, PETSC_DECIDE, 3, 3);
#endif
        MatSetFromOptions(jacobian);
        
        ComputeNumericalJacobian(snes, solution_guess, &jacobian, &preconditioner, &mat_structure, &backward_euler_structure);
        
        double true_jacobian[3][3] = {{ 0, 1,-1},
                                      { 0, 0, 1},
                                      { 0,-2, 2}};
        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                int row_as_array[1];
                row_as_array[0] = row;
                int col_as_array[1];
                col_as_array[0] = col;
                
                double ret_array[1];
                
                MatGetValues(jacobian, 1, row_as_array, 1, col_as_array, ret_array);
                
                TS_ASSERT_DELTA(ret_array[0],true_jacobian[row][col], 1e-3);
            }
        }
    }
    
    void testBackwardEulerSystemOf3EquationsWithEvents()
    {
        OdeThirdOrderWithEvents ode_system_with_events;
        
        double h_value=0.01;
        
        //Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver;
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
    
    void testBackwardEulerAnotherNonlinearEquationAnalytic() throw(Exception)
    {
        Ode5Jacobian ode_system;
        
        double h_value = 0.01;
        double end_time = 1.0;
        
        //Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver;
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
    
    void testBetterBackwardEulerIvpOdeSolver() throw(Exception)
    {
        Ode5Jacobian ode_system;
        
        double h_value = 0.01;
        double end_time = 1.0;
        
        //Euler solver solution worked out
        BetterBackwardEulerIvpOdeSolver<1> backward_euler_solver;
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
