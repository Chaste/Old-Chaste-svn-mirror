#ifndef TESTBACKWARDEULERIVPODESOLVER_HPP_
#define TESTBACKWARDEULERIVPODESOLVER_HPP_

#include <cxxtest/TestSuite.h>

#include "OdeThirdOrder.hpp"
//#include "AnotherOde.hpp"
#include "BackwardEulerIvpOdeSolver.cpp"

#include "PetscSetupAndFinalize.hpp"

class TestBackwardEulerIvpOdeSolver: public CxxTest::TestSuite
{
public:
  /*      void NO_testComputeResidualAndComputeJacobian()
    {
        AnotherOde ode_system;

        // create BackwardEulerStructure
        BackwardEulerStructure backward_euler_structure;
        backward_euler_structure.TimeStep = 0.1;
        backward_euler_structure.Time = 0;
        backward_euler_structure.pAbstractOdeSystem = &ode_system;

        std::vector<double> current_y_value;
        current_y_value.push_back(1.0);
        current_y_value.push_back(2.0);
        current_y_value.push_back(3.0);
        backward_euler_structure.currentYValue = current_y_value;

        // create solution guess
        Vec initial_guess;
        VecCreate(PETSC_COMM_WORLD, &initial_guess);
        VecSetSizes(initial_guess, PETSC_DECIDE,3);
        VecSetFromOptions(initial_guess);
        for (int i=0;i<3; i++)
        {
            VecSetValue(initial_guess, i, current_y_value[i] ,INSERT_VALUES);

        }
        VecAssemblyBegin(initial_guess);
        VecAssemblyEnd(initial_guess);

        // compute residual
        Vec residual;
        VecDuplicate(initial_guess, &residual);
        SNES snes;
        PetscErrorCode ierr;
        ierr = ComputeResidual(snes, initial_guess, residual, &backward_euler_structure);


        // check results
        PetscScalar *p_residual_array;
        ierr = VecGetArray(residual, &p_residual_array);

        TS_ASSERT_DELTA(p_residual_array[0],-2.0, 1e-6);
        TS_ASSERT_DELTA(p_residual_array[1], 1.0, 1e-6);
        TS_ASSERT_DELTA(p_residual_array[2],-1.0, 1e-6);

        ierr = VecRestoreArray(residual, &p_residual_array);

        // NOW TEST COMPUTE JACOBIAN

        Mat jacobian;
        Mat preconditioner;
        MatStructure mat_structure;
        ierr = MatCreate(PETSC_COMM_WORLD, &jacobian);
        ierr = MatSetSizes(jacobian, PETSC_DECIDE, PETSC_DECIDE, 3, 3);
        ierr = MatSetFromOptions(jacobian);

        ierr = ComputeJacobian(snes, initial_guess, &jacobian, &preconditioner, &mat_structure, &backward_euler_structure);

        double true_jacobian[3][3] = {{ 3,-1,-1},
                                      {-1, 1,-1},
                                      { 0, 2,-1}};
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

                TS_ASSERT_DELTA(ret_array[0],true_jacobian[row][col], 1e-6);
            }
        }
    }
  */

    void testBackwardEulerSystemOf3Equations() throw (Exception)
    {
        OdeThirdOrder ode_system;
        
        double h_value=0.01;
        
        //Euler solver solution worked out
        BackwardEulerIvpOdeSolver backward_euler_solver;
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
    
    void testComputeResidual() throw (Exception)
    {
    	double h_value=1.0;
    	OdeThirdOrder ode_system;
    	
        BackwardEulerStructure backward_euler_structure;
        backward_euler_structure.TimeStep = h_value;
        backward_euler_structure.Time = 0.0;
        backward_euler_structure.pAbstractOdeSystem = &ode_system;  
        std::vector<double> current_y_value;
        current_y_value.push_back(1.0);
        current_y_value.push_back(2.0);
        current_y_value.push_back(3.0);
        backward_euler_structure.currentYValue = current_y_value;        
        
        Vec solution_guess, residual;
        int indices[3] = {0,1,2};
        double values[3] = {1.0, 2.0, 3.0};
        SNES snes;
        
        VecCreate(PETSC_COMM_WORLD,solution_guess);
        VecSetSizes(solution_guess,PETSC_DECIDE,3);
        VecSetFromOptions(solution_guess);
        VecDuplicate(solution_guess,&residual);
        
        VecSetValues(solution_guess,3,indices,values,INSERT_VALUES);
        VecAssemblyBegin(solution_guess);
        VecAssemblyEnd(solution_guess);
        
        PetscErrorCode compute_residual = ComputeResidual(snes,solution_guess,residual,&backward_euler_structure);
        
        PetscScalar *p_residual_array;
        ierr = VecGetArray(residual, &p_residual_array);

        TS_ASSERT_DELTA(p_residual_array[0],-2.0, 1e-6);
        TS_ASSERT_DELTA(p_residual_array[1], 1.0, 1e-6);
        TS_ASSERT_DELTA(p_residual_array[2],-1.0, 1e-6);

        VecRestoreArray(residual, &p_residual_array);
        
  
    }   
    
    
};

#endif /*TESTBACKWARDEULERIVPODESOLVER_HPP_*/
