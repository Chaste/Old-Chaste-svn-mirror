#ifndef _TESTABSTRACTANALYTICJACOBIAN_HPP_
#define _TESTABSTRACTANALYTICJACOBIAN_HPP_

// TestAbstractAnalyticJacobian.hpp

#include <cmath>
#include <iostream>
#include <vector>
#include "Jacobian1.hpp"
#include "Jacobian2.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"

// Tolerance for tests
double tol=0.01;


class TestAbstractAnalyticJacobian : public CxxTest::TestSuite
{
public:

    void TestJacobianOne(void)
    {
        // pointer to TestOde1 class
        Jacobian1 jacobian1;
        
        // Make a PetSc vector of Y values...        
        Vec solution_guess;
        int indices[1] = {0};
        double values[1] = {2.0};
        VecCreate(PETSC_COMM_WORLD,&solution_guess);
        VecSetSizes(solution_guess,PETSC_DECIDE,1);
        VecSetFromOptions(solution_guess);
        VecSetValues(solution_guess,1,indices,values,INSERT_VALUES);
        VecAssemblyBegin(solution_guess);
        VecAssemblyEnd(solution_guess);
        
        // Set up a Jacobian matrix for function to put values in
        Mat jacobian;
        //MatStructure mat_structure;
        #if (PETSC_VERSION_MINOR == 2) //Old API
                MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 1, 1, &jacobian);
        #else
                MatCreate(PETSC_COMM_WORLD, &jacobian);
                MatSetSizes(jacobian, PETSC_DECIDE, PETSC_DECIDE, 1, 1);
        #endif
        MatSetFromOptions(jacobian);
        
        // This is the function we are testing...
        jacobian1.AnalyticJacobian(solution_guess, &jacobian, 1.0, 0.01);
        
        
        // Get the answer out! 
        int row = 0;
        int col = 0;
        int row_as_array[1];
        row_as_array[0] = row;
        int col_as_array[1];
        col_as_array[0] = col;
        double ret_array[1];
        MatGetValues(jacobian, 1, row_as_array, 1, col_as_array, ret_array);
                
        TS_ASSERT_DELTA(ret_array[0], 0.96, tol);      
    }
    
    void TestJacobianTwo(void)
    {
        // pointer to TestOde1 class
        Jacobian2 jacobian2;
//        const int numEqns = jacobian2.GetNumberOfStateVariables();
        
        // Make a PetSc vector of Y values...        
        Vec solution_guess;
        int indices[2] = {0,1};
        double values[2] = {1.0,2.0};
        VecCreate(PETSC_COMM_WORLD,&solution_guess);
        VecSetSizes(solution_guess,PETSC_DECIDE,2);
        VecSetFromOptions(solution_guess);
        VecSetValues(solution_guess,2,indices,values,INSERT_VALUES);
        VecAssemblyBegin(solution_guess);
        VecAssemblyEnd(solution_guess);
        
        // Set up a Jacobian matrix for function to put values in
        Mat jacobian;
        //MatStructure mat_structure;
        #if (PETSC_VERSION_MINOR == 2) //Old API
                MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2, 2, &jacobian);
        #else
                MatCreate(PETSC_COMM_WORLD, &jacobian);
                MatSetSizes(jacobian, PETSC_DECIDE, PETSC_DECIDE, 2, 2);
        #endif
        MatSetFromOptions(jacobian);
        
        // This is the function we are testing...
        jacobian2.AnalyticJacobian(solution_guess, &jacobian, 1.0, 0.01);
                
        // Get the answer out!
        double big_vector[2][2]; 
        for(int row = 0; row<2; row++)
        {
            for(int col = 0; col<2 ; col++)
            {
                int row_as_array[1];
                row_as_array[0] = row;
                int col_as_array[1];
                col_as_array[0] = col;
                double ret_array[1];
                MatGetValues(jacobian, 1, row_as_array, 1, col_as_array, ret_array);
                big_vector[row][col]=ret_array[0];
            }
        }        
        TS_ASSERT_DELTA(big_vector[0][0], 0.98, tol);
        TS_ASSERT_DELTA(big_vector[0][1], -0.04, tol);
        TS_ASSERT_DELTA(big_vector[1][0], -0.02, tol);
        TS_ASSERT_DELTA(big_vector[1][1], 0.92, tol);
              
    }
        
};



#endif //_TESTABSTRACTANALYTICJACOBIAN_HPP_
