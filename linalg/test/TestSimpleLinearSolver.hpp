#ifndef _TESTSIMPLELINEARSOLVER_HPP_
#define _TESTSIMPLELINEARSOLVER_HPP_

#include "SimpleLinearSolver.hpp"
#include <cxxtest/TestSuite.h>
#include <petsc.h>
 
#include "PetscSetupAndFinalize.hpp"

class TestSimpleLinearSolver : public CxxTest::TestSuite 
{
public:
	void TestLinearSolverEasy( void )
	{
		// Solve Ax=b. 2x2 matrix
	
		SimpleLinearSolver solver;
	
		// Set rhs vector
		Vec rhs_vector;
		VecCreate(PETSC_COMM_WORLD, &rhs_vector);
		VecSetSizes(rhs_vector,PETSC_DECIDE,2);
		//VecSetType(rhs_vector, VECSEQ);
	   	VecSetFromOptions(rhs_vector);
	   	
	   	VecSetValue(rhs_vector, 0, (PetscReal) 1, INSERT_VALUES);
	   	VecSetValue(rhs_vector, 1, (PetscReal) 1, INSERT_VALUES);
	   	
	   	//Set Matrix
	   	Mat lhs_matrix;
#if (PETSC_VERSION_MINOR == 2) //Old API
	   	MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2, 2, &lhs_matrix);
#else
        MatCreate(PETSC_COMM_WORLD,&lhs_matrix);
        MatSetSizes(lhs_matrix, PETSC_DECIDE, PETSC_DECIDE,2,2);
#endif
	   	//MatSetType(lhs_matrix, MATSEQDENSE);
	   	MatSetType(lhs_matrix, MATMPIDENSE);
	   	
	   	// Set Matrix to Identity matrix
	   	MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 0, (PetscReal) 1, INSERT_VALUES);
		MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 1, (PetscReal) 0, INSERT_VALUES);
		MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 0, (PetscReal) 0, INSERT_VALUES);
		MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 1, (PetscReal) 1, INSERT_VALUES);
		
		// Assemble matrix
	   	MatAssemblyBegin(lhs_matrix, MAT_FINAL_ASSEMBLY);
	   	MatAssemblyEnd(lhs_matrix, MAT_FINAL_ASSEMBLY);
	   	
		// Call solver
		Vec lhs_vector;
		TS_ASSERT_THROWS_NOTHING(lhs_vector = solver.Solve(lhs_matrix, rhs_vector, 2));
		
		// Check result
		PetscScalar *p_lhs_elements_array;
		VecGetArray(lhs_vector, &p_lhs_elements_array);
        int lo, hi;
        VecGetOwnershipRange(lhs_vector, &lo, &hi);
        
        for(int global_index=0; global_index<2; global_index++)
        {
            int local_index = global_index-lo;
            if(lo<=global_index && global_index<hi)
            {    
                TS_ASSERT_DELTA(p_lhs_elements_array[local_index], 1.0, 0.000001);
            }
        }

		// Free memory
		VecDestroy(rhs_vector);
		VecDestroy(lhs_vector);
		MatDestroy(lhs_matrix);
	}
	
	void TestLinearSolverThrowsIfDoesNotConverge( void )
	{
		// Solve Ax=b. 2x2 matrix
		SimpleLinearSolver solver;
	
		// Set rhs vector
		Vec rhs_vector;
		VecCreate(PETSC_COMM_WORLD, &rhs_vector);
		VecSetSizes(rhs_vector,PETSC_DECIDE,2);
		//VecSetType(rhs_vector, VECSEQ);
	   	VecSetFromOptions(rhs_vector);
	   	VecSetValue(rhs_vector, 0, (PetscReal) 1, INSERT_VALUES);
	   	VecSetValue(rhs_vector, 1, (PetscReal) 1, INSERT_VALUES);
	   	
	   	//Set Matrix
	   	Mat lhs_matrix;
#if (PETSC_VERSION_MINOR == 2) //Old API
        MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2, 2, &lhs_matrix);
#else
        MatCreate(PETSC_COMM_WORLD,&lhs_matrix);
        MatSetSizes(lhs_matrix, PETSC_DECIDE, PETSC_DECIDE,2,2);
#endif	   	
    //MatSetType(lhs_matrix, MATSEQDENSE);
	   	MatSetType(lhs_matrix, MATMPIDENSE);
	   	
	   	// Set Matrix to Zero matrix
	   	MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 0, (PetscReal) 0, INSERT_VALUES);
		MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 1, (PetscReal) 0, INSERT_VALUES);
		MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 0, (PetscReal) 0, INSERT_VALUES);
		MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 1, (PetscReal) 0, INSERT_VALUES);
		
		// Assemble matrix
	   	MatAssemblyBegin(lhs_matrix, MAT_FINAL_ASSEMBLY);
	   	MatAssemblyEnd(lhs_matrix, MAT_FINAL_ASSEMBLY);
	   	
		// Call solver
		Vec lhs_vector;
		
		TS_ASSERT_THROWS_ANYTHING(lhs_vector = solver.Solve(lhs_matrix, rhs_vector, 2));
		
		// Free memory
		VecDestroy(rhs_vector);
		MatDestroy(lhs_matrix);
	}
	
	void TestLinearSolverHarder( void )
	{
		// Solve Ax=b. 2x2 matrix
		SimpleLinearSolver solver;
	
		// Set rhs vector
		Vec rhs_vector;
		VecCreate(PETSC_COMM_WORLD, &rhs_vector);
		VecSetSizes(rhs_vector,PETSC_DECIDE,2);
		//VecSetType(rhs_vector, VECSEQ);
		VecSetFromOptions(rhs_vector);
		
	   	VecSetValue(rhs_vector, 0, (PetscReal) 17, INSERT_VALUES);
	   	VecSetValue(rhs_vector, 1, (PetscReal) 39, INSERT_VALUES);
	   	
	   	//Set Matrix
	   	Mat lhs_matrix;
#if (PETSC_VERSION_MINOR == 2) //Old API
        MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2, 2, &lhs_matrix);
#else
        MatCreate(PETSC_COMM_WORLD,&lhs_matrix);
        MatSetSizes(lhs_matrix, PETSC_DECIDE, PETSC_DECIDE,2,2);
#endif
	   	//MatSetType(lhs_matrix, MATSEQDENSE);
	   	MatSetType(lhs_matrix, MATMPIDENSE);
	   	
	   	// Set Matrix values
	   	MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 0, (PetscReal) 1, INSERT_VALUES);
		MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 1, (PetscReal) 2, INSERT_VALUES);
		MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 0, (PetscReal) 3, INSERT_VALUES);
		MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 1, (PetscReal) 4, INSERT_VALUES);
		
		// Assemble matrix
	   	MatAssemblyBegin(lhs_matrix, MAT_FINAL_ASSEMBLY);
	   	MatAssemblyEnd(lhs_matrix, MAT_FINAL_ASSEMBLY);
	   	
		// Call solver
		Vec lhs_vector;
		TS_ASSERT_THROWS_NOTHING(lhs_vector = solver.Solve(lhs_matrix, rhs_vector, 2));
		
		// Check result
		PetscScalar *p_lhs_elements_array;
		VecGetArray(lhs_vector, &p_lhs_elements_array);
        int lo, hi;
        VecGetOwnershipRange(lhs_vector, &lo, &hi);


        // p_lhs_elements_array should be equal to [5,6]
        for(int global_index=0; global_index<2; global_index++)
        {
            int local_index = global_index-lo;
            if(lo<=global_index && global_index<hi)
            {    
                TS_ASSERT_DELTA(p_lhs_elements_array[local_index], global_index+5.0, 0.000001);
            }
        }

		// Free memory
		VecDestroy(rhs_vector);
		VecDestroy(lhs_vector);
		MatDestroy(lhs_matrix);
	}
    
    // This test illustrate what happens when solving a singular 
    // linear system in a couple of cases
    void TestLinearSolverWithSingularMatrix( void )
    {
        // Solve Ax=b. 2x2 matrix
        //
        // A = 6 0, b = 3
        //     0 0      3
        SimpleLinearSolver solver;
    
        // Set rhs vector
        Vec rhs_vector;
        VecCreate(PETSC_COMM_WORLD, &rhs_vector);
        VecSetSizes(rhs_vector,PETSC_DECIDE,2);
        VecSetFromOptions(rhs_vector);
        
        VecSetValue(rhs_vector, 0, (PetscReal) 3, INSERT_VALUES);
        VecSetValue(rhs_vector, 1, (PetscReal) 0, INSERT_VALUES);
        
        //Set Matrix
        Mat lhs_matrix;
#if (PETSC_VERSION_MINOR == 2) //Old API
        MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2, 2, &lhs_matrix);
#else
        MatCreate(PETSC_COMM_WORLD,&lhs_matrix);
        MatSetSizes(lhs_matrix, PETSC_DECIDE, PETSC_DECIDE,2,2);
#endif
        MatSetType(lhs_matrix, MATMPIDENSE);
        
        // Set Matrix values
        MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 0, (PetscReal) 6, INSERT_VALUES);
        MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 1, (PetscReal) 0, INSERT_VALUES);
        MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 0, (PetscReal) 0, INSERT_VALUES);
        MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 1, (PetscReal) 0, INSERT_VALUES);
        
        // Assemble matrix
        MatAssemblyBegin(lhs_matrix, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(lhs_matrix, MAT_FINAL_ASSEMBLY);
        
        // Call solver
        Vec lhs_vector;
        TS_ASSERT_THROWS_NOTHING(lhs_vector = solver.Solve(lhs_matrix, rhs_vector, 2));
        
        // Check result
        PetscScalar *p_lhs_elements_array;
        VecGetArray(lhs_vector, &p_lhs_elements_array);
        int lo, hi;
        VecGetOwnershipRange(lhs_vector, &lo, &hi);

        // p_lhs_elements_array[0] should be equal to 0.5
        int global_index = 0;
        int local_index = global_index-lo;
        
        if(lo<=global_index && global_index<hi)
        { 
            TS_ASSERT_DELTA( p_lhs_elements_array[local_index], 0.5, 1e-6);
        }    
            
        // Now change the rhs vector to from (3,0) to (3,3). There are no solutions
        // so check a petsc error occurs
        VecSetValue(rhs_vector, 1, (PetscReal) 3, INSERT_VALUES);
        TS_ASSERT_THROWS_ANYTHING(lhs_vector = solver.Solve(lhs_matrix, rhs_vector, 2));

        // Free memory
        VecDestroy(rhs_vector);
        VecDestroy(lhs_vector);
        MatDestroy(lhs_matrix);
    }
    
    void testLinearSolverWithMatrixIsConstantAndNullSpace()
    {
        // Solve Ax=b. 5x5 matrix
        SimpleLinearSolver solver;
    
        // Set rhs vector
        Vec rhs_vector;
        VecCreate(PETSC_COMM_WORLD, &rhs_vector);
        VecSetSizes(rhs_vector,PETSC_DECIDE,5);
        //VecSetType(rhs_vector, VECSEQ);
        VecSetFromOptions(rhs_vector);
        
        //Set Matrix
        Mat lhs_matrix;
#if (PETSC_VERSION_MINOR == 2) //Old API
        MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 5, 5, &lhs_matrix);
#else
        MatCreate(PETSC_COMM_WORLD,&lhs_matrix);
        MatSetSizes(lhs_matrix, PETSC_DECIDE, PETSC_DECIDE,5,5);
#endif
        MatSetType(lhs_matrix, MATMPIDENSE);

        Vec null_basis_vector;
        VecDuplicate(rhs_vector, &null_basis_vector);

        for(int i=0; i<5; i++)
        {
            VecSetValue(rhs_vector, i, (PetscReal) i, INSERT_VALUES);
            VecSetValue(null_basis_vector, i, 1.0, INSERT_VALUES); 
            for(int j=0; j<5; j++)
            {
                double val=0;
                if(i==j)
                {
                    val = -2;
                }
                else if( (i==j+1) || (i==j-1) )
                {
                    val = 1;
                }
                MatSetValue(lhs_matrix, i, j, val, INSERT_VALUES);
            }
        }
        
        // Assemble matrix
        MatAssemblyBegin(lhs_matrix, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(lhs_matrix, MAT_FINAL_ASSEMBLY);
        
        // set matrix is constant
        solver.SetMatrixIsConstant();
        
        MatNullSpace mat_null_space;
        MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, 1, &null_basis_vector, &mat_null_space);
        
        // Call solver
        Vec lhs_vector;
        TS_ASSERT_THROWS_NOTHING(lhs_vector = solver.Solve(lhs_matrix, rhs_vector, 5, mat_null_space));
        
        // Check result
        PetscScalar *p_lhs_elements_array;
        VecGetArray(lhs_vector, &p_lhs_elements_array);
        int lo, hi;
        VecGetOwnershipRange(lhs_vector, &lo, &hi);
        
        double answers[] = {-3.33, -6.66, -9.0, -9.33, -6.66};
        for (int global_index = lo; global_index<hi; global_index++)
        {
            int local_index = global_index-lo;        
            TS_ASSERT_DELTA(p_lhs_elements_array[local_index], answers[local_index], 0.1);
        }    
    }
        
};

#endif //_TESTSIMPLELINEARSOLVER_HPP_
