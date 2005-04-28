#ifndef _TESTSIMPLELINEARSOLVER_HPP_
#define _TESTSIMPLELINEARSOLVER_HPP_

#include "SimpleLinearSolver.hpp"
#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "petscmat.h"
  
class TestSimpleLinearSolver : public CxxTest::TestSuite 
{
public:
    void setUp()
    {
		int FakeArgc=0;
		char *FakeArgv0="testrunner";
		char **FakeArgv=&FakeArgv0;
   		PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
    }   
        

	void testLinearSolverEasy( void )
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
   	MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2, 2, &lhs_matrix);
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
    TS_ASSERT_THROWS_NOTHING(lhs_vector = solver.Solve(lhs_matrix, rhs_vector));
    
    // Check result
    PetscScalar *lhs_elements;
    VecGetArray(lhs_vector, &lhs_elements);
    TS_ASSERT_DELTA(lhs_elements[0], 1.0, 0.000001);
    TS_ASSERT_DELTA(lhs_elements[1], 1.0, 0.000001);
    
    }
    
	void testLinearSolverThrowsIfDoesNotConverge( void )
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
   	MatCreate(PETSC_COMM_WORLD,  PETSC_DECIDE, PETSC_DECIDE, 2, 2, &lhs_matrix);
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
    
    TS_ASSERT_THROWS_ANYTHING(lhs_vector = solver.Solve(lhs_matrix, rhs_vector));
        
    }
    
	void testLinearSolverHarder( void )
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
   	MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2, 2, &lhs_matrix);
   	//MatSetType(lhs_matrix, MATSEQDENSE);
   	MatSetType(lhs_matrix, MATMPIDENSE);
   	
   	// Set Matrix to Zero matrix
   	MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 0, (PetscReal) 1, INSERT_VALUES);
	MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 1, (PetscReal) 2, INSERT_VALUES);
	MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 0, (PetscReal) 3, INSERT_VALUES);
	MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 1, (PetscReal) 4, INSERT_VALUES);
	
	// Assemble matrix
   	MatAssemblyBegin(lhs_matrix, MAT_FINAL_ASSEMBLY);
   	MatAssemblyEnd(lhs_matrix, MAT_FINAL_ASSEMBLY);
   	
    // Call solver
    Vec lhs_vector;
    TS_ASSERT_THROWS_NOTHING(lhs_vector = solver.Solve(lhs_matrix, rhs_vector));
    
    // Check result
    PetscScalar *lhs_elements;
    VecGetArray(lhs_vector, &lhs_elements);
    TS_ASSERT_DELTA(lhs_elements[0], 5.0, 0.000001);
    TS_ASSERT_DELTA(lhs_elements[1], 6.0, 0.000001);
    //TS_TRACE("here simp lin\n");
    }
    
};

#endif //_TESTSIMPLELINEARSOLVER_HPP_
