#ifndef _TESTSIMPLENONLINEARSOLVER_HPP_
#define _TESTSIMPLENONLINEARSOLVER_HPP_

#include "SimpleNonLinearSolver.hpp"
#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "petscmat.h"
  
class TestSimpleNonlinearSolver : public CxxTest::TestSuite 
{
//public:
//    void setUp()
//    {
//       	PetscInitialize(&sFakeArgc, &sFakeArgv, PETSC_NULL, 0);
//    }   
//        
//    /* What We need to do:
//     * 
//     * 1. initialize pesky vectors
//     * 2. instantiate nonlinearEllipticPde
//     * 3. instantiate nonlinearIntegrator
//     * 
//     * 
//     * 
//     * 4. output solution
//     * 5. check against exact solution
//     *  
//    */	
//
//	void donttestLinearSolverEasy( void )
//    {
//    // Solve Ax=b. 2x2 matrix
//
//	SimpleLinearSolver solver;
//
//	// Set rhs vector
//	Vec rhs_vector;
//	VecCreate(PETSC_COMM_WORLD, &rhs_vector);
//	VecSetSizes(rhs_vector,PETSC_DECIDE,2);
//	VecSetType(rhs_vector, VECSEQ);
//   	VecSetValue(rhs_vector, 0, (PetscReal) 1, INSERT_VALUES);
//   	VecSetValue(rhs_vector, 1, (PetscReal) 1, INSERT_VALUES);
//   	
//   	//Set Matrix
//   	Mat lhs_matrix;
//   	MatCreate(PETSC_COMM_WORLD, 2, 2, PETSC_DETERMINE, PETSC_DETERMINE, &lhs_matrix);
//   	MatSetType(lhs_matrix, MATSEQDENSE);
//   	
//   	// Set Matrix to Identity matrix
//   	MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 0, (PetscReal) 1, INSERT_VALUES);
//	MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 1, (PetscReal) 0, INSERT_VALUES);
//	MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 0, (PetscReal) 0, INSERT_VALUES);
//	MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 1, (PetscReal) 1, INSERT_VALUES);
//	
//	// Assemble matrix
//   	MatAssemblyBegin(lhs_matrix, MAT_FINAL_ASSEMBLY);
//   	MatAssemblyEnd(lhs_matrix, MAT_FINAL_ASSEMBLY);
//   	
//    // Call solver
//    Vec lhs_vector;
//    TS_ASSERT_THROWS_NOTHING(lhs_vector = solver.Solve(lhs_matrix, rhs_vector));
//    
//    // Check result
//    PetscScalar *lhs_elements;
//    VecGetArray(lhs_vector, &lhs_elements);
//    TS_ASSERT_DELTA(lhs_elements[0], 1.0, 0.000001);
//    TS_ASSERT_DELTA(lhs_elements[1], 1.0, 0.000001);
//    
//    }
//    
//	void testLinearSolverThrowsIfDoesNotConverge( void )
//    {
//    // Solve Ax=b. 2x2 matrix
//	SimpleLinearSolver solver;
//
//	// Set rhs vector
//	Vec rhs_vector;
//	VecCreate(PETSC_COMM_WORLD, &rhs_vector);
//	VecSetSizes(rhs_vector,PETSC_DECIDE,2);
//	VecSetType(rhs_vector, VECSEQ);
//   	VecSetValue(rhs_vector, 0, (PetscReal) 1, INSERT_VALUES);
//   	VecSetValue(rhs_vector, 1, (PetscReal) 1, INSERT_VALUES);
//   	
//   	//Set Matrix
//   	Mat lhs_matrix;
//   	MatCreate(PETSC_COMM_WORLD, 2, 2, PETSC_DETERMINE, PETSC_DETERMINE, &lhs_matrix);
//   	MatSetType(lhs_matrix, MATSEQDENSE);
//   	
//   	// Set Matrix to Zero matrix
//   	MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 0, (PetscReal) 0, INSERT_VALUES);
//	MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 1, (PetscReal) 0, INSERT_VALUES);
//	MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 0, (PetscReal) 0, INSERT_VALUES);
//	MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 1, (PetscReal) 0, INSERT_VALUES);
//	
//	// Assemble matrix
//   	MatAssemblyBegin(lhs_matrix, MAT_FINAL_ASSEMBLY);
//   	MatAssemblyEnd(lhs_matrix, MAT_FINAL_ASSEMBLY);
//   	
//    // Call solver
//    Vec lhs_vector;
//    
//    TS_ASSERT_THROWS_ANYTHING(lhs_vector = solver.Solve(lhs_matrix, rhs_vector));
//        
//    }
//    
//	void donttestLinearSolverHarder( void )
//    {
//    // Solve Ax=b. 2x2 matrix
//	SimpleLinearSolver solver;
//
//	// Set rhs vector
//	Vec rhs_vector;
//	VecCreate(PETSC_COMM_WORLD, &rhs_vector);
//	VecSetSizes(rhs_vector,PETSC_DECIDE,2);
//	VecSetType(rhs_vector, VECSEQ);
//   	VecSetValue(rhs_vector, 0, (PetscReal) 17, INSERT_VALUES);
//   	VecSetValue(rhs_vector, 1, (PetscReal) 39, INSERT_VALUES);
//   	
//   	//Set Matrix
//   	Mat lhs_matrix;
//   	MatCreate(PETSC_COMM_WORLD, 2, 2, PETSC_DETERMINE, PETSC_DETERMINE, &lhs_matrix);
//   	MatSetType(lhs_matrix, MATSEQDENSE);
//   	
//   	// Set Matrix to Zero matrix
//   	MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 0, (PetscReal) 1, INSERT_VALUES);
//	MatSetValue(lhs_matrix, (PetscInt) 0, (PetscInt) 1, (PetscReal) 2, INSERT_VALUES);
//	MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 0, (PetscReal) 3, INSERT_VALUES);
//	MatSetValue(lhs_matrix, (PetscInt) 1, (PetscInt) 1, (PetscReal) 4, INSERT_VALUES);
//	
//	// Assemble matrix
//   	MatAssemblyBegin(lhs_matrix, MAT_FINAL_ASSEMBLY);
//   	MatAssemblyEnd(lhs_matrix, MAT_FINAL_ASSEMBLY);
//   	
//    // Call solver
//    Vec lhs_vector;
//    TS_ASSERT_THROWS_NOTHING(lhs_vector = solver.Solve(lhs_matrix, rhs_vector));
//    
//    // Check result
//    PetscScalar *lhs_elements;
//    VecGetArray(lhs_vector, &lhs_elements);
//    TS_ASSERT_DELTA(lhs_elements[0], 5.0, 0.000001);
//    TS_ASSERT_DELTA(lhs_elements[1], 6.0, 0.000001);
//    
//    }
    
};

#endif //_TESTSIMPLENONLINEARSOLVER_HPP_
