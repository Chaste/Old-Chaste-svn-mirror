#ifndef _TESTPRACTICALONE_HPP_
#define _TESTPRACTICALONE_HPP_

#include "SimpleLinearSolver.hpp"
#include "petscmat.h"
#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "SimpleNonlinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include <vector>
#include <iostream>
#include "Node.hpp" 
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"

#include "Practical1Question1Pde.hpp"
  
class TestPracticalOne : public CxxTest::TestSuite 
{
public:
	
	void testQuestion1(void)
	{
		int FakeArgc=0;
		char *FakeArgv0="testrunner";
		char **FakeArgv=&FakeArgv0;
		PetscInitialize(&FakeArgc, &FakeArgv, PETSC_NULL, 0);
		
		// Create mesh from mesh reader
		
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/practical1_1d_mesh");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		// Instantiate PDE object
		Practical1Question1Pde<1> pde;  
		
		// Boundary conditions
		// u'(0)=0, u(1)=1
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition1 = new ConstBoundaryCondition<1>(0.0);
        
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
        bcc.AddNeumannBoundaryCondition(*iter,pBoundaryCondition1);
        
        ConstBoundaryCondition<1>* pBoundaryCondition2 = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(100), pBoundaryCondition2);
        
		// Linear solver
		SimpleLinearSolver solver;
		
		// Assembler
		SimpleLinearEllipticAssembler<1,1> assembler;
		
		Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
		
		// Check result
		double *res;
		int ierr = VecGetArray(result, &res);
		// Solution should be u = 0.5*x*(3-x)
		for (int i=0; i < mesh.GetNumElements()+1; i++)
		{
			double x = 0.0 + 0.01*i;
			double u = 0.5*(x*x+1);
			TS_ASSERT_DELTA(res[i], u, 0.001);
		}
		VecRestoreArray(result, &res);
	}

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

#endif //_TESTPRACTICALONE_HPP_
