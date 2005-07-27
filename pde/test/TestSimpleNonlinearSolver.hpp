#ifndef _TESTSIMPLENONLINEARSOLVER_HPP_
#define _TESTSIMPLENONLINEARSOLVER_HPP_

#include "SimpleNonlinearSolver.hpp"
#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "petscmat.h"
#include <iostream>
#include <cmath>

#include "PetscSetupAndFinalize.hpp"

PetscErrorCode ComputeTestResidual(SNES snes,Vec solutionGuess,Vec residual,void *pContext);  
PetscErrorCode ComputeTestJacobian(SNES snes,Vec input,Mat *pJacobian ,Mat *pPreconditioner,MatStructure *pMatStructure ,void *pContext);
PetscErrorCode ComputeTestResidual3d(SNES snes,Vec solutionGuess,Vec residual,void *pContext);  
PetscErrorCode ComputeTestJacobian3d(SNES snes,Vec input,Mat *pJacobian ,Mat *pPreconditioner,MatStructure *pMatStructure ,void *pContext);
  
class TestSimpleNonlinearSolver : public CxxTest::TestSuite 
{
public:
    
    void TestOn2dNonlinearProblem(void)
    {
    	SimpleNonlinearSolver solver;
    	
    	// Set up solution guess for residuals
    	int length=2;
		    	
    	// Set up initial Guess
    	Vec initialGuess;
    	VecCreate(PETSC_COMM_WORLD, &initialGuess);
    	VecSetSizes(initialGuess, PETSC_DECIDE,length);
    	//VecSetType(initialGuess, VECSEQ);
    	VecSetFromOptions(initialGuess);
    	VecSetValue(initialGuess, 0, 1.0 ,INSERT_VALUES);
		VecSetValue(initialGuess, 1, 1.0 ,INSERT_VALUES);
    	VecAssemblyBegin(initialGuess);
		VecAssemblyEnd(initialGuess);
		
    	Vec answer;
    	Vec residual;
    	VecDuplicate(initialGuess, &residual);
    	    	
 		answer = solver.Solve(&ComputeTestResidual, &ComputeTestJacobian,
 							  residual, initialGuess, NULL);
    	
    	PetscScalar *answerElements;
		VecGetArray(answer, &answerElements);
		double x = answerElements[0];
		double y = answerElements[1];
		VecRestoreArray(answer,&answerElements);
    	
    	double tol = 1e-6;
    	TS_ASSERT_DELTA(x,1/sqrt(2),tol);
    	TS_ASSERT_DELTA(y,1/sqrt(2),tol);
    	
    	VecDestroy(initialGuess);
    	VecDestroy(residual);
    	VecDestroy(answer);
    }
    
    void TestOn3dNonlinearProblem(void)
    {
    	SimpleNonlinearSolver solver;
    	
    	// Set up solution guess for residuals
    	int length=3;
		    	
    	// Set up initial Guess
    	Vec initialGuess;
    	VecCreate(PETSC_COMM_WORLD, &initialGuess);
    	VecSetSizes(initialGuess, PETSC_DECIDE,length);
    	//VecSetType(initialGuess, VECSEQ);
    	VecSetFromOptions(initialGuess);
    	VecSetValue(initialGuess, 0, 1.0 ,INSERT_VALUES);
		VecSetValue(initialGuess, 1, 1.0 ,INSERT_VALUES);
		VecSetValue(initialGuess, 2, 1.0 ,INSERT_VALUES);
    	VecAssemblyBegin(initialGuess);
		VecAssemblyEnd(initialGuess);
		
    	Vec answer;
    	Vec residual;
    	VecDuplicate(initialGuess,&residual);
    	    	
 		answer = solver.Solve(&ComputeTestResidual3d, &ComputeTestJacobian3d,
 							  residual, initialGuess, NULL);
    	
    	PetscScalar *answerElements;
		VecGetArray(answer, &answerElements);
		double x = answerElements[0];
		double y = answerElements[1];
		double z = answerElements[2];
		VecRestoreArray(answer,&answerElements);
    	
    	double tol = 1e-6;
    	TS_ASSERT_DELTA(x,1/sqrt(3),tol);
    	TS_ASSERT_DELTA(y,1/sqrt(3),tol);
    	TS_ASSERT_DELTA(z,1/sqrt(3),tol);
    	
    	VecDestroy(initialGuess);
    	VecDestroy(residual);
    	VecDestroy(answer);
	}

};

PetscErrorCode ComputeTestResidual(SNES snes,Vec solutionGuess,Vec residual,void *pContext)
{
	double x,y;
	
	PetscScalar *solutionGuessElements;
	VecGetArray(solutionGuess, &solutionGuessElements);
	x = solutionGuessElements[0];
	y = solutionGuessElements[1];
	
	VecRestoreArray(solutionGuess,&solutionGuessElements);
	VecSetValue(residual,0,x*x+y*y-1,INSERT_VALUES);
	VecSetValue(residual,1,x-y,INSERT_VALUES);
	VecAssemblyBegin(residual);
	VecAssemblyEnd(residual);
	return 0;
}

PetscErrorCode ComputeTestJacobian(SNES snes,Vec input,Mat *pJacobian ,Mat *pPreconditioner,MatStructure *pMatStructure ,void *pContext)
{
	double x, y;
		
	PetscScalar *inputElements;
	VecGetArray(input, &inputElements);
	x = inputElements[0];
	y = inputElements[1];
	VecRestoreArray(input,&inputElements);
		
	MatSetValue(*pJacobian, 0 , 0 , 2.0*x , INSERT_VALUES);
	MatSetValue(*pJacobian, 0 , 1 , 2.0*y, INSERT_VALUES);
	MatSetValue(*pJacobian, 1 , 0 , 1.0, INSERT_VALUES);
	MatSetValue(*pJacobian, 1 , 1 , -1.0, INSERT_VALUES);
    MatAssemblyBegin(*pJacobian,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*pJacobian,MAT_FINAL_ASSEMBLY);
	
	return 0;
}

PetscErrorCode ComputeTestResidual3d(SNES snes,Vec solutionGuess,Vec residual,void *pContext)
{
	double x,y,z;
	
	PetscScalar *solutionGuessElements;
	VecGetArray(solutionGuess, &solutionGuessElements);
	x = solutionGuessElements[0];
	y = solutionGuessElements[1];
	z = solutionGuessElements[2];
	
	VecRestoreArray(solutionGuess,&solutionGuessElements);
	VecSetValue(residual,0,x*x+y*y+z*z-1,INSERT_VALUES);
	VecSetValue(residual,1,x-y,INSERT_VALUES);
	VecSetValue(residual,2,y-z,INSERT_VALUES);
	VecAssemblyBegin(residual);
	VecAssemblyEnd(residual);
	return 0;
}

PetscErrorCode ComputeTestJacobian3d(SNES snes,Vec input,Mat *pJacobian ,Mat *pPreconditioner,MatStructure *pMatStructure ,void *pContext)
{
	double x, y, z;
		
	PetscScalar *inputElements;
	VecGetArray(input, &inputElements);
	x = inputElements[0];
	y = inputElements[1];
	z = inputElements[2];
	VecRestoreArray(input,&inputElements);
		
	MatSetValue(*pJacobian, 0 , 0 , 2.0*x , INSERT_VALUES);
	MatSetValue(*pJacobian, 0 , 1 , 2.0*y, INSERT_VALUES);
	MatSetValue(*pJacobian, 0 , 2 , 2.0*z, INSERT_VALUES);
	MatSetValue(*pJacobian, 1 , 0 , 1.0, INSERT_VALUES);
	MatSetValue(*pJacobian, 1 , 1 , -1.0, INSERT_VALUES);
	MatSetValue(*pJacobian, 1 , 2 , 0.0, INSERT_VALUES);
	MatSetValue(*pJacobian, 2 , 0 , 0.0, INSERT_VALUES);
	MatSetValue(*pJacobian, 2 , 1 , 1.0, INSERT_VALUES);
	MatSetValue(*pJacobian, 2 , 2 , -1.0, INSERT_VALUES);
    MatAssemblyBegin(*pJacobian,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*pJacobian,MAT_FINAL_ASSEMBLY);
	
	return 0;
}

#endif //_TESTSIMPLENONLINEARSOLVER_HPP_
