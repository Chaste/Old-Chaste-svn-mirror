#ifndef _TESTSIMPLENONLINEARSOLVER_HPP_
#define _TESTSIMPLENONLINEARSOLVER_HPP_

#include "SimpleNonlinearSolver.hpp"
#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "petscmat.h"
#include <iostream>
  
PetscErrorCode ComputeTestResidual(SNES snes,Vec solutionGuess,Vec residual,void *pContext);  
PetscErrorCode ComputeTestJacobian(SNES snes,Vec input,Mat *pJacobian ,Mat *pPreconditioner,MatStructure *pMatStructure ,void *pContext);
  
class TestSimpleNonlinearSolver : public CxxTest::TestSuite 
{
public:
    void setUp()
    {
		PetscInitialize(&sFakeArgc, &sFakeArgv, PETSC_NULL, 0);
    }   
        
    
    void testOn2dNonlinearProblem(void)
    {
    	SimpleNonlinearSolver *solver=new SimpleNonlinearSolver();
    	
    	// Set up solution guess for residuals
    	int length=2;
		    	
    	// Set up initial Guess
    	Vec initialGuess;
    	VecCreate(PETSC_COMM_WORLD, &initialGuess);
    	VecSetSizes(initialGuess, PETSC_DECIDE,length);
    	VecSetType(initialGuess, VECSEQ);
    	VecSetValue(initialGuess, 0, 1.0 ,INSERT_VALUES);
		VecSetValue(initialGuess, 1, 1.0 ,INSERT_VALUES);
    	VecAssemblyBegin(initialGuess);
		VecAssemblyEnd(initialGuess);
		
    	Vec answer;
    	Vec residual;
    	VecDuplicate(initialGuess,&residual);
    	VecDuplicate(initialGuess,&answer);
    	    	
 		answer=solver->Solve(&ComputeTestResidual, &ComputeTestJacobian, residual, initialGuess, NULL);
    	
    	PetscScalar *answerElements;
		VecGetArray(answer, &answerElements);
		double x = answerElements[0];
		double y = answerElements[1];
		VecRestoreArray(answer,&answerElements);
    	
    	double tol = 1e-5;
    	TS_ASSERT_DELTA(x,0.707107,tol);
    	TS_ASSERT_DELTA(y,0.707107,tol);
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

#endif //_TESTSIMPLENONLINEARSOLVER_HPP_
