#ifndef _TESTSIMPLENONLINEARSOLVER_HPP_
#define _TESTSIMPLENONLINEARSOLVER_HPP_

#include "SimpleNonlinearSolver.hpp"
#include <cxxtest/TestSuite.h>
#include <petsc.h>
#include <petsc.h>
#include <iostream>
#include <cmath>

#include "PetscSetupAndFinalize.hpp"

PetscErrorCode ComputeTestResidual(SNES snes,Vec solution_guess,Vec residual,void *pContext);  
PetscErrorCode ComputeTestJacobian(SNES snes,Vec input,Mat *pJacobian ,Mat *pPreconditioner,MatStructure *pMatStructure ,void *pContext);
PetscErrorCode ComputeTestResidual3d(SNES snes,Vec solution_guess,Vec residual,void *pContext);  
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
    	Vec initial_guess;
    	VecCreate(PETSC_COMM_WORLD, &initial_guess);
    	VecSetSizes(initial_guess, PETSC_DECIDE,length);
    	//VecSetType(initial_guess, VECSEQ);
    	VecSetFromOptions(initial_guess);
    	VecSetValue(initial_guess, 0, 1.0 ,INSERT_VALUES);
		VecSetValue(initial_guess, 1, 1.0 ,INSERT_VALUES);
    	VecAssemblyBegin(initial_guess);
		VecAssemblyEnd(initial_guess);
		
    	Vec answer;
    	Vec residual;
    	VecDuplicate(initial_guess, &residual);
    	 
           	
 		answer = solver.Solve(&ComputeTestResidual, &ComputeTestJacobian,
 							  residual, initial_guess, NULL);
                              
    	
    	PetscScalar *answer_elements;
		VecGetArray(answer, &answer_elements);
        int lo, hi;
        VecGetOwnershipRange(answer, &lo, &hi);
        
	    double tol = 1e-6;

        for(int global_index=0; global_index<2; global_index++)
        {
            int local_index = global_index-lo;
            if(lo<=global_index && global_index<hi)
            {    
                double x = answer_elements[local_index];
                TS_ASSERT_DELTA(x,1/sqrt(2),tol);                
            }
        }
//    	if (lo<=0 && 0<hi)
//        {
//          double x = answer_elements[0-lo];
//	      TS_ASSERT_DELTA(x,1/sqrt(2),tol);
//        }
//        if (lo<=1 && 1<hi)
//        {
//        	double y = answer_elements[1-lo];	
//        	TS_ASSERT_DELTA(y,1/sqrt(2),tol);
//        }
        VecRestoreArray(answer,&answer_elements);
    	
    	VecDestroy(initial_guess);
    	VecDestroy(residual);
    	VecDestroy(answer);
    }
    
    void TestOn3dNonlinearProblem(void)
    {
    	SimpleNonlinearSolver solver;
    	
    	// Set up solution guess for residuals
    	int length=3;
		    	
    	// Set up initial Guess
    	Vec initial_guess;
    	VecCreate(PETSC_COMM_WORLD, &initial_guess);
    	VecSetSizes(initial_guess, PETSC_DECIDE,length);
    	//VecSetType(initial_guess, VECSEQ);
    	VecSetFromOptions(initial_guess);
    	VecSetValue(initial_guess, 0, 1.0 ,INSERT_VALUES);
		VecSetValue(initial_guess, 1, 1.0 ,INSERT_VALUES);
		VecSetValue(initial_guess, 2, 1.0 ,INSERT_VALUES);
    	VecAssemblyBegin(initial_guess);
		VecAssemblyEnd(initial_guess);
		
    	Vec answer;
    	Vec residual;
    	VecDuplicate(initial_guess,&residual);
    	    	
 		answer = solver.Solve(&ComputeTestResidual3d, &ComputeTestJacobian3d,
 							  residual, initial_guess, NULL);
    	
    	PetscScalar *answer_elements;
		VecGetArray(answer, &answer_elements);
        int lo, hi;
        VecGetOwnershipRange(answer, &lo, &hi);
        
        double tol = 1e-6;
        
        for(int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if(lo<=global_index && global_index<hi)
            {    
                double x = answer_elements[local_index];
                TS_ASSERT_DELTA(x,1/sqrt(3),tol);                
            }
        }        
        
        VecRestoreArray(answer,&answer_elements);
    	
    	VecDestroy(initial_guess);
    	VecDestroy(residual);
    	VecDestroy(answer);
	}

};

PetscErrorCode ComputeTestResidual(SNES snes,Vec solution_guess,Vec residual,void *pContext)
{
	double x,y;
	
	PetscScalar *solution_guess_elements;
	VecGetArray(solution_guess, &solution_guess_elements);
    
    int lo, hi;
	VecGetOwnershipRange(solution_guess, &lo, &hi);
    
    double all_solution_guess[2], all_solution_guess_replicated[2];
    all_solution_guess[0]=0.0;
    all_solution_guess[1]=0.0;
        
        
    for(int global_index=0; global_index<2; global_index++)
    {
        int local_index = global_index-lo;
        if(lo<=global_index && global_index<hi)
        { 
            all_solution_guess[global_index]=solution_guess_elements[local_index];
        }
    } 



    VecRestoreArray(solution_guess,&solution_guess_elements);
    
    MPI_Allreduce(all_solution_guess, all_solution_guess_replicated, 2, MPI_DOUBLE,
                             MPI_SUM, PETSC_COMM_WORLD);
    
    x = all_solution_guess_replicated[0];
	y = all_solution_guess_replicated[1];
	
	VecSetValue(residual,0,x*x+y*y-1,INSERT_VALUES);
	VecSetValue(residual,1,x-y,INSERT_VALUES);
	VecAssemblyBegin(residual);
	VecAssemblyEnd(residual);
	return 0;
}

PetscErrorCode ComputeTestJacobian(SNES snes,Vec input,Mat *pJacobian ,Mat *pPreconditioner,MatStructure *pMatStructure ,void *pContext)
{
	double x, y;
		
	PetscScalar *input_elements;
	VecGetArray(input, &input_elements);
	int lo, hi;
    VecGetOwnershipRange(input, &lo, &hi);
    
    double all_input[2], all_input_replicated[2];
    all_input[0]=0.0;
    all_input[1]=0.0;
    
    for(int global_index=0; global_index<2; global_index++)
    {
        int local_index = global_index-lo;
        if(lo<=global_index && global_index<hi)
        {    
            all_input[global_index]=input_elements[local_index];
        }
    }    
    
    
    VecRestoreArray(input,&input_elements);
    MPI_Allreduce(all_input, all_input_replicated, 2, MPI_DOUBLE,
                             MPI_SUM, PETSC_COMM_WORLD);
    x = all_input_replicated[0];
	y = all_input_replicated[1];
		
	MatSetValue(*pJacobian, 0 , 0 , 2.0*x , INSERT_VALUES);
	MatSetValue(*pJacobian, 0 , 1 , 2.0*y, INSERT_VALUES);
	MatSetValue(*pJacobian, 1 , 0 , 1.0, INSERT_VALUES);
	MatSetValue(*pJacobian, 1 , 1 , -1.0, INSERT_VALUES);
    MatAssemblyBegin(*pJacobian,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*pJacobian,MAT_FINAL_ASSEMBLY);
	
	return 0;
}

PetscErrorCode ComputeTestResidual3d(SNES snes,Vec solution_guess,Vec residual,void *pContext)
{
	double x,y,z;
	
	PetscScalar *solution_guess_elements;
	VecGetArray(solution_guess, &solution_guess_elements);
    int lo, hi;
    VecGetOwnershipRange(solution_guess, &lo, &hi);
    
    double all_solution_guess[3], all_solution_guess_replicated[3];
    all_solution_guess[0]=0.0;
    all_solution_guess[1]=0.0;
    all_solution_guess[2]=0.0;
    

    for(int global_index=0; global_index<3; global_index++)
    {
        int local_index = global_index-lo;
        if(lo<=global_index && global_index<hi)
        {    
            all_solution_guess[global_index]=solution_guess_elements[local_index];
        }
    }  

    VecRestoreArray(solution_guess,&solution_guess_elements);
    
    MPI_Allreduce(all_solution_guess, all_solution_guess_replicated, 3, MPI_DOUBLE,
                             MPI_SUM, PETSC_COMM_WORLD);
    
    x = all_solution_guess_replicated[0];
    y = all_solution_guess_replicated[1];
    z = all_solution_guess_replicated[2];
	
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
		
	PetscScalar *input_elements;
	VecGetArray(input, &input_elements);
    int lo, hi;
    VecGetOwnershipRange(input, &lo, &hi);
    
    double all_input[3], all_input_replicated[3];
    all_input[0]=0.0;
    all_input[1]=0.0;
    all_input[2]=0.0;
    

    for(int global_index=0; global_index<3; global_index++)
    {
        int local_index = global_index-lo;
        if(lo<=global_index && global_index<hi)
        {    
            all_input[global_index]=input_elements[local_index];
        }
    }  

    
    VecRestoreArray(input,&input_elements);
    MPI_Allreduce(all_input, all_input_replicated, 3, MPI_DOUBLE,
                             MPI_SUM, PETSC_COMM_WORLD);
    x = all_input_replicated[0];
    y = all_input_replicated[1];
    z = all_input_replicated[2];
		
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
