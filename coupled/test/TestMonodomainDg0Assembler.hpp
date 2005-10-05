#ifndef _TESTMONODOMAINDG0ASSEMBLER_HPP_
#define _TESTMONODOMAINDG0ASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>

#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Node.hpp"
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"  
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "ColumnDataWriter.hpp"

#include "MonodomainPde.hpp"
#include "MockEulerIvpOdeSolver.hpp"
#include "FischerPde.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"


class TestMonodomainDg0Assembler : public CxxTest::TestSuite 
{   
private:
	
	/**
	 * Refactor code to set up a PETSc vector holding the initial condition.
	 */
	Vec CreateInitialConditionVec(int size)
	{
    	Vec initial_condition;
    	VecCreate(PETSC_COMM_WORLD, &initial_condition);
    	VecSetSizes(initial_condition, PETSC_DECIDE, size);
    	VecSetFromOptions(initial_condition);
    	return initial_condition;
	}
    
public:
	void TestMonodomainDg0AssemblerWithFischer1DAgainstSimpleDg0Assembler()
	{
		
		double tStart = 0;
		double tFinal = 1;
		double tBigStep = 0.01;
		// Create mesh from mesh reader 
		TrianglesMeshReader mesh_reader("mesh/test/data/heart_FHN_mesh");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
        
		// Instantiate PDE object
		FischerPde<1> pde;
        
		// Boundary conditions: zero neumann on entire boundary (2 elements)
		BoundaryConditionsContainer<1,1> bcc(1, mesh.GetNumNodes());
		ConstBoundaryCondition<1>* pNeumannBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
		ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
		bcc.AddNeumannBoundaryCondition(*iter, pNeumannBoundaryCondition);
		iter++;
		bcc.AddNeumannBoundaryCondition(*iter, pNeumannBoundaryCondition);

		// Linear solver
		SimpleLinearSolver linearSolver;

		// Assembler
		MonodomainDg0Assembler<1,1> monodomain_assembler;
		SimpleDg0ParabolicAssembler<1,1> simple_assembler;
        
		// initial condition;   
		Vec initial_condition_1, initial_condition_2;
		initial_condition_1 = CreateInitialConditionVec(mesh.GetNumNodes());
		VecDuplicate(initial_condition_1, &initial_condition_2);
  
		double* init_array;
	    int lo, hi;
		VecGetOwnershipRange(initial_condition_1, &lo, &hi);
		int ierr = VecGetArray(initial_condition_1, &init_array); 
		for (int global_index=lo; global_index<hi; global_index++)
		{
			double x=mesh.GetNodeAt(global_index)->GetPoint()[0];
			init_array[global_index-lo] = exp(-(x*x)/100);
		}
		VecRestoreArray(initial_condition_1, &init_array);      
		VecAssemblyBegin(initial_condition_1);
		VecAssemblyEnd(initial_condition_1);
		
		VecCopy(initial_condition_1, initial_condition_2); // Both assemblers use same initial cond'n

		// Vars to hold current solutions at each iteration
		Vec current_solution_1, current_solution_2;

		double tCurrent = tStart;
		while( tCurrent < tFinal )
		{
			monodomain_assembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
			simple_assembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
			
			monodomain_assembler.SetInitialCondition( initial_condition_1 );
			simple_assembler.SetInitialCondition( initial_condition_2 );

			current_solution_1 = monodomain_assembler.Solve(mesh, &pde, bcc, &linearSolver);
            
			current_solution_2 = simple_assembler.Solve(mesh, &pde, bcc, &linearSolver);
			
     		// Next iteration uses current solution as initial condition
     		VecDestroy(initial_condition_1); // Free old initial condition
     		VecDestroy(initial_condition_2);
     		initial_condition_1 = current_solution_1;
     		initial_condition_2 = current_solution_2;
     
			tCurrent += tBigStep;
		}
		
		// Compare the results
		double *res1, *res2;
		VecGetArray(current_solution_1, &res1);
		VecGetArray(current_solution_2, &res2);
		for (int global_index=lo; global_index<hi; global_index++)
		{
			TS_ASSERT_DELTA(res1[global_index-lo], res2[global_index-lo], 1e-3);
            if (global_index==10) TS_ASSERT_DELTA(res1[global_index-lo], 2.8951e-7, 1e-9);
            if (global_index==25) TS_ASSERT_DELTA(res1[global_index-lo], 0.0060696, 1e-5);
            if (global_index==50) TS_ASSERT_DELTA(res1[global_index-lo], 0.992834, 1e-5);
            if (global_index==75) TS_ASSERT_DELTA(res1[global_index-lo], 0.0060696, 1e-5);
        }
		VecRestoreArray(current_solution_1, &res1);
		VecRestoreArray(current_solution_2, &res2);
		VecDestroy(current_solution_1);
		VecDestroy(current_solution_2);
	}
    
	void TestMonodomainDg01D()
	{
        MonodomainProblem<1> monodomainProblem("mesh/test/data/1D_0_to_1_100_elements",
                                               5, 
                                               -600.0, 
                                               "testoutput/MonoDg01d",
                                               "NewMonodomainLR91_1d");

        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.currentVoltage, &currentVoltageArray); 
        
        for(int global_index=monodomainProblem.lo; global_index<monodomainProblem.hi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[global_index-monodomainProblem.lo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[global_index-monodomainProblem.lo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomainProblem.monodomain_pde->GetOdeVarsAtNode(global_index);           
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }
               
                if (global_index==25)
                {
                    TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.lo], 16.972, 0.001);
                }
                if (global_index==30)
                {
                    TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.lo], 18.8826, 0.001);
                }
                if (global_index==40)
                {
                    TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.lo], 18.2916, 0.001);
                }
                if (global_index==50)
                {
                    TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.lo], -82.6575, 0.001);
                }
                if (global_index==51)
                {
                    TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.lo], -83.4065, 0.001);
                }
                if (global_index==75)
                {
                    TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.lo], -84.5504, 0.001);
                }
            }
        }
        VecRestoreArray(monodomainProblem.currentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.currentVoltage);
        VecAssemblyEnd(monodomainProblem.currentVoltage);
        VecDestroy(monodomainProblem.currentVoltage);
    }
    
 
    
    void testMonodomainDg02D( void )
    {   
        MonodomainProblem<2> monodomainProblem("mesh/test/data/square_128_elements",
                                               0.1, 
                                               -80.0, 
                                               "testoutput/MonoDg02d",
                                               "NewMonodomainLR91_2d");

        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.currentVoltage, &currentVoltageArray); 
        
        for(int global_index=monodomainProblem.lo; global_index<monodomainProblem.hi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[global_index-monodomainProblem.lo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[global_index-monodomainProblem.lo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomainProblem.monodomain_pde->GetOdeVarsAtNode(global_index);           
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }
            }
        }
        VecRestoreArray(monodomainProblem.currentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.currentVoltage);
        VecAssemblyEnd(monodomainProblem.currentVoltage);
        VecDestroy(monodomainProblem.currentVoltage);
    }   


    void testMonodomainDg03D( void )
    {
        MonodomainProblem<3> monodomainProblem("mesh/test/data/slab_138_elements",
                                               0.1, 
                                               -80.0, 
                                               "testoutput/MonoDg03d",
                                               "NewMonodomainLR91_3d");

        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.currentVoltage, &currentVoltageArray); 
        
        for(int global_index=monodomainProblem.lo; global_index<monodomainProblem.hi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[global_index-monodomainProblem.lo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[global_index-monodomainProblem.lo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomainProblem.monodomain_pde->GetOdeVarsAtNode(global_index);           
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }
            }
        }
        VecRestoreArray(monodomainProblem.currentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.currentVoltage);
        VecAssemblyEnd(monodomainProblem.currentVoltage);
        VecDestroy(monodomainProblem.currentVoltage);
    }   
};
#endif //_TESTMONODOMAINDG0ASSEMBLER_HPP_
