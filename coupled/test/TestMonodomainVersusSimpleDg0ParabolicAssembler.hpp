#ifndef _TESTMONODOMAINVERSUSSIMPLEDG0PARABOLICASSEMBLER_HPP_
#define _TESTMONODOMAINVERSUSSIMPLEDG0PARABOLICASSEMBLER_HPP_

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
#include "AbstractLinearParabolicPde.hpp"

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
             if (global_index==10) TS_ASSERT_DELTA(res1[global_index-lo], 5.8028e-07, 1e-9);
            if (global_index==25) TS_ASSERT_DELTA(res1[global_index-lo], 0.00648079, 1e-5);
            if (global_index==50) TS_ASSERT_DELTA(res1[global_index-lo], 0.992718, 1e-5);
            if (global_index==75) TS_ASSERT_DELTA(res1[global_index-lo], 0.00648079, 1e-5);
        }
  
  
        VecRestoreArray(current_solution_1, &res1);
        VecRestoreArray(current_solution_2, &res2);
        VecDestroy(current_solution_1);
        VecDestroy(current_solution_2);
    }
    void TestMonodomainDg0AssemblerWithFischer2DAgainstSimpleDg0Assembler()
    {
        double tStart = 0;
        double tFinal = 1;
        double tBigStep = 0.01;
        // Create mesh from mesh reader 
        TrianglesMeshReader mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        FischerPde<2> pde;
        
        // Boundary conditions: zero neumann on entire boundary (2 elements)
        BoundaryConditionsContainer<2,2> bcc(1, mesh.GetNumNodes());
        ConstBoundaryCondition<2>* p_neumann_boundary_condition = new ConstBoundaryCondition<2>(0.0);
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetFirstBoundaryElement();
       
        while(iter < mesh.GetLastBoundaryElement())
        {
            bcc.AddNeumannBoundaryCondition(*iter, p_neumann_boundary_condition);
            iter++;
        }

        // Linear solver
        SimpleLinearSolver linearSolver;

        // Assembler
        MonodomainDg0Assembler<2,2> monodomain_assembler;
        SimpleDg0ParabolicAssembler<2,2> simple_assembler;
        
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
//            if (global_index==10) TS_ASSERT_DELTA(res1[global_index-lo], 2.8951e-7, 1e-9);
//            if (global_index==25) TS_ASSERT_DELTA(res1[global_index-lo], 0.0060696, 1e-5);
//            if (global_index==50) TS_ASSERT_DELTA(res1[global_index-lo], 0.992834, 1e-5);
//            if (global_index==75) TS_ASSERT_DELTA(res1[global_index-lo], 0.0060696, 1e-5);
        }
        VecRestoreArray(current_solution_1, &res1);
        VecRestoreArray(current_solution_2, &res2);
        VecDestroy(current_solution_1);
        VecDestroy(current_solution_2);
    }
};

#endif //_TESTMONODOMAINVERSUSSIMPLEDG0PARABOLICASSEMBLER_HPP_
