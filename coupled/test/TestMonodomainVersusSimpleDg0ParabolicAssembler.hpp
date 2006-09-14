#ifndef _TESTMONODOMAINVERSUSSIMPLEDG0PARABOLICASSEMBLER_HPP_
#define _TESTMONODOMAINVERSUSSIMPLEDG0PARABOLICASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petsc.h>
#include <vector>

#include "SimpleLinearSolver.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "AbstractLinearParabolicPde.hpp"


template <int SPACE_DIM>
class ZeroStimCellFactory : public AbstractCardiacCellFactory<SPACE_DIM>
{
public:
    ZeroStimCellFactory() : AbstractCardiacCellFactory<SPACE_DIM>(0.01)
    {
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        return new LuoRudyIModel1991OdeSystem(this->mpSolver, this->mTimeStep, this->mpZeroStimulus);
    }
    
    virtual unsigned GetNumberOfCells()
    {
        return 1;
    }
};




/*
 * A simple parabolic PDE used in this test.
 * 
 * NOTE: The fischer pde is a parabolic linear pde, however we want to get the 
 * MonodomainDg0Assembler, which expects a MonodomainPde passed it, to solve it.
 * Therefore, we get this pde to inherit from MonodomainPde, and overload all the
 * main functions. For this reason it has to take in a cell class, although this is
 * not used.
 */
template <int SPACE_DIM>
class FischerPde : public MonodomainPde<SPACE_DIM>
{
public:
    double ComputeLinearSourceTerm(Point<SPACE_DIM> )
    {
        return 0.0;
    }
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> , double u)
    {
        return u*(1-u);
    }
    
    double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& , double u)
    {
        return u*(1-u);
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(Point<SPACE_DIM> )
    {
        return identity_matrix<double>(SPACE_DIM);
    }
    
    double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> )
    {
        return 1;
    }
    
    FischerPde() : MonodomainPde<SPACE_DIM>(new ZeroStimCellFactory<SPACE_DIM>)
    {
    }    
};




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
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/heart_FHN_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        FischerPde<1> pde;
        
        // Boundary conditions: zero neumann on entire boundary (2 elements)
        BoundaryConditionsContainer<1,1,1> bcc(mesh.GetNumNodes());
        ConstBoundaryCondition<1>* p_zero_condition = new ConstBoundaryCondition<1>(0.0);
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc.AddNeumannBoundaryCondition(*iter, p_zero_condition);
            iter++;
        }
        
        // Assembler
        MonodomainDg0Assembler<1,1> monodomain_assembler(&mesh,&pde);
        SimpleDg0ParabolicAssembler<1,1> simple_assembler(&mesh,&pde,&bcc);
        
        // initial condition;
        Vec initial_condition_1, initial_condition_2;
        initial_condition_1 = CreateInitialConditionVec(mesh.GetNumNodes());
        VecDuplicate(initial_condition_1, &initial_condition_2);
        
        double* p_init_array;
        int lo, hi;
        VecGetOwnershipRange(initial_condition_1, &lo, &hi);
        VecGetArray(initial_condition_1, &p_init_array);
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index-lo;
            double x=mesh.GetNodeAt(global_index)->GetPoint()[0];
            p_init_array[local_index] = exp(-(x*x)/100);
        }
        VecRestoreArray(initial_condition_1, &p_init_array);
        VecAssemblyBegin(initial_condition_1);
        VecAssemblyEnd(initial_condition_1);
        
        VecCopy(initial_condition_1, initial_condition_2); // Both assemblers use same initial cond'n
        
        // Vars to hold current solutions at each iteration
        Vec current_solution_1, current_solution_2;

        double tCurrent = tStart;
        while ( tCurrent < tFinal )
        {
            monodomain_assembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            simple_assembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            
            monodomain_assembler.SetInitialCondition( initial_condition_1 );
            simple_assembler.SetInitialCondition( initial_condition_2 );
            
            current_solution_1 = monodomain_assembler.Solve();

            current_solution_2 = simple_assembler.Solve();
            
            // Next iteration uses current solution as initial condition
            VecDestroy(initial_condition_1); // Free old initial condition
            VecDestroy(initial_condition_2);
            initial_condition_1 = current_solution_1;
            initial_condition_2 = current_solution_2;
            
            tCurrent += tBigStep;
        }
        
        // Compare the results
        double *p_current_solution1_array, *p_current_solution2_array;
        VecGetArray(current_solution_1, &p_current_solution1_array);
        VecGetArray(current_solution_2, &p_current_solution2_array);
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index-lo;
            TS_ASSERT_DELTA(p_current_solution1_array[global_index-lo], p_current_solution2_array[local_index], 1e-3);
            if (global_index==10) TS_ASSERT_DELTA(p_current_solution1_array[local_index], 5.8028e-07, 1e-9);
            if (global_index==25) TS_ASSERT_DELTA(p_current_solution1_array[local_index], 0.00648079, 1e-5);
            if (global_index==50) TS_ASSERT_DELTA(p_current_solution1_array[local_index], 0.992718, 1e-5);
            if (global_index==75) TS_ASSERT_DELTA(p_current_solution1_array[local_index], 0.00648079, 1e-5);
        }
        
        
        VecRestoreArray(current_solution_1, &p_current_solution1_array);
        VecRestoreArray(current_solution_2, &p_current_solution2_array);
        VecDestroy(current_solution_1);
        VecDestroy(current_solution_2);
    }
    void TestMonodomainDg0AssemblerWithFischer2DAgainstSimpleDg0Assembler()
    {
        double tStart = 0;
        double tFinal = 1;
        double tBigStep = 0.01;
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        FischerPde<2> pde;
        
        // Boundary conditions: zero neumann on entire boundary (2 elements)
        BoundaryConditionsContainer<2,2,1> bcc(mesh.GetNumNodes());
        ConstBoundaryCondition<2>* p_neumann_boundary_condition = new ConstBoundaryCondition<2>(0.0);
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc.AddNeumannBoundaryCondition(*iter, p_neumann_boundary_condition);
            iter++;
        }
        
        // Assembler
        MonodomainDg0Assembler<2,2> monodomain_assembler(&mesh,&pde);
        SimpleDg0ParabolicAssembler<2,2> simple_assembler(&mesh,&pde,&bcc);
                
        // initial condition;
        Vec initial_condition_1, initial_condition_2;
        initial_condition_1 = CreateInitialConditionVec(mesh.GetNumNodes());
        VecDuplicate(initial_condition_1, &initial_condition_2);
        
        double* p_init_array;
        int lo, hi;
        VecGetOwnershipRange(initial_condition_1, &lo, &hi);
        VecGetArray(initial_condition_1, &p_init_array);
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index-lo;
            double x=mesh.GetNodeAt(global_index)->GetPoint()[0];
            p_init_array[local_index] = exp(-(x*x)/100);
        }
        VecRestoreArray(initial_condition_1, &p_init_array);
        VecAssemblyBegin(initial_condition_1);
        VecAssemblyEnd(initial_condition_1);
        
        VecCopy(initial_condition_1, initial_condition_2); // Both assemblers use same initial cond'n
        
        // Vars to hold current solutions at each iteration
        Vec current_solution_1, current_solution_2;
        
        double tCurrent = tStart;
        while ( tCurrent < tFinal )
        {
            monodomain_assembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            simple_assembler.SetTimes(tCurrent, tCurrent+tBigStep, tBigStep);
            
            monodomain_assembler.SetInitialCondition( initial_condition_1 );
            simple_assembler.SetInitialCondition( initial_condition_2 );
            
            current_solution_1 = monodomain_assembler.Solve();
            
            current_solution_2 = simple_assembler.Solve();
            
            // Next iteration uses current solution as initial condition
            VecDestroy(initial_condition_1); // Free old initial condition
            VecDestroy(initial_condition_2);
            initial_condition_1 = current_solution_1;
            initial_condition_2 = current_solution_2;
            
            tCurrent += tBigStep;
        }
        
        // Compare the results
        double *p_current_solution1_array, *p_current_solution2_array;
        VecGetArray(current_solution_1, &p_current_solution1_array);
        VecGetArray(current_solution_2, &p_current_solution2_array);
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index-lo;
            TS_ASSERT_DELTA(p_current_solution1_array[local_index], p_current_solution2_array[local_index], 1e-3);
//            if (global_index==10) TS_ASSERT_DELTA(p_current_solution1_array[local_index], 2.8951e-7, 1e-9);
//            if (global_index==25) TS_ASSERT_DELTA(p_current_solution1_array[local_index], 0.0060696, 1e-5);
//            if (global_index==50) TS_ASSERT_DELTA(p_current_solution1_array[local_index], 0.992834, 1e-5);
//            if (global_index==75) TS_ASSERT_DELTA(p_current_solution1_array[local_index], 0.0060696, 1e-5);
        }
        VecRestoreArray(current_solution_1, &p_current_solution1_array);
        VecRestoreArray(current_solution_2, &p_current_solution2_array);
        VecDestroy(current_solution_1);
        VecDestroy(current_solution_2);
    }
};

#endif //_TESTMONODOMAINVERSUSSIMPLEDG0PARABOLICASSEMBLER_HPP_
