#ifndef TESTFLAGGEDMESHASSEMBLER_HPP_
#define TESTFLAGGEDMESHASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petsc.h>
#include <vector>
#include <cmath>
#include <iostream>
#include "TrianglesMeshReader.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "FlaggedMeshAssembler.hpp"
#include "TimeDependentDiffusionEquationPde.hpp"
//#include "ConstBoundaryCondition.hpp"

class TestFlaggedMeshAssembler : public CxxTest::TestSuite
{
private :
    Vec CreateInitialConditionVec(int size)
    {
        Vec initial_condition;
        VecCreate(PETSC_COMM_WORLD, &initial_condition);
        VecSetSizes(initial_condition, PETSC_DECIDE, size);
        VecSetFromOptions(initial_condition);
        return initial_condition;
    }
    
    Vec CreateConstantConditionVec(int size, double value)
    {
        Vec initial_condition = CreateInitialConditionVec(size);
        
#if (PETSC_VERSION_MINOR == 2) //Old API
        VecSet(&value, initial_condition);
#else
        VecSet(initial_condition, value);
#endif
        
        VecAssemblyBegin(initial_condition);
        VecAssemblyEnd(initial_condition);
        return initial_condition;
    }
    
    
    
public :
    void TestAssembleSystem() throw(Exception)
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Flag four middle elements
        ConformingTetrahedralMesh<1,1>::ElementIterator iter
        = mesh.GetElementIteratorBegin();
        
        while (iter!=mesh.GetElementIteratorEnd())
        {
            if ((4<=(*iter)->GetIndex()) && ((*iter)->GetIndex()<=7))
            {
                (*iter)->Flag();
                TS_ASSERT((*iter)->GetOwnership());
                TS_ASSERT((*iter)->IsFlagged());
            }
            
            iter++;
        }
        
        // Instantiate PDE object
        TimeDependentDiffusionEquationPde<1> pde;
        
        // Boundary conditions - zero dirichlet at first and last node of flagged region
        FlaggedMeshBoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition =
            new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(4), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(8), p_boundary_condition);
        
        // Assembler
        FlaggedMeshAssembler<1> assembler(&mesh,&pde,&bcc);
        assembler.SetTimes(0.0, 1.0, 0.01);
        
        const size_t full_size = 11u;
        const size_t smasrm_size = 5u;
        Vec initial_condition = CreateConstantConditionVec(full_size, 0.0);
        
        assembler.AssembleSystem(initial_condition, 0.0);
        
        const unsigned size_of_linear_system = assembler.mpLinearSystem->GetSize();
        TS_ASSERT_EQUALS(size_of_linear_system, smasrm_size);
        
        std::cout << "smasrm vector:\n";
        assembler.mpLinearSystem->DisplayRhs();
        std::cout << "smasrm matrix:\n";
        assembler.mpLinearSystem->DisplayMatrix();
        
        
        // Test SMASRM values
        // It should have -8.3.., 26.6.., -8.3.. on the centre diagonals, except for the boundaries
        // These numbers come from looking at the matrix of an equivalent problem
        // using the simple dg0 parabolic assembler
        PetscInt lo, hi;
        assembler.mpLinearSystem->GetOwnershipRange(lo, hi);
        for (unsigned i=1; i<smasrm_size-1; i++)
        {
            if ((unsigned)lo <= i && i < (unsigned)hi)
            {
                TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(i, i-1), -8.3333333333, 1e-8);
                TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(i, i),   26.6666666666, 1e-8);
                TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(i, i+1), -8.3333333333, 1e-8);
            }
        }
        
        // Dirichlet boundaries
        unsigned i=0;
        if ((unsigned)lo <= i && i < (unsigned)hi)
        {
            TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(0, 0), 1.0, 1e-8);
            TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(0, 1), 0.0, 1e-8);
        }
        i=4;
        if ((unsigned)lo <= i && i < (unsigned)hi)
        {
            TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(4, 3), 0.0, 1e-8);
            TS_ASSERT_DELTA(assembler.mpLinearSystem->GetMatrixElement(4, 4), 1.0, 1e-8);
        }
        
        // Check that if we add a boundary condition to an unflagged node
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(3), p_boundary_condition);
        // and then assemble the system
        TS_ASSERT_THROWS_ANYTHING(assembler.AssembleSystem(initial_condition, 0.0));
    }
};
#endif /*TESTFLAGGEDMESHASSEMBLER_HPP_*/
