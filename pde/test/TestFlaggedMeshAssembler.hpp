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

        // flag four middle elements         
        ConformingTetrahedralMesh<1,1>::ElementIterator iter
           = mesh.GetElementIteratorBegin();
          
        while(iter!=mesh.GetElementIteratorEnd())
        {
            if((4<=(*iter)->GetIndex()) && ((*iter)->GetIndex()<=7))
            {
                (*iter)->Flag();
            }
            
            iter++;
        }


        // Instantiate PDE object
        TimeDependentDiffusionEquationPde<1> pde;
        
        // Boundary conditions - zero dirichlet at first and last node;
        BoundaryConditionsContainer<1,1,1> bcc(mesh.GetNumNodes());
        ConstBoundaryCondition<1>* p_boundary_condition =
            new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode( mesh.GetNumNodes()-1 ), p_boundary_condition);
        
        // Assembler
        FlaggedMeshAssembler<1> assembler(&mesh,&pde,&bcc);
// need this?
        assembler.SetMatrixIsConstant();

// change this from 11 to 5        
        Vec initial_condition = CreateConstantConditionVec(11,0.0);
        
        assembler.AssembleSystem(initial_condition,0.0);
        
//        unsigned size_of_linear_system = assembler.mpLinearSystem->GetSize();
        
//        TS_ASSERT_EQUALS(size_of_linear_system,5u); 
        
        std::cout << "smasrm vector:\n";
        assembler.mpLinearSystem->DisplayRhs();
        std::cout << "smasrm matrix:\n";
        assembler.mpLinearSystem->DisplayMatrix();                                  
                               
    }
};
#endif /*TESTFLAGGEDMESHASSEMBLER_HPP_*/
