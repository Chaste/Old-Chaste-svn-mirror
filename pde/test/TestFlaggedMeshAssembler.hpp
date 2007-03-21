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
#include "RefinedTetrahedralMesh.cpp"
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
    
    void TestInterpolateBoundaryConditionsFromCourseToFine()
    {
        ConformingTetrahedralMesh<3,3> fine_mesh;
        
        fine_mesh.ConstructCuboid(6, 6, 6);
        double sixth=1.0L/6.0L;
        fine_mesh.Scale(sixth, sixth, sixth);
        
        // create coarse mesh as RTM
        RefinedTetrahedralMesh<3,3> coarse_mesh;
        
        coarse_mesh.ConstructCuboid(3, 3, 3);
        double third=1.0L/3.0L;
        coarse_mesh.Scale(third, third, third);
        
        // give fine mesh to coarse mesh
        coarse_mesh.SetFineMesh(&fine_mesh);
        
        // Flag some elements
        ConformingTetrahedralMesh<3,3>::ElementIterator iter
        = fine_mesh.GetElementIteratorBegin();
        
        while (iter!=fine_mesh.GetElementIteratorEnd())
        {
            if((*iter)->CalculateCentroid()[0] > 0.5)
            {
                (*iter)->Flag();
            }
            iter++;
        }
        
        // set up petsc vector of the solution on the coarse mesh 
        unsigned num_coarse_nodes = coarse_mesh.GetNumNodes();
        Vec solution_vector;
        int lo, hi;
        VecCreate(PETSC_COMM_WORLD, &solution_vector);
        VecSetSizes(solution_vector, PETSC_DECIDE, num_coarse_nodes);
        VecSetFromOptions(solution_vector);
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        
        double *p_solution_vector;
        
        VecGetArray(solution_vector, &p_solution_vector);
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            
            c_vector<double,3> posn=coarse_mesh.GetNode(global_index)->rGetLocation();
            p_solution_vector[local_index] = posn[0] + 2*posn[1] - posn[2];
        }

        VecRestoreArray(solution_vector, &p_solution_vector);
        VecAssemblyBegin(solution_vector);
        VecAssemblyEnd(solution_vector);
        
        
        // interpolate boundary conditions        
        FlaggedMeshBoundaryConditionsContainer<3,1> bcc(coarse_mesh, solution_vector);
        
        // get the boundary of the flagged region
        std::set<unsigned> boundary = fine_mesh.CalculateBoundaryOfFlaggedRegion();
        
        std::set<unsigned>::iterator it = boundary.begin();
        while(it!=boundary.end())
        {
            unsigned node_index = *it;

            // an assertion would fail here if there is no boundary condition for 
            // this node
            double value = bcc.GetDirichletBCValue(fine_mesh.GetNode(node_index));
            
            c_vector<double,3> posn=fine_mesh.GetNode(node_index)->rGetLocation();
            double true_value = posn[0] + 2*posn[1] - posn[2];
            
            TS_ASSERT_DELTA(value, true_value, 1e-12); 
            it++;
        }
    }
};
#endif /*TESTFLAGGEDMESHASSEMBLER_HPP_*/
