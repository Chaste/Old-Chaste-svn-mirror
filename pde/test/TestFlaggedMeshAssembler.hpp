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
#include "SimpleDg0ParabolicAssembler.hpp"
#include "TimeDependentDiffusionEquationPde.hpp"
#include "TimeDependentDiffusionEquationWithSourceTermPde.hpp"
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
        
        //std::cout << "smasrm vector:\n";
        //assembler.mpLinearSystem->DisplayRhs();
        //std::cout << "smasrm matrix:\n";
        //assembler.mpLinearSystem->DisplayMatrix();
        
        
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
    
    
    void TestInterpolateAndExtrapolateBoundaryConditionsFromCourseToFine()
    {
        TrianglesMeshReader<2,2> fine_mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructFromMeshReader(fine_mesh_reader,1);
        
        TrianglesMeshReader<2,2> coarse_mesh_reader("mesh/test/data/DecimatedDisk");
        RefinedTetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructFromMeshReader(coarse_mesh_reader,1);
        
        
        coarse_mesh.SetFineMesh(&fine_mesh);
        
        // Flag the right semicircle of the coarse mesh
        ConformingTetrahedralMesh<2, 2>::ElementIterator i_coarse_element;
        for (i_coarse_element = coarse_mesh.GetElementIteratorBegin();
             i_coarse_element != coarse_mesh.GetElementIteratorEnd();
             i_coarse_element++)
        {
            Element<2,2> &element = **i_coarse_element;
            Point<2> centroid = Point<2>(element.CalculateCentroid());
            if (centroid[0] > 0)
            {
                element.Flag();
            }
            else
            {
                element.Unflag();
            }
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
            
            c_vector<double,2> posn=coarse_mesh.GetNode(global_index)->rGetLocation();
            p_solution_vector[local_index] = 5*posn[0] + 7*posn[1];
        }

        VecRestoreArray(solution_vector, &p_solution_vector);
        VecAssemblyBegin(solution_vector);
        VecAssemblyEnd(solution_vector);
        
        // throws as flags not transfered to fine mesh yet       
        typedef FlaggedMeshBoundaryConditionsContainer<2,1> FlaggedMeshBccTwoDimOneUnknown;
        TS_ASSERT_THROWS_ANYTHING(FlaggedMeshBccTwoDimOneUnknown bad_bcc(coarse_mesh, solution_vector));
        
        // Flag the corresponding region of the fine mesh
        coarse_mesh.TransferFlags();
        
        
        
        // interpolate boundary conditions        
        FlaggedMeshBoundaryConditionsContainer<2,1> bcc(coarse_mesh, solution_vector);
        
        // get the boundary of the flagged region
        std::set<unsigned> boundary = fine_mesh.CalculateBoundaryOfFlaggedRegion();
        
        std::set<unsigned>::iterator it = boundary.begin();
        while(it!=boundary.end())
        {
            unsigned node_index = *it;

            // an assertion would fail here if there is no boundary condition for 
            // this node
            double value = bcc.GetDirichletBCValue(fine_mesh.GetNode(node_index));
            
            c_vector<double,2> posn=fine_mesh.GetNode(node_index)->rGetLocation();
            double true_value = 5*posn[0] + 7*posn[1];
            
            TS_ASSERT_DELTA(value, true_value, 1e-12); 
            it++;
        }
    }
    
    void TestCoarseAndFineDiffusion()
    {
        ConformingTetrahedralMesh<2,2> fine_mesh;
        
        unsigned num_elem = 48;
        fine_mesh.ConstructRectangularMesh(num_elem,num_elem);
        fine_mesh.Scale(1.0/num_elem, 1.0/num_elem);
        
        // create coarse mesh as RTM
        RefinedTetrahedralMesh<2,2> coarse_mesh;
        
        num_elem = 12;
        coarse_mesh.ConstructRectangularMesh(num_elem, num_elem);
        coarse_mesh.Scale(1.0/num_elem, 1.0/num_elem);
        
        // give fine mesh to coarse mesh
        coarse_mesh.SetFineMesh(&fine_mesh);        
     
        // Instantiate PDE object
        TimeDependentDiffusionEquationWithSourceTermPde<2> pde;
        
        // Boundary conditions - zero dirichlet everywhere on boundary
        BoundaryConditionsContainer<2,2,1> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&coarse_mesh);
        
        // Assembler
        SimpleDg0ParabolicAssembler<2,2> assembler(&coarse_mesh,&pde,&bcc);
        
        Vec initial_condition_coarse = CreateInitialConditionVec(coarse_mesh.GetNumNodes());
        
        double* p_initial_condition_coarse;
        VecGetArray(initial_condition_coarse, &p_initial_condition_coarse);
        
        int lo, hi;
        VecGetOwnershipRange(initial_condition_coarse, &lo, &hi);
        
        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            double x = coarse_mesh.GetNode(global_index)->GetPoint()[0];
            double y = coarse_mesh.GetNode(global_index)->GetPoint()[1];
            if(sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)) < 0.3)
            {
                p_initial_condition_coarse[local_index] = 1.0;
            }
            else
            {
                p_initial_condition_coarse[local_index] = 0.0;
            }
        }
        VecRestoreArray(initial_condition_coarse, &p_initial_condition_coarse);
        
        ReplicatableVector ic_coarse_replicated(initial_condition_coarse);
        
                
        // Create initial_condition_fine from initial_condition_coarse by interpolation        
        Vec initial_condition_fine = CreateInitialConditionVec(fine_mesh.GetNumNodes());
        DistributedVector::SetProblemSize(fine_mesh.GetNumNodes());
        coarse_mesh.InterpolateOnUnflaggedRegion(initial_condition_coarse, initial_condition_fine);
        
        // Solve on coarse mesh        
        assembler.SetTimes(0, 0.01, 0.01);
        assembler.SetInitialCondition(initial_condition_coarse);
        
        Vec result = assembler.Solve();  
        ReplicatableVector result_replicated(result);
        
        for(unsigned i=0; i<result_replicated.size(); i++)
        {
            double x = coarse_mesh.GetNode(i)->rGetLocation()[0];
            double y = coarse_mesh.GetNode(i)->rGetLocation()[1];
            std::cout << i << " " << x << " " << y << " " << result_replicated[i] << "\n";
        }
        TS_ASSERT_DELTA(result_replicated[14],  0.0117061, 1e-6);
        
        // Flag the some elements of the coarse mesh
        ConformingTetrahedralMesh<2, 2>::ElementIterator i_coarse_element;
        for (i_coarse_element = coarse_mesh.GetElementIteratorBegin();
             i_coarse_element != coarse_mesh.GetElementIteratorEnd();
             i_coarse_element++)
        {
            Element<2,2> &element = **i_coarse_element;
            for(unsigned i=0; i<element.GetNumNodes(); i++)
            {
                if(result_replicated[element.GetNodeGlobalIndex(i)]>0.4)
                {
                    element.Flag();
                    std::cout << "Flagging element " << element.GetIndex() << "\n";
                }
            }
        }
        
        // Flag the corresponding region of the fine mesh
        coarse_mesh.TransferFlags();
        
        // interpolate boundary conditions        
        FlaggedMeshBoundaryConditionsContainer<2,1> flagged_bcc(coarse_mesh, result);

        // Assembler
        FlaggedMeshAssembler<2> flagged_assembler(&fine_mesh,&pde,&flagged_bcc);
        flagged_assembler.SetTimes(0.0, 0.01, 0.01);
        flagged_assembler.SetInitialCondition(initial_condition_fine);

        flagged_assembler.Solve();
        
        Vec result_fine = flagged_assembler.Solve();  
        ReplicatableVector result_fine_replicated(result_fine);

        std::map<unsigned, unsigned> map = flagged_assembler.GetSmasrmIndexMap();
        for(unsigned i=0; i<fine_mesh.GetNumNodes(); i++)
        {
            std::map<unsigned, unsigned>::iterator iter = map.find(i);
            if(iter!=map.end())
            {
                unsigned smasrm_index = iter->second;

                double x = fine_mesh.GetNode(i)->rGetLocation()[0];
                double y = fine_mesh.GetNode(i)->rGetLocation()[1];
                
                std::cout << i << " " << x << " " << y << " " << result_fine_replicated[smasrm_index] << "\n";
                if (i == 1992)
                {
                    TS_ASSERT_DELTA(result_fine_replicated[smasrm_index], 0.170209, 1e-6);
                }
            }
        }
    }
};
#endif /*TESTFLAGGEDMESHASSEMBLER_HPP_*/
