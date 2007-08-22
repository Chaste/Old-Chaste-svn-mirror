#ifndef TESTFLAGGEDMESHBOUNDARYCONDITIONSCONTAINER_HPP_
#define TESTFLAGGEDMESHBOUNDARYCONDITIONSCONTAINER_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petsc.h>
#include <vector>
#include <cmath>
#include <iostream>
#include "TrianglesMeshReader.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "FlaggedMeshBoundaryConditionsContainer.hpp"

class TestFlaggedMeshBoundaryConditionsContainer : public CxxTest::TestSuite
{
public:
    void TestWithDefaultConstructor()
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
        
        // Boundary conditions - zero dirichlet at first and last node of flagged region
        FlaggedMeshBoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition =
            new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(4), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(8), p_boundary_condition);
        
        LinearSystem linear_system(3);
        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                // LHS matrix is all 3s
                linear_system.SetMatrixElement(i,j,3);
            }
            // RHS vector is all 2s
            linear_system.SetRhsVectorElement(i,2);
        }
        
        linear_system.AssembleIntermediateLinearSystem();

        std::map<unsigned, unsigned> map;
        map[4]=0; // node 4 = row 0, say
        map[8]=2; // node 8 = row 2, say

        bcc.ApplyDirichletToLinearProblem(linear_system, map, true);
    
        linear_system.AssembleFinalLinearSystem();

        // altered row
        TS_ASSERT_DELTA( linear_system.GetMatrixElement(0,0), 1.0, 1e-9);
        TS_ASSERT_DELTA( linear_system.GetMatrixElement(0,1), 0.0, 1e-9);
        TS_ASSERT_DELTA( linear_system.GetMatrixElement(0,2), 0.0, 1e-9);

        // unaltered row
        TS_ASSERT_DELTA( linear_system.GetMatrixElement(1,0), 3.0, 1e-9);
        TS_ASSERT_DELTA( linear_system.GetMatrixElement(1,1), 3.0, 1e-9);
        TS_ASSERT_DELTA( linear_system.GetMatrixElement(1,2), 3.0, 1e-9);
        
        // altered row
        TS_ASSERT_DELTA( linear_system.GetMatrixElement(2,0), 0.0, 1e-9);
        TS_ASSERT_DELTA( linear_system.GetMatrixElement(2,1), 0.0, 1e-9);
        TS_ASSERT_DELTA( linear_system.GetMatrixElement(2,2), 1.0, 1e-9);
    
        // vector
        TS_ASSERT_DELTA( linear_system.GetRhsVectorElement(0), 0.0, 1e-9);
        TS_ASSERT_DELTA( linear_system.GetRhsVectorElement(1), 2.0, 1e-9);
        TS_ASSERT_DELTA( linear_system.GetRhsVectorElement(2), 0.0, 1e-9);
    } 
    
//    
//    void TestWithMixedMeshConstructor()
//    {
//        ConformingTetrahedralMesh<2,2> fine_mesh;
//        
//        unsigned num_elem = 48;
//        fine_mesh.ConstructRectangularMesh(num_elem,num_elem);
//        fine_mesh.Scale(1.0/num_elem, 1.0/num_elem);
//        
//        // create coarse mesh as RTM
//        RefinedTetrahedralMesh<2,2> coarse_mesh;
//        
//        num_elem = 12;
//        coarse_mesh.ConstructRectangularMesh(num_elem, num_elem);
//        coarse_mesh.Scale(1.0/num_elem, 1.0/num_elem);
//        
//        // give fine mesh to coarse mesh
//        coarse_mesh.SetFineMesh(&fine_mesh);
//        
//        
//        Vec initial_condition_coarse = CreateInitialConditionVec(coarse_mesh.GetNumNodes());
//        
//        double* p_initial_condition_coarse;
//        VecGetArray(initial_condition_coarse, &p_initial_condition_coarse);
//        
//        int lo, hi;
//        VecGetOwnershipRange(initial_condition_coarse, &lo, &hi);
//        
//        for (int global_index = lo; global_index < hi; global_index++)
//        {
//            int local_index = global_index - lo;
//            double x = coarse_mesh.GetNode(global_index)->GetPoint()[0];
//            double y = coarse_mesh.GetNode(global_index)->GetPoint()[1];
//            if(sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)) < 0.3)
//            {
//                p_initial_condition_coarse[local_index] = 1.0;
//            }
//            else
//            {
//                p_initial_condition_coarse[local_index] = 0.0;
//            }
//        }
//        VecRestoreArray(initial_condition_coarse, &p_initial_condition_coarse);
//                
//            Vec result = assembler.Solve();
//
//            // Reset problem size to fine mesh - will have been changed by assembler to match coarse mesh
//            DistributedVector::SetProblemSize(fine_mesh.GetNumNodes());
//            
//            coarse_mesh.UnflagAllElements();
//            
//            // Flag the same elements of the coarse mesh each time step, for now
//            ConformingTetrahedralMesh<2, 2>::ElementIterator i_coarse_element;
//            for (i_coarse_element = coarse_mesh.GetElementIteratorBegin();
//                 i_coarse_element != coarse_mesh.GetElementIteratorEnd();
//                 i_coarse_element++)
//            {
//                Element<2,2> &element = **i_coarse_element;
//                for(unsigned i=0; i<element.GetNumNodes(); i++)
//                {
//                    if(result_replicated[element.GetNodeGlobalIndex(i)]>0.4)
//                    {
//                        element.Flag();
//                    }
//                }
//            }
//
//            coarse_mesh.TransferFlags();
//            FlaggedMeshBoundaryConditionsContainer<2,1> flagged_bcc(coarse_mesh, result);
//        
//        
//    }
};


#endif /*TESTFLAGGEDMESHBOUNDARYCONDITIONSCONTAINER_HPP_*/
