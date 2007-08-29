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
#include "PetscTools.hpp"

#define ONEDIM_FLAGGEDBCC FlaggedMeshBoundaryConditionsContainer<1,1>

class TestFlaggedMeshBoundaryConditionsContainer : public CxxTest::TestSuite
{
public:
    // a test where the bcs are added to the bcc explicity
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
        
        double value = 32.23432;
        
        // Boundary conditions - zero dirichlet at first and last node of flagged region
        FlaggedMeshBoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition =
            new ConstBoundaryCondition<1>(value);
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

        if(PetscTools::IsSequential())
        {
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
            TS_ASSERT_DELTA( linear_system.GetRhsVectorElement(0), value, 1e-9);
            TS_ASSERT_DELTA( linear_system.GetRhsVectorElement(1), 2.0, 1e-9);
            TS_ASSERT_DELTA( linear_system.GetRhsVectorElement(2), value, 1e-9);
        } 
    }

    // a test where the bcc computes the bcs (for a mixed mesh problem)
    void TestInterpolateBoundaryConditionsFromCourseToFine()
    {
        ConformingTetrahedralMesh<3,3> fine_mesh;
        
        fine_mesh.ConstructCuboid(6, 6, 6);
        double sixth=1.0L/6.0L;
        fine_mesh.Scale(sixth, sixth, sixth);
        
        // create coarse mesh as RTM
        MixedTetrahedralMesh<3,3> coarse_mesh;
        
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
        std::vector<double> soln(coarse_mesh.GetNumNodes());
        for (unsigned i=0; i<coarse_mesh.GetNumNodes(); i++)
        {
            c_vector<double,3> posn=coarse_mesh.GetNode(i)->rGetLocation();
            soln[i] = posn[0] + 2*posn[1] - posn[2];
        }
        Vec solution_vector = PetscTools::CreateVec(soln);
        
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
        
        VecDestroy(solution_vector);
    }

    
    // a test where the bcc computes the bcs (for a mixed mesh problem), including extrapolation
    void TestInterpolateAndExtrapolateBoundaryConditionsFromCourseToFine()
    {
        TrianglesMeshReader<2,2> fine_mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructFromMeshReader(fine_mesh_reader);
        
        TrianglesMeshReader<2,2> coarse_mesh_reader("mesh/test/data/DecimatedDisk");
        MixedTetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructFromMeshReader(coarse_mesh_reader);
        
        coarse_mesh.SetFineMesh(&fine_mesh);
        
        // Flag the right semicircle of the coarse mesh
        ConformingTetrahedralMesh<2, 2>::ElementIterator i_coarse_element;
        for (i_coarse_element = coarse_mesh.GetElementIteratorBegin();
             i_coarse_element != coarse_mesh.GetElementIteratorEnd();
             i_coarse_element++)
        {
            Element<2,2> &element = **i_coarse_element;
            ChastePoint<2> centroid = ChastePoint<2>(element.CalculateCentroid());
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
        std::vector<double> soln(coarse_mesh.GetNumNodes());
        for (unsigned i=0; i<coarse_mesh.GetNumNodes(); i++)
        {
            c_vector<double,3> posn=coarse_mesh.GetNode(i)->rGetLocation();
            soln[i] = 5*posn[0] + 7*posn[1];
        }
        Vec solution_vector = PetscTools::CreateVec(soln);
        
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
        
        VecDestroy(solution_vector);
    }
    
    void TestConstantDirichletConstructor()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        typedef FlaggedMeshBoundaryConditionsContainer<1,1> FlaggedMeshBccOneDim;
        TS_ASSERT_THROWS_ANYTHING(FlaggedMeshBccOneDim bad_bcc(mesh, 1.0));
        
        // flag elements 0-3, therefore nodes 0 and 4 are the boundary 
        mesh.GetElement(0)->Flag();
        mesh.GetElement(1)->Flag();
        mesh.GetElement(2)->Flag();
        mesh.GetElement(3)->Flag();

        double value = 23.2342;
        FlaggedMeshBoundaryConditionsContainer<1,1> bcc(mesh, value);
        
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
        map[0]=0; // node 4 = row 0, say
        map[4]=2; // node 8 = row 2, say

        bcc.ApplyDirichletToLinearProblem(linear_system, map, true);
    
        linear_system.AssembleFinalLinearSystem();

        if(PetscTools::IsSequential())
        {
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
            TS_ASSERT_DELTA( linear_system.GetRhsVectorElement(0), value, 1e-9);
            TS_ASSERT_DELTA( linear_system.GetRhsVectorElement(1), 2.0, 1e-9);
            TS_ASSERT_DELTA( linear_system.GetRhsVectorElement(2), value, 1e-9);
        } 
    }
};


#endif /*TESTFLAGGEDMESHBOUNDARYCONDITIONSCONTAINER_HPP_*/
