#ifndef TESTELASTICITYBOUNDARYCONDITIONSCONTAINER_HPP_
#define TESTELASTICITYBOUNDARYCONDITIONSCONTAINER_HPP_

#include <cxxtest/TestSuite.h>

#include "ElasticityBoundaryConditionsContainer.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "PetscSetupAndFinalize.hpp"

class TestBoundaryConditionsContainer : public CxxTest::TestSuite
{
public:

    void test_FixNode_and_SetDisplacement_1d() throw(Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        ElasticityBoundaryConditionsContainer<1> bcc(mesh.GetNumNodes());
        
        Node<1>* p_node_0 = mesh.GetNodeAt(0);
        Node<1>* p_node_1 = mesh.GetNodeAt(10);
        
        bcc.FixNode(p_node_0);
        
        TS_ASSERT_DELTA( bcc.GetDirichletBCValue(p_node_0,0), 0.0, 1e-10 );
        
        c_vector<double,1> displacement;
        displacement(0) = 3.14;
        
        bcc.SetDisplacement(p_node_1, displacement);
        
        TS_ASSERT_DELTA( bcc.GetDirichletBCValue(p_node_1,0), 3.14, 1e-10 );
    }
    
    
    void test_FixNode_and_SetDisplacement_2d() throw(Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        ElasticityBoundaryConditionsContainer<2> bcc(mesh.GetNumNodes());
        
        Node<2>* p_node_0 = mesh.GetNodeAt(0);
        Node<2>* p_node_1 = mesh.GetNodeAt(1);
        
        bcc.FixNode(p_node_0);
        
        TS_ASSERT_DELTA( bcc.GetDirichletBCValue(p_node_0,0), 0.0, 1e-10 );
        TS_ASSERT_DELTA( bcc.GetDirichletBCValue(p_node_0,1), 0.0, 1e-10 );
        
        c_vector<double,2> displacement;
        displacement(0) = 3.14;
        displacement(1) = 2.81;
        
        bcc.SetDisplacement(p_node_1, displacement);
        
        TS_ASSERT_DELTA( bcc.GetDirichletBCValue(p_node_1,0), 3.14, 1e-10 );
        TS_ASSERT_DELTA( bcc.GetDirichletBCValue(p_node_1,1), 2.81, 1e-10 );
    }
    
    
    void test_FixNode_and_SetDisplacement_3d() throw(Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        ElasticityBoundaryConditionsContainer<3> bcc(mesh.GetNumNodes());
        
        Node<3>* p_node_0 = mesh.GetNodeAt(0);
        Node<3>* p_node_1 = mesh.GetNodeAt(1);
        
        bcc.FixNode(p_node_0);
        
        TS_ASSERT_DELTA( bcc.GetDirichletBCValue(p_node_0,0), 0.0, 1e-10 );
        TS_ASSERT_DELTA( bcc.GetDirichletBCValue(p_node_0,1), 0.0, 1e-10 );
        TS_ASSERT_DELTA( bcc.GetDirichletBCValue(p_node_0,2), 0.0, 1e-10 );
        
        c_vector<double,3> displacement;
        displacement(0) = 3.14;
        displacement(1) = 2.81;
        displacement(2) = -1;
        
        bcc.SetDisplacement(p_node_1, displacement);
        
        TS_ASSERT_DELTA( bcc.GetDirichletBCValue(p_node_1,0),  3.14, 1e-10 );
        TS_ASSERT_DELTA( bcc.GetDirichletBCValue(p_node_1,1),  2.81, 1e-10 );
        TS_ASSERT_DELTA( bcc.GetDirichletBCValue(p_node_1,2), -1.00, 1e-10 );
    }
    
};
#endif /*TESTELASTICITYBOUNDARYCONDITIONSCONTAINER_HPP_*/
