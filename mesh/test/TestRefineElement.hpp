#ifndef TESTREFINEELEMENT_HPP_
#define TESTREFINEELEMENT_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "Exception.hpp"


class TestRefineElement : public CxxTest::TestSuite
{
public:

    void Test1DRefineElement()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        
        ConformingTetrahedralMesh<1,1> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Point<1> new_point(0.01);
        Element<1,1>* p_first_element=mesh.GetElement(0);
        
        TS_ASSERT_THROWS_NOTHING(mesh.RefineElement(p_first_element,new_point));
        // Instead of a element with nodes at 0 and 0.1
        // there should be an element with nodes at 0 and 0.01 and
        // an element with nodes at 0.01 and 0.1. Other elements should stay the same
        
        const Node<1>* p_first_node=p_first_element->GetNode(0);
        const Node<1>* p_second_node=p_first_element->GetNode(1);
        
        TS_ASSERT_EQUALS(p_first_node->GetPoint().rGetLocation()(0), 0);
        TS_ASSERT_EQUALS(p_second_node->GetPoint().rGetLocation()(0), 0.01);
        
        // test second element
        
        Element<1,1> *p_second_element=mesh.GetElement(1);
        
        p_first_node=p_second_element->GetNode(0);
        p_second_node=p_second_element->GetNode(1);
        
        TS_ASSERT_EQUALS(p_first_node->GetPoint().rGetLocation()(0), 0.1);
        TS_ASSERT_EQUALS(p_second_node->GetPoint().rGetLocation()(0), 0.2);
        
        // test last element
        
        Element<1,1>* p_last_element=mesh.GetElement(10);
        
        p_first_node=p_last_element->GetNode(0);
        p_second_node=p_last_element->GetNode(1);
        
        TS_ASSERT_EQUALS(p_first_node->GetPoint().rGetLocation()(0), 0.01);
        TS_ASSERT_EQUALS(p_second_node->GetPoint().rGetLocation()(0), 0.1);
        
        // test jacobians
        TS_ASSERT_EQUALS(p_first_element->GetJacobianDeterminant(), 0.01);
        TS_ASSERT_EQUALS(p_second_element->GetJacobianDeterminant(), 0.1);
        TS_ASSERT_DELTA(p_last_element->GetJacobianDeterminant(), 0.09,1e-6);
        
        // test mesh length
        TS_ASSERT_EQUALS(mesh.CalculateMeshVolume(), 1);
    }
    
    void Test1DRefineElementWithPointTooFarRightFails() throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        
        ConformingTetrahedralMesh<1,1> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Point<1> new_point(0.11);
        Element<1,1>* p_first_element=mesh.GetElement(0);
        
        // Trying to add Point(0.11) to Element(0)
        // This point is contained in Element(1)
        TS_ASSERT_THROWS_ANYTHING(mesh.RefineElement(p_first_element, new_point));
        
    }
    
    void Test1DRefineElementWithPointTooFarLeftFails() throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        
        ConformingTetrahedralMesh<1,1> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Point<1> new_point(-0.1);
        Element<1,1>* p_first_element=mesh.GetElement(0);
        
        // Trying to add Point(-0.1) to Element(0)
        // This point is to the left of Element(0)
        TS_ASSERT_THROWS_ANYTHING(mesh.RefineElement(p_first_element, new_point));
        
    }
    
    void Test1DRefineElementManyNodes() throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        
        ConformingTetrahedralMesh<1,1> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Element<1,1>* p_first_element=mesh.GetElement(0);
        
        //There's space on the node vector for 10 new points
        // but more than 10 should still work
        for (int i=1; i<=20; i++)
        {
            Point<1> new_point(0.1 - i*0.0005);
            mesh.RefineElement(p_first_element,new_point);
        }
        
    }
    
    void Test2DRefineElement() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        
        ConformingTetrahedralMesh<2,2> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01,1e-6);
        TS_ASSERT_DELTA(mesh.CalculateMeshSurface(), 0.4,1e-6);
        
        
        // Refine an element in the bottom right corner
        Element<2,2>* p_corner_element=mesh.GetElement(18);
        
        TS_ASSERT_DELTA(p_corner_element->GetJacobianDeterminant(),0.0001 , 1e-6);
        
        // Point to be inserted in the bottom right corner
        Point<2> new_point(0.095,0.003);
        
        TS_ASSERT_THROWS_NOTHING(mesh.RefineElement(p_corner_element,new_point));
        
        //Testing the JacobianDeterminant of the element that has changed
        TS_ASSERT_DELTA(p_corner_element->GetJacobianDeterminant(),3e-05 , 1e-6);
        
        //Testing invariants
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01,1e-6);
        TS_ASSERT_DELTA(mesh.CalculateMeshSurface(), 0.4,1e-6);
        
        
        // Refine an element in the middle of the mesh
        Element<2,2>* p_middle_element=mesh.GetElement(108);
        
        TS_ASSERT_DELTA(p_middle_element->GetJacobianDeterminant(),0.0001 , 1e-6);
        
        // Point to be inserted in the middle of the mesh
        Point<2> new_point1(0.045,0.053);
        
        TS_ASSERT_THROWS_NOTHING(mesh.RefineElement(p_middle_element,new_point1));
        
        //Testing the JacobianDeterminant of the element that has changed
        TS_ASSERT_DELTA(p_middle_element->GetJacobianDeterminant(),3e-05 , 1e-6);
        
        //Testing invariants
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 0.01,1e-6);
        TS_ASSERT_DELTA(mesh.CalculateMeshSurface(), 0.4,1e-6);
    }
    
    void Test2DRefineElementFails() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        
        ConformingTetrahedralMesh<2,2> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Refine an element in the bottom right corner
        Element<2,2>* p_corner_element=mesh.GetElement(18);
        
        // Point to be inserted on the edge of the element
        Point<2> new_point(0.095,0.005);
        
        //Shouldn't be able to insert this point at the edge of the element
        TS_ASSERT_THROWS_ANYTHING(mesh.RefineElement(p_corner_element,new_point));
        
    }
    
    void Test3DRefineElement() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        
        ConformingTetrahedralMesh<3,3> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.CalculateMeshSurface(), 6.0, 1e-6);
        
        
        // Refine an element in the top corner (1, 1, 1)
        Element<3,3>* p_corner_element=mesh.GetElement(64);
        
        TS_ASSERT_DELTA(p_corner_element->GetJacobianDeterminant(),0.03125 , 1e-6);
        
        // Point to be inserted in the top corner
        Point<3> new_point(0.9,0.75,0.9);
        
        TS_ASSERT_THROWS_NOTHING(mesh.RefineElement(p_corner_element,new_point));
        
        //Testing the JacobianDeterminant of the element that has changed
        TS_ASSERT_DELTA(p_corner_element->GetJacobianDeterminant(), 0.0125 , 1e-6);
        
        //Testing invariants
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.CalculateMeshSurface(), 6.0, 1e-6);
        
        
        
        // Refine an element which includes the middle node
        Element<3,3>* p_middle_element=mesh.GetElement(49);
        
        TS_ASSERT_DELTA(p_middle_element->GetJacobianDeterminant(),0.0625 , 1e-6);
        
        // Point to be inserted near node 22 (middle of cube)
        Point<3> new_point1(0.49, 0.47, 0.6);
        
        TS_ASSERT_THROWS_NOTHING(mesh.RefineElement(p_middle_element, new_point1));
        
        //Testing the JacobianDeterminant of the element that has changed
        TS_ASSERT_DELTA(p_middle_element->GetJacobianDeterminant(), 0.01125 , 1e-6);
        
        //Testing invariants
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.CalculateMeshSurface(), 6.0, 1e-6);
    }
    
    void Test3DRefineElementFails() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        
        ConformingTetrahedralMesh<3,3> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.CalculateMeshSurface(), 6.0, 1e-6);
        
        
        // Refine an element in the top corner (1, 1, 1)
        Element<3,3>* p_corner_element=mesh.GetElement(64);
        
        TS_ASSERT_DELTA(p_corner_element->GetJacobianDeterminant(),0.03125 , 1e-6);
        
        // Point to be inserted in wrong place
        Point<3> new_point(0.9,0.75,1.0);
        
        TS_ASSERT_THROWS_ANYTHING(mesh.RefineElement(p_corner_element,new_point));
    }
    

};

#endif /*TESTREFINEELEMENT_HPP_*/
