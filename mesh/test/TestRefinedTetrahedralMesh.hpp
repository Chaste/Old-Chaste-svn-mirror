#ifndef TESTREFINEDTETRAHEDRALMESH_HPP_
#define TESTREFINEDTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>

#include "RefinedTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"


class TestRefinedTetrahedralMesh : public CxxTest::TestSuite
{
public:

    void TestMeshConstructionFromMeshReader(void)
    {
        // create fine mesh as CTM
        
        ConformingTetrahedralMesh<3,3> fine_mesh;
        
        fine_mesh.ConstructCuboid(6, 6, 6);
        double sixth=1.0L/6.0L;
        fine_mesh.Scale(sixth, sixth, sixth);
        
        // create coarse mesh as RTM
        RefinedTetrahedralMesh<3,3> coarse_mesh;
        
        coarse_mesh.ConstructCuboid(3, 3, 3);
        double third=1.0L/3.0L;
        coarse_mesh.Scale(third, third, third);
       
        // give fine mesh to coarse mesh and calculate node map
        coarse_mesh.SetFineMesh(&fine_mesh);
        
        NodeMap &node_map = coarse_mesh.rGetCoarseFineNodeMap();
        
        TS_ASSERT_EQUALS(node_map.GetNewIndex(0), 0u);
        TS_ASSERT_EQUALS(node_map.GetNewIndex(1), 2u);
        TS_ASSERT_EQUALS(node_map.GetNewIndex(4), 14u);
        
        TS_ASSERT_EQUALS(node_map.GetNewIndex(63), 342u);
        //Top node is 4^3-1 and 7^3-1 respectively
    }
        
    void TestCoarseFineElementsMap(void)
    {   
        ConformingTetrahedralMesh<2,2> fine_mesh;        
        fine_mesh.ConstructRectangularMesh(2, 2, false);
        double half=1.0L/2.0L;
        fine_mesh.Scale(half, half, 0.0);
        
        // create coarse mesh as RTM
        RefinedTetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructRectangularMesh(1, 1, false);
       
        // give fine mesh to coarse mesh and calculate node map
        coarse_mesh.SetFineMesh(&fine_mesh);
        
        std::set<Element<2,2>*> expected_elements;
        expected_elements.insert(fine_mesh.GetElement(6));
        expected_elements.insert(fine_mesh.GetElement(0));
        expected_elements.insert(fine_mesh.GetElement(3));
        expected_elements.insert(fine_mesh.GetElement(2));
        // Elements 0,2,3 and 6 form the upper right half of the space
        // (added in a funny order since set equality ought to cope with this)
        
        
        TS_ASSERT(expected_elements ==
                         coarse_mesh.GetFineElementsForCoarseElementIndex(0)); 
    }
    
    void TestFineMeshIncorrect(void)
    {
        ConformingTetrahedralMesh<3,3> fine_mesh;
        
        fine_mesh.ConstructCuboid(5, 5, 5);
        double fifth=1.0L/5.0L;
        fine_mesh.Scale(fifth, fifth, fifth);
        
        // create coarse mesh as RTM
        RefinedTetrahedralMesh<3,3> coarse_mesh;
        
        coarse_mesh.ConstructCuboid(3, 3, 3);
        double third=1.0L/3.0L;
        coarse_mesh.Scale(third, third, third);
       
        // give fine mesh to coarse mesh and calculate node map
        // should throw because not every coarse node has a coincident fine node
        TS_ASSERT_THROWS_ANYTHING(coarse_mesh.SetFineMesh(&fine_mesh));
    }
    
};





#endif /*TESTREFINEDTETRAHEDRALMESH_HPP_*/
