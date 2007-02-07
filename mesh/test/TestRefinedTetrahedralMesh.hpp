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
};





#endif /*TESTREFINEDTETRAHEDRALMESH_HPP_*/
