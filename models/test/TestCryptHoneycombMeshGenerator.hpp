#ifndef TESTCRYPTHONECOMBMESHGENERATOR_HPP_
#define TESTCRYPTHONECOMBMESHGENERATOR_HPP_
#include <cxxtest/TestSuite.h>
#include "CryptHoneycombMeshGenerator.hpp"

class TestCryptHoneycombMeshGenerator : public CxxTest::TestSuite
{
public:
    void testCryptHoneycombMeshGenerator() throw(Exception)
    {
        int num_cells_width = 8;
        int num_cells_depth = 22;
        CryptHoneycombMeshGenerator generator(num_cells_width, num_cells_depth);
        
        // check the mesh
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),12*26);

        // zeroth node
        TS_ASSERT_DELTA(p_mesh->GetNodeAt(0)->GetPoint()[0],-2.0, 1e-12); 
        TS_ASSERT_DELTA(p_mesh->GetNodeAt(0)->GetPoint()[1],-sqrt(3)/2,1e-8); 
        
        // first real node
        int index = num_cells_width+4+2; // 4 here is the number of ghost nodes in a row
        TS_ASSERT_DELTA(p_mesh->GetNodeAt(index)->GetPoint()[0], 0.5,1e-12); 
        TS_ASSERT_DELTA(p_mesh->GetNodeAt(index)->GetPoint()[1], 0.0,1e-12); 

        // last real node
        index = num_cells_depth*(num_cells_width+4)+9;
        TS_ASSERT_DELTA(p_mesh->GetNodeAt(index)->GetPoint()[0], 7.0,1e-12); 
        TS_ASSERT_DELTA(p_mesh->GetNodeAt(index)->GetPoint()[1], 21.0*sqrt(3)/2.0,1e-4); 

        // last node
        int last_node = p_mesh->GetNumNodes()-1;
        TS_ASSERT_DELTA(p_mesh->GetNodeAt(last_node)->GetPoint()[0], 9.5,1e-12); 
        TS_ASSERT_DELTA(p_mesh->GetNodeAt(last_node)->GetPoint()[1], 24.0*sqrt(3)/2.0,1e-6); 

        // check the ghost nodes
        std::vector<int> ghost_node_indices = generator.GetGhostNodeIndices();
        for(int i=0; i<13; i++)
        {
           TS_ASSERT(ghost_node_indices[i]==i);
        }
        
        delete p_mesh;
    } 
};


#endif /*TESTCRYPTHONECOMBMESHGENERATOR_HPP_*/
