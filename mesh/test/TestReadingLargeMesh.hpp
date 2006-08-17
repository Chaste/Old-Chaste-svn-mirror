#ifndef TESTREADINGLARGEMESH_HPP_
#define TESTREADINGLARGEMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"

#include <vector>

class TestConformingTetrahedralMesh : public CxxTest::TestSuite
{
public:
    // This test is mainly here for performance testing, to check that loading a
    // (relatively) large mesh doesn't take too long.
    // It's a nightly test because it takes 10 hours to run under MemoryTesting!
    void TestLoadingLargeMesh(void)
    {
        TrianglesMeshReader<3,3> meshReader("mesh/test/data/heart");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(meshReader, 1);
        
        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 63885);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 322267);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 41812);
    }
};

#endif /*TESTREADINGLARGEMESH_HPP_*/
