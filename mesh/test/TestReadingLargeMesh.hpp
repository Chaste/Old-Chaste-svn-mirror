/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TESTREADINGLARGEMESH_HPP_
#define TESTREADINGLARGEMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"

#include <vector>

class TestReadingLargeConformingTetrahedralMesh : public CxxTest::TestSuite
{
public:
    // This test is mainly here for performance testing, to check that loading a
    // (relatively) large mesh doesn't take too long.
    // It's a nightly test because it takes 10 hours to run under MemoryTesting!
    void TestLoadingLargeMesh(void)
    {
        TrianglesMeshReader<3,3> meshReader("heart/test/data/heart");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(meshReader);
        
        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 63885U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 322267U);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 41812U);
    }
};

#endif /*TESTREADINGLARGEMESH_HPP_*/
