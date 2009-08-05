/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef TESTMIXEDTETRAHEDRALMESH_HPP_
#define TESTMIXEDTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>

#include "MixedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestMixedTetrahedralMesh : public CxxTest::TestSuite
{
private:
    double GetMinOrMax(TetrahedralMesh<2,2>* pMesh, bool min, unsigned component)
    {
        assert(component<2);
        double ret = min ? DBL_MAX : -DBL_MAX;
        for(unsigned i=0; i<pMesh->GetNumNodes(); i++)
        {
            double val = pMesh->GetNode(i)->rGetLocation()[component];
            if( (min && (val<ret)) || (!min && (val>ret)) )
            {
                ret = val;
            }
        }
        return ret;
    }

public:
//   void TestMeshConstructionFromMeshReader(void)
// is not yet implemented

    void TestConstructionFromCuboidMeshes3D()
    {
        // create fine mesh as CTM

        TetrahedralMesh<3,3> fine_mesh;

        fine_mesh.ConstructCuboid(6, 6, 6);
        double sixth=1.0L/6.0L;
        fine_mesh.Scale(sixth, sixth, sixth);

        // create coarse mesh as RTM
        MixedTetrahedralMesh<3,3> coarse_mesh;

        coarse_mesh.ConstructCuboid(3, 3, 3);
        double third=1.0L/3.0L;
        coarse_mesh.Scale(third, third, third);

        // give fine mesh to coarse mesh and calculate node map
        coarse_mesh.SetFineMesh(&fine_mesh);

        TS_ASSERT_EQUALS(coarse_mesh.GetFineMesh(), &fine_mesh);

        const NodeMap &node_map = coarse_mesh.rGetCoarseFineNodeMap();

        TS_ASSERT_EQUALS(node_map.GetNewIndex(0), 0u);
        TS_ASSERT_EQUALS(node_map.GetNewIndex(1), 2u);
        TS_ASSERT_EQUALS(node_map.GetNewIndex(4), 14u);
        TS_ASSERT_EQUALS(node_map.GetNewIndex(63), 342u);
        //Top node is 4^3-1 and 7^3-1 respectively

        TS_ASSERT_EQUALS(coarse_mesh.GetFineNodeIndexForCoarseNode(0), 0u);
        TS_ASSERT_EQUALS(coarse_mesh.GetFineNodeIndexForCoarseNode(1), 2u);
        TS_ASSERT_EQUALS(coarse_mesh.GetFineNodeIndexForCoarseNode(4), 14u);
        TS_ASSERT_EQUALS(coarse_mesh.GetFineNodeIndexForCoarseNode(63), 342u);

        // We're not allowed to call SetFineMesh twice
        TS_ASSERT_THROWS_ANYTHING(coarse_mesh.SetFineMesh(&fine_mesh));
    }

//    void TestCoarseFineElementsMap2D(void)
//    {
//        TetrahedralMesh<2,2> fine_mesh;
//        fine_mesh.ConstructRectangularMesh(2, 2, false);
//        double half=1.0L/2.0L;
//        fine_mesh.Scale(half, half, 0.0);
//
//        // create coarse mesh as RTM
//        MixedTetrahedralMesh<2,2> coarse_mesh;
//        coarse_mesh.ConstructRectangularMesh(1, 1, false);
//
//        // give fine mesh to coarse mesh and calculate node map
//        coarse_mesh.SetFineMesh(&fine_mesh);
//
//        std::set<Element<2,2>*> expected_elements;
//        expected_elements.insert(fine_mesh.GetElement(6));
//        expected_elements.insert(fine_mesh.GetElement(0));
//        expected_elements.insert(fine_mesh.GetElement(3));
//        expected_elements.insert(fine_mesh.GetElement(2));
//        // Elements 0,2,3 and 6 form the upper right half of the space
//        // (added in a funny order since set equality ought to cope with this)
//
//        TS_ASSERT(expected_elements ==
//                  coarse_mesh.GetFineElementsForCoarseElementIndex(0));
//    }


    void TestFineNodesCoarseElementsMap2D(void)
    {
        TetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructRectangularMesh(2, 2, false);
        double half=1.0L/2.0L;
        fine_mesh.Scale(half, half, 0.0);

        // create coarse mesh as RTM
        MixedTetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructRectangularMesh(1, 1, false);

        // give fine mesh to coarse mesh and calculate node map
        coarse_mesh.SetFineMesh(&fine_mesh);

        //node 1 is on the top edge of the fine mesh
        TS_ASSERT(coarse_mesh.GetElement(0) ==
                  coarse_mesh.GetACoarseElementForFineNodeIndex(1));
        //node 3 is on the left edge of the fine mesh
        TS_ASSERT(coarse_mesh.GetElement(1) ==
                  coarse_mesh.GetACoarseElementForFineNodeIndex(3));
    }

    void TestFineMeshIncorrect3D(void)
    {
        TetrahedralMesh<3,3> fine_mesh;

        fine_mesh.ConstructCuboid(5, 5, 5);
        double fifth=1.0L/5.0L;
        fine_mesh.Scale(fifth, fifth, fifth);

        // create coarse mesh as RTM
        MixedTetrahedralMesh<3,3> coarse_mesh;

        coarse_mesh.ConstructCuboid(3, 3, 3);
        double third=1.0L/3.0L;
        coarse_mesh.Scale(third, third, third);

        // give fine mesh to coarse mesh and calculate node map
        // should throw because not every coarse node has a coincident fine node
        TS_ASSERT_THROWS_ANYTHING(coarse_mesh.SetFineMesh(&fine_mesh));
    }

    void TestFineAndCoarseDisc(void)
    {
        //Covers cases where there are fine nodes outside the convex hull of the coarse mesh
        TrianglesMeshReader<2,2> fine_mesh_reader("mesh/test/data/disk_984_elements");
        TetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructFromMeshReader(fine_mesh_reader);

        TrianglesMeshReader<2,2> coarse_mesh_reader("projects/PhilTransA2009/test/data/DecimatedDisk");
        MixedTetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructFromMeshReader(coarse_mesh_reader);

        coarse_mesh.SetFineMesh(&fine_mesh);
    }

    void TestConstructRectangularMeshes()
    {
        MixedTetrahedralMesh<2,2> mixed_mesh;
        mixed_mesh.ConstructRectangularMeshes(2.2, 1.3, 4, 12);

        TS_ASSERT_EQUALS(mixed_mesh.GetNumNodes(), 25u);
        TS_ASSERT_EQUALS(mixed_mesh.GetFineMesh()->GetNumNodes(), 169u);

        TS_ASSERT_DELTA(GetMinOrMax(&mixed_mesh, true,  0), 0.0, 1e-8);
        TS_ASSERT_DELTA(GetMinOrMax(&mixed_mesh, false, 0), 2.2, 1e-8);
        TS_ASSERT_DELTA(GetMinOrMax(&mixed_mesh, true,  1), 0.0, 1e-8);
        TS_ASSERT_DELTA(GetMinOrMax(&mixed_mesh, false, 1), 1.3, 1e-8);

        TS_ASSERT_DELTA(GetMinOrMax(mixed_mesh.GetFineMesh(), true,  0), 0.0, 1e-8);
        TS_ASSERT_DELTA(GetMinOrMax(mixed_mesh.GetFineMesh(), false, 0), 2.2, 1e-8);
        TS_ASSERT_DELTA(GetMinOrMax(mixed_mesh.GetFineMesh(), true,  1), 0.0, 1e-8);
        TS_ASSERT_DELTA(GetMinOrMax(mixed_mesh.GetFineMesh(), false, 1), 1.3, 1e-8);
    }
};





#endif /*TESTMIXEDTETRAHEDRALMESH_HPP_*/
