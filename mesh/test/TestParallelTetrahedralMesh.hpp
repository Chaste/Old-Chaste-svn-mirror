/*

Copyright (C) University of Oxford, 2008

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

#ifndef TESTPARALLELTETRAHEDRALMESH_HPP_
#define TESTPARALLELTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "PetscSetupAndFinalize.hpp"
#include "ParallelTetrahedralMesh.hpp"

class TestParallelTetrahedralMesh : public CxxTest::TestSuite
{
    
public:

    void TestSimple()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        
        ParallelTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 543U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984U);

        unsigned num_local_nodes = mesh.GetNumLocalNodes();

        unsigned nodes_reduction;        

        MPI_Allreduce(&num_local_nodes, &nodes_reduction, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);        

        TS_ASSERT_EQUALS(nodes_reduction, mesh.GetNumNodes());        
        
        // Check some node co-ordinates
//        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0],  0.9980267283, 1e-6);
//        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], -0.0627905195, 1e-6);
//        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 1.0, 1e-6);
//        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 0.0, 1e-6);
//
//        // Check first element has the right nodes
//        TetrahedralMesh<2,2>::ElementIterator it = mesh.GetElementIteratorBegin();
//        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(0), 309U);
//        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(1), 144U);
//        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(2), 310U);
//        TS_ASSERT_EQUALS((*it)->GetNode(1), mesh.GetNode(144));
    }    
    
};
#endif /*TESTPARALLELTETRAHEDRALMESH_HPP_*/
