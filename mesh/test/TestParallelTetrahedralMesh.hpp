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

#ifndef TESTPARALLELTETRAHEDRALMESH_HPP_
#define TESTPARALLELTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ParallelTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

class TestParallelTetrahedralMesh : public CxxTest::TestSuite
{
    
public:

    void TestConstructFromMeshReader1D()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_with_attributes");
        
        ParallelTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 11U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 10U);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1U);
        for (unsigned i=0; i<10; i++)
        {
            try
            {
                unsigned region = mesh.GetElement(i)->GetRegion();
                TS_ASSERT_EQUALS(region, i%5+1);
            }
            catch(Exception& e)
            {
                // I don't own this element                
            }
        }
    }

    void TestConstructFromMeshReader2D()
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
    
        mesh_reader.Reset();        
        TetrahedralMesh<2,2> seq_mesh;
        seq_mesh.ConstructFromMeshReader(mesh_reader);
        
        for (ParallelTetrahedralMesh<2,2>::ElementIterator it=mesh.GetElementIteratorBegin(); 
             it!=mesh.GetElementIteratorEnd(); 
             ++it)
        {
            Element<2,2>* p_para_element = *it;
            unsigned element_index = p_para_element->GetIndex();
            
            Element<2,2>* p_sequ_element = seq_mesh.GetElement(element_index);            
            TS_ASSERT_EQUALS(element_index, p_sequ_element->GetIndex());
            
            for (unsigned node_local_index=0; node_local_index < p_para_element->GetNumNodes(); node_local_index++)
            {
                TS_ASSERT_EQUALS(p_para_element->GetNodeGlobalIndex(node_local_index), 
                                 p_sequ_element->GetNodeGlobalIndex(node_local_index));                                 

                TS_ASSERT_EQUALS(p_para_element->GetNode(node_local_index)->GetPoint()[0], 
                                 p_sequ_element->GetNode(node_local_index)->GetPoint()[0]);                                 
            }
        }

        for (ParallelTetrahedralMesh<2,2>::BoundaryElementIterator it=mesh.GetBoundaryElementIteratorBegin(); 
             it!=mesh.GetBoundaryElementIteratorEnd(); 
             ++it)
        {
            BoundaryElement<1,2>* p_para_boundary_element = *it;
            unsigned boundary_element_index = p_para_boundary_element->GetIndex();
            
            BoundaryElement<1,2>* p_sequ_boundary_element = seq_mesh.GetBoundaryElement(boundary_element_index);            
            TS_ASSERT_EQUALS(boundary_element_index, p_sequ_boundary_element->GetIndex());
            
            for (unsigned node_local_index=0; node_local_index < p_para_boundary_element->GetNumNodes(); node_local_index++)
            {
                TS_ASSERT_EQUALS(p_para_boundary_element->GetNodeGlobalIndex(node_local_index), 
                                 p_sequ_boundary_element->GetNodeGlobalIndex(node_local_index));                                 

                TS_ASSERT_EQUALS(p_para_boundary_element->GetNode(node_local_index)->GetPoint()[0], 
                                 p_sequ_boundary_element->GetNode(node_local_index)->GetPoint()[0]);                                 
            }
            
            
        }
         
    }

//    void TestAllNodesAssigned()
//    {
//        
//    }
//    
//    void TestAllElementsAssigned()
//    {
//        
//    }    

    void TestConstructFromMeshReader3D()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 51U);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 136U);
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 96U);

        ParallelTetrahedralMesh<3,3> mesh;        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 51U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 136U);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 96U);

        try
        {
            ChastePoint<3> coords = mesh.GetNode(0)->GetPoint();
            TS_ASSERT_DELTA(coords[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(coords[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(coords[2], 0.0, 1e-6);
        }
        catch(Exception& e)
        {
            //I'm not the owner of node 0
        }

        try
        {
            ChastePoint<3> coords = mesh.GetNode(19)->GetPoint();            
            TS_ASSERT_DELTA(coords[0], 0.75, 1e-6);
            TS_ASSERT_DELTA(coords[1], 0.25, 1e-6);
            TS_ASSERT_DELTA(coords[2], 0.0, 1e-6);
        }
        catch(Exception& e)
        {
            //I'm not the owner of node 19
        }
    }
};
#endif /*TESTPARALLELTETRAHEDRALMESH_HPP_*/
