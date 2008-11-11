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
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ParallelTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

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
    
        mesh_reader.Reset();        
        TetrahedralMesh<2,2> seq_mesh;
        seq_mesh.ConstructFromMeshReader(mesh_reader);
        
        Element<2,2>* p_para_element;
        Element<2,2>* p_sequ_element;
        unsigned element_index;
        
        for (ParallelTetrahedralMesh<2,2>::ElementIterator it=mesh.GetElementIteratorBegin(); 
             it!=mesh.GetElementIteratorEnd(); 
             ++it)
        {
            p_para_element = *it;
            element_index = p_para_element->GetIndex();
            
            p_sequ_element = mesh.GetElement(element_index);            
            TS_ASSERT_EQUALS(element_index, p_sequ_element->GetIndex());
            
            for (unsigned node_local_index=0; node_local_index < p_para_element->GetNumNodes(); node_local_index++)
            {
                TS_ASSERT_EQUALS(p_para_element->GetNodeGlobalIndex(node_local_index), 
                                 p_sequ_element->GetNodeGlobalIndex(node_local_index));                                 

                TS_ASSERT_EQUALS(p_para_element->GetNode(node_local_index)->GetPoint()[0], 
                                 p_sequ_element->GetNode(node_local_index)->GetPoint()[0]);                                 
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

            
};
#endif /*TESTPARALLELTETRAHEDRALMESH_HPP_*/
