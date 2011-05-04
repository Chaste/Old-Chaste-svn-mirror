/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTMIXEDDIMENSIONMESH_HPP_
#define TESTMIXEDDIMENSIONMESH_HPP_

#include <cxxtest/TestSuite.h>

#include "MixedDimensionMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "Exception.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestMixedDimensionMesh : public CxxTest::TestSuite
{
public:
    void TestReadingSquareMesh() throw (Exception)
    {
        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);
        MixedDimensionMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 121u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 200u);
        
        TS_ASSERT_EQUALS(mesh.GetNumCableElements(), 10u);
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(mesh.GetNumLocalCableElements(), 10u);
            
            for (unsigned i=0; i<10u; i++)
            {
                Element<1,2>* p_cable_elt = mesh.GetCableElement(i);
                TS_ASSERT_EQUALS(p_cable_elt->GetNumNodes(), 2u);
                TS_ASSERT_EQUALS(p_cable_elt->GetNodeGlobalIndex(0u), 55u + i);
                TS_ASSERT_EQUALS(p_cable_elt->GetNodeGlobalIndex(1u), 56u + i);
                TS_ASSERT_EQUALS(p_cable_elt->GetNode(0u), mesh.GetNode(55u + i));
                TS_ASSERT_EQUALS(p_cable_elt->GetNode(1u), mesh.GetNode(56u + i));
                TS_ASSERT_EQUALS(p_cable_elt->GetRegion(), i+1);
                TS_ASSERT( mesh.CalculateDesignatedOwnershipOfCableElement(i) );
            }
            
            for (unsigned i=0; i<200u; i++)
            {
                Element<2,2>* p_elt = mesh.GetElement(i);
                TS_ASSERT_EQUALS(p_elt->GetNumNodes(), 3u);
                TS_ASSERT_EQUALS(p_elt->GetNode(0u), mesh.GetNode(p_elt->GetNodeGlobalIndex(0u)));
                TS_ASSERT_EQUALS(p_elt->GetNode(1u), mesh.GetNode(p_elt->GetNodeGlobalIndex(1u)));
                TS_ASSERT_EQUALS(p_elt->GetNode(2u), mesh.GetNode(p_elt->GetNodeGlobalIndex(2u)));
            }
        }
        else if (PetscTools::GetNumProcs() == 2)
        {
            if (PetscTools::GetMyRank() == 0)
            {
                // For a dumb partition, process 0 owns nodes 0 to 59
                TS_ASSERT_EQUALS(mesh.GetNumLocalCableElements(), 5u);

                for (unsigned i=0; i<5u; i++)
                {
                    Element<1,2>* p_cable_elt = mesh.GetCableElement(i);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNumNodes(), 2u);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNodeGlobalIndex(0u), 55u + i);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNodeGlobalIndex(1u), 56u + i);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNode(0u), mesh.GetNode(55u + i));
                    TS_ASSERT_EQUALS(p_cable_elt->GetNode(1u), mesh.GetNodeOrHaloNode(56u + i));
                    TS_ASSERT_EQUALS(p_cable_elt->GetRegion(), i+1);
                    TS_ASSERT( mesh.CalculateDesignatedOwnershipOfCableElement(i) ); // Designated owner of all these five, since we own node 0 (lowest index)
                }
            }
            else
            {
                // For a dumb partition, process 1 owns nodes 60 to 120
                TS_ASSERT_EQUALS(mesh.GetNumLocalCableElements(), 6u);

                for (unsigned i=0; i<6u; i++)
                {
                    Element<1,2>* p_cable_elt = mesh.GetCableElement(i);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNumNodes(), 2u);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNodeGlobalIndex(0u), 59u + i);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNodeGlobalIndex(1u), 60u + i);
                    TS_ASSERT_EQUALS(p_cable_elt->GetNode(0u), mesh.GetNodeOrHaloNode(59u + i));
                    TS_ASSERT_EQUALS(p_cable_elt->GetNode(1u), mesh.GetNode(60u + i));
                    TS_ASSERT_EQUALS(p_cable_elt->GetRegion(), i+5);
                    
                    // Not designated owner of the first of these as node 0 is owned by process 0
                    if (i==0)
                    {
                        TS_ASSERT( ! mesh.CalculateDesignatedOwnershipOfCableElement(i) );
                    }
                    else
                    {
                        TS_ASSERT( mesh.CalculateDesignatedOwnershipOfCableElement(i) );
                    }
                        
                }
            }
        }
        else if (PetscTools::GetNumProcs() == 3)
        {
            if (PetscTools::GetMyRank() == 0)
            {
                // For a dumb partition, process 0 owns nodes 0 to 39
                TS_ASSERT_EQUALS(mesh.GetNumLocalCableElements(), 0u);
            }
            if (PetscTools::GetMyRank() == 1)
            {
                // For a dumb partition, process 1 owns nodes 40 to 79
                TS_ASSERT_EQUALS(mesh.GetNumLocalCableElements(), 10u);
            }
            if (PetscTools::GetMyRank() == 2)
            {
                // For a dumb partition, process 0 owns nodes 80 to 120
                TS_ASSERT_EQUALS(mesh.GetNumLocalCableElements(), 0u);
            }
        }
        else
        {
            TS_ASSERT_LESS_THAN(mesh.GetNumLocalCableElements(), 11u);
        }
    }
    
    void TestCableElementIterator() throw (Exception)
    {
        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);
        MixedDimensionMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        unsigned count = 0;
        for (MixedDimensionMesh<2,2>::CableElementIterator iter = mesh.GetCableElementIteratorBegin();
             iter != mesh.GetCableElementIteratorEnd();
             ++iter)
        {
            unsigned index = (*iter)->GetIndex();

            TS_ASSERT_EQUALS((*iter)->GetNumNodes(), 2u);
            TS_ASSERT_EQUALS((*iter)->GetNodeGlobalIndex(0u), 55u + index);
            TS_ASSERT_EQUALS((*iter)->GetNodeGlobalIndex(1u), 56u + index);
            TS_ASSERT_EQUALS((*iter)->GetNode(0u), mesh.GetNodeOrHaloNode(55u + index));
            TS_ASSERT_EQUALS((*iter)->GetNode(1u), mesh.GetNodeOrHaloNode(56u + index));
            TS_ASSERT_EQUALS((*iter)->GetRegion(), index+1);
            
            count++;
        }

        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(count, 10u);                
        }
        if (PetscTools::GetNumProcs() == 2)
        {            
            if (PetscTools::GetMyRank() == 0)
            {
                TS_ASSERT_EQUALS(count, 5u);
            }
            else
            {
                TS_ASSERT_EQUALS(count, 6u);
            }
        }
    }
    
    void TestReadingMeshWithNoCables() throw (Exception)
    {
        std::string mesh_base("mesh/test/data/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);
        MixedDimensionMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 121u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 200u);
        TS_ASSERT_EQUALS(mesh.GetNumCableElements(), 0u);
    }
    
    void TestExceptions() throw (Exception)
    {
        // Only TrianglesMeshReader supports cables
        MemfemMeshReader<3,3> memfem_reader("mesh/test/data/Memfem_slab");
        TS_ASSERT_EQUALS(memfem_reader.GetNumCableElements(), 0u);
        TS_ASSERT_EQUALS(memfem_reader.GetNumCableElementAttributes(), 0u);
        TS_ASSERT_THROWS_THIS(memfem_reader.GetNextCableElementData(), "Cable elements are not supported by this mesh format.");
        
        
    }
};

#endif /*TESTMIXEDDIMENSIONMESH_HPP_*/
