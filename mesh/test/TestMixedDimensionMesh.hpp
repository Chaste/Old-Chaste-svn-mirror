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
        MixedDimensionMesh<2,2> mesh; ///\todo Make a dumb partition 
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
        else
        {
            TS_ASSERT_LESS_THAN(mesh.GetNumLocalCableElements(), 11u);
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
