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
#ifndef _TESTQUADRATICMESH_HPP_
#define _TESTQUADRATICMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "QuadraticMesh.hpp"

class TestQuadraticMesh : public CxxTest::TestSuite 
{
public:
    void TestQuadraticMesh1d() throw(Exception)
    {
        QuadraticMesh<1> mesh("mesh/test/data/1D_0_to_1_10_elements_quadratic");
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 10u);

        // node 2 (ie middle) of element 0 
        TS_ASSERT_EQUALS(mesh.GetElementNode(0, 2), 11u);
        TS_ASSERT_DELTA(mesh.GetNode(11)->rGetLocation()[0], 0.05, 1e-12);
        
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            bool is_boundary_node = mesh.GetNode(i)->IsBoundaryNode();
            TS_ASSERT_EQUALS(is_boundary_node, ((x==0)||(x==1)));
        }
    }
    
    void TestQuadraticMesh2d() throw(Exception)
    {
        QuadraticMesh<2> mesh("mesh/test/data/square_128_elements_quadratic");
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 289u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 128u);

        // node 3 (ie fourth) of element 0 
        TS_ASSERT_EQUALS(mesh.GetElementNode(0, 3), 82u);
        // node 4 (ie fifth) of element 0 
        TS_ASSERT_EQUALS(mesh.GetElementNode(0, 4), 83u);
        // node 5 (ie last) of element 0 
        TS_ASSERT_EQUALS(mesh.GetElementNode(0, 5), 81u);

        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            bool is_boundary_node = mesh.GetNode(i)->IsBoundaryNode();
            TS_ASSERT_EQUALS(is_boundary_node, ((x==0)||(x==1)||(y==0)||(y==1)));
        }
    }
    
    void TestQuadraticMesh3d() throw(Exception)
    {
        QuadraticMesh<3> mesh("mesh/test/data/3D_Single_tetrahedron_element_quadratic");
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);

        // check getting global numbers of nodes 4-9 (in non-vertices)
        for(unsigned i=4; i<10; i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElementNode(0, i), i);
        }

        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetNode(i)->IsBoundaryNode(), true);
        }
  
        QuadraticMesh<3> mesh2("mesh/test/data/cube_136_elements_quadratic");
        TS_ASSERT_EQUALS(mesh2.GetNumNodes(), 285u);
        TS_ASSERT_EQUALS(mesh2.GetNumElements(), 136u);

        for(unsigned i=0; i<mesh2.GetNumNodes(); i++)
        {
            double x = mesh2.GetNode(i)->rGetLocation()[0];
            double y = mesh2.GetNode(i)->rGetLocation()[1];
            double z = mesh2.GetNode(i)->rGetLocation()[2];
            bool is_boundary_node = mesh2.GetNode(i)->IsBoundaryNode();
            TS_ASSERT_EQUALS(is_boundary_node,  ((x==0)||(x==1)||(y==0)||(y==1)||(z==0)||(z==1)));
        }
    }
};

#endif // _TESTQUADRATICMESH_HPP_
