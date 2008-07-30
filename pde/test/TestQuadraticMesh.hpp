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
    void testQuadraticMesh() throw(Exception)
    {
        QuadraticMesh mesh("mesh/test/data/square_128_elements_quadratics");
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 289u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 128u);

        // node 3 (ie fourth) of element 0 
        TS_ASSERT_EQUALS(mesh.GetElementNode(0, 3), 82u);
        // node 4 (ie fifth) of element 0 
        TS_ASSERT_EQUALS(mesh.GetElementNode(0, 4), 83u);
        // node 5 (ie last) of element 0 
        TS_ASSERT_EQUALS(mesh.GetElementNode(0, 5), 81u);
    }
};

#endif // _TESTQUADRATICMESH_HPP_
