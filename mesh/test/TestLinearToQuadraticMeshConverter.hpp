/*

Copyright (C) University of Oxford, 2005-2010

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

#ifndef TESTLINEARTOQUADRATICMESHCONVERTER_HPP_
#define TESTLINEARTOQUADRATICMESHCONVERTER_HPP_

#include <cxxtest/TestSuite.h>
#include "LinearToQuadraticMeshConverter.hpp"


class TestLinearToQuadraticMeshConverter : public CxxTest::TestSuite
{
public:
    void TestLinearToQuadraticMeshConverter2d() throw(Exception)
    {
        LinearToQuadraticMeshConverter<2> converter("mesh/test/data/", "square_128_elements", "TestLinearToQuadraticMeshConverter");

        // read in the new mesh to check it worked
        OutputFileHandler handler("TestLinearToQuadraticMeshConverter",false);
        std::string full_new_mesh = handler.GetOutputDirectoryFullPath() + "square_128_elements_quadratic";        

        QuadraticMesh<2> quad_mesh;
        TrianglesMeshReader<2,2> reader(full_new_mesh, 2, 2);
        quad_mesh.ConstructFromMeshReader(reader);
        
        TS_ASSERT_EQUALS(quad_mesh.GetNumNodes(), 17*17u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumElements(), 128u);
    }

    void TestLinearToQuadraticMeshConverter3d() throw(Exception)
    {
        LinearToQuadraticMeshConverter<3> converter("mesh/test/data/", "cube_136_elements", "TestLinearToQuadraticMeshConverter");

        // read in the new mesh to check it worked
        OutputFileHandler handler("TestLinearToQuadraticMeshConverter",false);
        std::string full_new_mesh = handler.GetOutputDirectoryFullPath() + "cube_136_elements_quadratic";        

        QuadraticMesh<3> quad_mesh;
        TrianglesMeshReader<3,3> reader(full_new_mesh, 2, 2);
        quad_mesh.ConstructFromMeshReader(reader);

        TS_ASSERT_EQUALS(quad_mesh.GetNumNodes(), 285u);
        TS_ASSERT_EQUALS(quad_mesh.GetNumElements(), 136u);
    }
};


#endif /*TESTLINEARTOQUADRATICMESHCONVERTER_HPP_*/
