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

#ifndef TESTINVENTORVORONOIWRITER_HPP_
#define TESTINVENTORVORONOIWRITER_HPP_


#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "VoronoiCell.hpp"
#include "VoronoiTessellation.cpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Exception.hpp"
#include "InventorVoronoiWriter.hpp"
#include "TrianglesMeshReader.cpp"

#include <cmath>
#include <vector>

class TestInventorVoronoiWriter : public CxxTest::TestSuite
{
public:
    void TestWriteTet() throw (Exception)
    {
        
        // Create conforming tetrahedral mesh which is Delaunay
        std::vector<Node<3> *> nodes;
        
        nodes.push_back(new Node<3>(0, true,  0.0,  0.0,  0.0));
        nodes.push_back(new Node<3>(1, true,  1.0,  1.0,  0.0));
        nodes.push_back(new Node<3>(2, true,  1.0,  0.0,  1.0));
        nodes.push_back(new Node<3>(3, true,  0.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(4, false, 0.5,  0.5,  0.5));
        
        ConformingTetrahedralMesh<3,3> mesh(nodes);
        
        // Create Voronoi Tesselation
        VoronoiTessellation<3> tessellation(mesh);
        
        InventorVoronoiWriter inventor_writer("InventorWriter", "SimpleTet");
        inventor_writer.Write(tessellation);
        
        // then compare against known good file
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "InventorWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/SimpleTet.iv mesh/test/data/InventorWriter/SimpleTet.iv").c_str()), 0); 
    }
    
    void TestWriteComplexCube() throw (Exception)
    {
        // Create conforming tetrahedral mesh which is Delaunay
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_136_elements");
        
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);
        
        TS_ASSERT(mesh.CheckVoronoi());
        
        // Create Voronoi Tesselation
        VoronoiTessellation<3> tessellation(mesh);
        
        InventorVoronoiWriter inventor_writer("InventorWriter", "Complex", false);
        inventor_writer.Write(tessellation);
        
        // then compare against known good file
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "InventorWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/Complex.iv mesh/test/data/InventorWriter/Complex.iv").c_str()), 0);        
    }
    
    void TestScaleAndWriteTet() throw (Exception)
    {
        // Create conforming tetrahedral mesh which is Delauny
        std::vector<Node<3> *> nodes;
        
        nodes.push_back(new Node<3>(0, true,  0.0,  0.0,  0.0));
        nodes.push_back(new Node<3>(1, true,  1.0,  1.0,  0.0));
        nodes.push_back(new Node<3>(2, true,  1.0,  0.0,  1.0));
        nodes.push_back(new Node<3>(3, true,  0.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(4, false, 0.5,  0.5,  0.5));
        
        ConformingTetrahedralMesh<3,3> mesh(nodes);

        // Create Voronoi Tesselation
        VoronoiTessellation<3> tessellation(mesh);
        
        InventorVoronoiWriter inventor_writer("InventorWriter", "ScaledSimpleTet");
        
        // fail because scale factor is not between 0 and 1
        TS_ASSERT_THROWS_ANYTHING(inventor_writer.ScaleAndWrite(tessellation,-1.0));
        TS_ASSERT_THROWS_ANYTHING(inventor_writer.ScaleAndWrite(tessellation, 2.0));
 
        inventor_writer.ScaleAndWrite(tessellation,2.0/3.0);

        // then compare against known good file
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "InventorWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/ScaledSimpleTet.iv mesh/test/data/InventorWriter/ScaledSimpleTet.iv").c_str()), 0);
    }

    // note that this is a different mesh to the non-scaled 'Complex' test
    void TestScaleAndWriteComplex() throw (Exception)
    {
        // Create conforming tetrahedral mesh which is Delaunay
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_152_elements");
        
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        // Create Voronoi Tesselation
        VoronoiTessellation<3> tessellation(mesh);
        
        InventorVoronoiWriter inventor_writer("InventorWriter", "ScaledComplex");
        inventor_writer.ScaleAndWrite(tessellation,0.5);
        
        // then compare against known good file
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "InventorWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/ScaledComplex.iv mesh/test/data/InventorWriter/ScaledComplex.iv").c_str()), 0);
    }

};
#endif /*TESTINVENTORVORONOIWRITER_HPP_*/
