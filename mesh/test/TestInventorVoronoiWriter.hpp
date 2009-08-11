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


#ifndef TESTINVENTORVORONOIWRITER_HPP_
#define TESTINVENTORVORONOIWRITER_HPP_


#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "VoronoiCell.hpp"
#include "VoronoiTessellation.hpp"
#include "MutableMesh.hpp"
#include "Exception.hpp"
#include "InventorVoronoiWriter.hpp"
#include "TrianglesMeshReader.hpp"

#include <cmath>
#include <vector>

class TestInventorVoronoiWriter : public CxxTest::TestSuite
{
public:

    void TestWriteTet() throw (Exception)
    {

        // Create mutable tetrahedral mesh which is Delaunay
        std::vector<Node<3> *> nodes;

        nodes.push_back(new Node<3>(0, true,  0.0,  0.0,  0.0));
        nodes.push_back(new Node<3>(1, true,  1.0,  1.0,  0.0));
        nodes.push_back(new Node<3>(2, true,  1.0,  0.0,  1.0));
        nodes.push_back(new Node<3>(3, true,  0.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(4, false, 0.5,  0.5,  0.5));

        MutableMesh<3,3> mesh(nodes);
        Element<3,3> *p_element = mesh.GetElement(0);

        //TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 0u);//Older tetgen
        //TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 3u);//Older tetgen
        //TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2), 2u);//Older tetgen
        //TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3), 4u);//Older tetgen
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2), 0u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3), 2u);

        // Create Voronoi tessellation
        VoronoiTessellation<3> tessellation(mesh);

        InventorVoronoiWriter inventor_writer("InventorWriter", "SimpleTet");
        inventor_writer.Write(tessellation);

        // Then compare against known good file
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "InventorWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/SimpleTet.iv mesh/test/data/InventorWriter/SimpleTet.iv").c_str()), 0);
    }

    void TestWriteComplexCube() throw (Exception)
    {
        // Create mutable tetrahedral mesh which is Delaunay
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_136_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        TS_ASSERT(mesh.CheckVoronoi());

        // Create Voronoi Tesselation
        VoronoiTessellation<3> tessellation(mesh);

        InventorVoronoiWriter inventor_writer("InventorWriter", "Complex", false);
        inventor_writer.Write(tessellation);

        // Then compare against known good file
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "InventorWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/Complex.iv mesh/test/data/InventorWriter/Complex.iv").c_str()), 0);
    }

    void TestScaleAndWriteTet() throw (Exception)
    {
        // Create a conforming tetrahedral mesh which is Delaunay
        std::vector<Node<3> *> nodes;

        nodes.push_back(new Node<3>(0, true,  0.0,  0.0,  0.0));
        nodes.push_back(new Node<3>(1, true,  1.0,  1.0,  0.0));
        nodes.push_back(new Node<3>(2, true,  1.0,  0.0,  1.0));
        nodes.push_back(new Node<3>(3, true,  0.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(4, false, 0.5,  0.5,  0.5));

        MutableMesh<3,3> mesh(nodes);

        // Create Voronoi tessellation
        VoronoiTessellation<3> tessellation(mesh);

        InventorVoronoiWriter inventor_writer("InventorWriter", "ScaledSimpleTet");

        // These call throw exceptions because the scale factor is not between 0 and 1
        TS_ASSERT_THROWS_THIS(inventor_writer.ScaleAndWrite(tessellation, -1.0),"scaleFactor should be between 0 and 1");
        TS_ASSERT_THROWS_THIS(inventor_writer.ScaleAndWrite(tessellation, 2.0),"scaleFactor should be between 0 and 1");

        inventor_writer.ScaleAndWrite(tessellation, 2.0/3.0);

        // Then compare against known good file
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "InventorWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/ScaledSimpleTet.iv mesh/test/data/InventorWriter/ScaledSimpleTet.iv").c_str()), 0);
    }

    // Note that this is a different mesh to the non-scaled 'Complex' test
    void TestScaleAndWriteComplex() throw (Exception)
    {
        // Create mutable tetrahedral mesh which is Delaunay
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_152_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        // Create Voronoi tessellation
        VoronoiTessellation<3> tessellation(mesh);

        InventorVoronoiWriter inventor_writer("InventorWriter", "ScaledComplex");
        inventor_writer.ScaleAndWrite(tessellation, 0.5);

        // Then compare against a known good file
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "InventorWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/ScaledComplex.iv mesh/test/data/InventorWriter/ScaledComplex.iv").c_str()), 0);
    }

};
#endif /*TESTINVENTORVORONOIWRITER_HPP_*/
