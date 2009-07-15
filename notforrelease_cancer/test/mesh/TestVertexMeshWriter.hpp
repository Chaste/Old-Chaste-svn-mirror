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
#ifndef TESTVERTEXMESHWRITER_HPP_
#define TESTVERTEXMESHWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <string>
#include <fstream>

#include "VertexMeshWriter.hpp"
#include "VtkWriter.hpp"

class TestVertexMeshWriter : public CxxTest::TestSuite
{
public:

    void TestMeshWriter() throw(Exception)
    {
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        basic_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        basic_nodes.push_back(new Node<2>(2, false, 1.5, 1.0));
        basic_nodes.push_back(new Node<2>(3, false, 1.0, 2.0));
        basic_nodes.push_back(new Node<2>(4, false, 0.0, 1.0));
        basic_nodes.push_back(new Node<2>(5, false, 2.0, 0.0));
        basic_nodes.push_back(new Node<2>(6, false, 2.0, 3.0));

        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;

        // Make two triangular elements out of these nodes
        nodes_elem_0.push_back(basic_nodes[0]);
        nodes_elem_0.push_back(basic_nodes[1]);
        nodes_elem_0.push_back(basic_nodes[2]);
        nodes_elem_0.push_back(basic_nodes[3]);
        nodes_elem_0.push_back(basic_nodes[4]);

        nodes_elem_1.push_back(basic_nodes[2]);
        nodes_elem_1.push_back(basic_nodes[5]);
        nodes_elem_1.push_back(basic_nodes[6]);

        std::vector<VertexElement<2,2>*> basic_vertex_elements;
        basic_vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        basic_vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        VertexMesh<2,2> basic_vertex_mesh(basic_nodes, basic_vertex_elements);

        // Create a vertex mesh writer
        VertexMeshWriter<2,2> vertex_mesh_writer("TestVertexMeshWriter", "vertex_mesh");
        vertex_mesh_writer.WriteFilesUsingMesh(basic_vertex_mesh);

        OutputFileHandler handler("TestVertexMeshWriter", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "vertex_mesh.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "vertex_mesh.cell";

        TS_ASSERT_EQUALS(system(("diff " + results_file1 + " notforrelease_cancer/test/data/TestVertexMesh/vertex_mesh.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_file2 + " notforrelease_cancer/test/data/TestVertexMesh/vertex_mesh.cell").c_str()), 0);
    }


    void TestVtkVertexMeshWriter() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        basic_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        basic_nodes.push_back(new Node<2>(2, false, 1.5, 1.0));
        basic_nodes.push_back(new Node<2>(3, false, 1.0, 2.0));
        basic_nodes.push_back(new Node<2>(4, false, 0.0, 1.0));
        basic_nodes.push_back(new Node<2>(5, false, 2.0, 0.0));
        basic_nodes.push_back(new Node<2>(6, false, 2.0, 3.0));

        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;

        // Make two triangular elements out of these nodes
        nodes_elem_0.push_back(basic_nodes[0]);
        nodes_elem_0.push_back(basic_nodes[1]);
        nodes_elem_0.push_back(basic_nodes[2]);
        nodes_elem_0.push_back(basic_nodes[3]);
        nodes_elem_0.push_back(basic_nodes[4]);

        nodes_elem_1.push_back(basic_nodes[2]);
        nodes_elem_1.push_back(basic_nodes[5]);
        nodes_elem_1.push_back(basic_nodes[6]);

        std::vector<VertexElement<2,2>*> basic_vertex_elements;
        basic_vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        basic_vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        VertexMesh<2,2> basic_vertex_mesh(basic_nodes, basic_vertex_elements);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0u)->IsDeleted(), false);
   
        //AbstractMesh<2,2> *p_abs_mesh=&basic_vertex_mesh;
        //TS_ASSERT_EQUALS(p_abs_mesh->GetNumElements(), 2u);
        
        //The next line segfaults since VertexMesh and AbstractTetrahedralMesh
        //both have a vector called mElements
        
        //***** TS_ASSERT_EQUALS(p_abs_mesh->GetElement(0u)->IsDeleted(), false);
   
        // Create a vertex mesh writer
        //VtkWriter<2> vertex_mesh_writer("TestVertexMeshWriter", "vertex_mesh");
        //vertex_mesh_writer.WriteFilesUsingMesh(basic_vertex_mesh);
#endif //CHASTE_VTK
    }
};


#endif /*TESTVERTEXMESHWRITER_HPP_*/
