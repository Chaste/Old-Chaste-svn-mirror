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
#ifndef TESTVERTEXMESHREADER2D_HPP_
#define TESTVERTEXMESHREADER2D_HPP_

#include <cxxtest/TestSuite.h>

#include "VertexMeshReader2d.hpp"

class TestVertexMeshReader2d : public CxxTest::TestSuite
{
public:

    /**
     * Check that input files are opened correctly.
     */
    void TestFilesOpen(void) throw(Exception)
    {
        VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/TestVertexMesh/vertex_mesh");
    }
    
    
    /**
     * Check that the nodes are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing nodes) then an exception is thrown.
     */
    void TestNodesDataRead(void) throw(Exception)
    {
        VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/TestVertexMesh/vertex_mesh");

        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 7u);

        VertexMeshReader2d mesh_reader2("notforrelease-cancer/test/data/baddata/vertex_mesh_bad_nodes");

        // Reads node 0 from file
        TS_ASSERT_THROWS_NOTHING(mesh_reader2.GetNextNode());

        // Reads node 3 from file when expecting number 1
        TS_ASSERT_THROWS_ANYTHING(mesh_reader2.GetNextNode());
    }


    /**
     * Check that the elements are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing elements) then an exception is thrown.
     */
    void TestElementsDataRead(void) throw(Exception)
    {
        VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/TestVertexMesh/vertex_mesh");
        
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 2u);
        
        // Read element 0 from file
        VertexElementData data = mesh_reader.GetNextElementData();

        TS_ASSERT_EQUALS(data.NodeIndices.size(), 5u);
        TS_ASSERT_EQUALS(data.NodeIndices[0], 0u);
        TS_ASSERT_EQUALS(data.NodeIndices[1], 1u);
        TS_ASSERT_EQUALS(data.NodeIndices[2], 2u);
        TS_ASSERT_EQUALS(data.NodeIndices[3], 3u);
        TS_ASSERT_EQUALS(data.NodeIndices[4], 4u);
        
        // Read element 1 from file
        VertexElementData data2 = mesh_reader.GetNextElementData();

        TS_ASSERT_EQUALS(data2.NodeIndices.size(), 3u);
        TS_ASSERT_EQUALS(data2.NodeIndices[0], 2u);
        TS_ASSERT_EQUALS(data2.NodeIndices[1], 5u);
        TS_ASSERT_EQUALS(data2.NodeIndices[2], 6u);
        
        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 0u);
        
        mesh_reader.Reset();
        for (unsigned i=1; i<mesh_reader.GetNumElements(); i++)
        {
            VertexElementData data = mesh_reader.GetNextElementData();
            TS_ASSERT_EQUALS(data.AttributeValue, 0u);
        }

        VertexMeshReader2d mesh_reader2("notforrelease-cancer/test/data/baddata/vertex_mesh_bad_elements");

        // Reads element 0 from file
        TS_ASSERT_THROWS_NOTHING(mesh_reader2.GetNextElementData());

        // Reads element 2 from file when expecting number 1                            
        TS_ASSERT_THROWS_ANYTHING(mesh_reader2.GetNextElementData());

    }
    

    /**
     * Checks that nodes in the input data file are numbered sequentially.
     * (In the input file nodes must appear in increasing order since the node
     * number is only stored as the index of the vector in which the coordinates
     * are stored.)
     */
    void TestPermutedNodesFail(void) throw(Exception)
    {
        VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/baddata/vertex_mesh_permuted_nodes");
        TS_ASSERT_THROWS_ANYTHING(for(unsigned i=0;i<mesh_reader.GetNumNodes();i++){mesh_reader.GetNextNode();})
    }


    /**
     * Check that GetNextNode() returns the coordinates of the correct node.
     * Compares the coordinates of the first two nodes with their known
     * values, checks that no errors are thrown for the remaining nodes and
     * that an error is thrown if we try to call the function too many times.
     */
    void TestGetNextNode(void) throw(Exception)
    {
        VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/TestVertexMesh/vertex_mesh");

        std::vector<double> first_node;
        first_node = mesh_reader.GetNextNode();

        TS_ASSERT_DELTA(first_node[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[1], 0.0, 1e-6)

        std::vector<double> next_node;
        next_node = mesh_reader.GetNextNode();

        TS_ASSERT_DELTA(next_node[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(next_node[1], 0.0, 1e-6);

        for (int i=0; i<5; i++)
        {
            TS_ASSERT_THROWS_NOTHING(next_node = mesh_reader.GetNextNode());
        }

        TS_ASSERT_THROWS_ANYTHING(next_node = mesh_reader.GetNextNode());
    }


    /**
     * Check that GetNextElementData() works. Checks that no errors are thrown for
     * all of the elements and that an error is thrown if we try to call the
     * function too many times.
     */
    void TestGetNextElementData(void) throw(Exception)
    {
        VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/TestVertexMesh/vertex_mesh");

        std::vector<unsigned> next_element;
        for (unsigned i=0; i<mesh_reader.GetNumElements(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(next_element = mesh_reader.GetNextElementData().NodeIndices);
        }

        TS_ASSERT_THROWS_ANYTHING(next_element = mesh_reader.GetNextElementData().NodeIndices);
    }


    void TestReadingElementAttributes() throw(Exception)
    {
        VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/TestVertexMeshReader2d/vertex_mesh_with_element_attributes");

        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 2u);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

        VertexElementData next_element_info = mesh_reader.GetNextElementData();
        std::vector<unsigned> nodes = next_element_info.NodeIndices;
        TS_ASSERT_EQUALS(nodes.size(), 5u);
        TS_ASSERT_EQUALS(next_element_info.AttributeValue, 97u);

        next_element_info = mesh_reader.GetNextElementData();
        nodes = next_element_info.NodeIndices;
        TS_ASSERT_EQUALS(nodes.size(), 3u);
        TS_ASSERT_EQUALS(next_element_info.AttributeValue, 152u)
    }

    void TestOtherExceptions() throw(Exception)
    {
        TS_ASSERT_THROWS_ANYTHING(VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/nonexistent_file"));
        TS_ASSERT_THROWS_ANYTHING(VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/baddata/vertex_mesh_without_element_file"));
    }

};
#endif /*TESTVERTEXMESHREADER2D_HPP_*/
