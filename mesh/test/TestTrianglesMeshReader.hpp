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


/**
 * Test suite for the TrianglesMeshReader class.
 *
 */

#ifndef _TESTTRIANGLESMESHREADER_HPP_
#define _TESTTRIANGLESMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.hpp"

// these typedefs are just because can't have lines such as
//  TS_ASSERT_THROWS_NOTHING(p_mesh_reader=new TrianglesMeshReader<2,2>(name));
// because the macro thinks the comma separates two arguments
typedef TrianglesMeshReader<1,1> READER_1D;
typedef TrianglesMeshReader<2,2> READER_2D;
typedef TrianglesMeshReader<3,3> READER_3D;
typedef TrianglesMeshReader<0,1> READER_0D_IN_1D;


class TestTrianglesMeshReader : public CxxTest::TestSuite
{
public:

    /**
     * Check that input files are opened correctly.
     */
    void TestFilesOpen(void) throw(Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
    }

    /**
     * Check that the nodes are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing nodes) then an exception is thrown.
     */
    void TestNodesDataRead(void) throw(Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements_indexed_from_1");

        TS_ASSERT_EQUALS( mesh_reader.GetNumNodes(), 543U);
        
        TrianglesMeshReader<2,2> mesh_reader2("mesh/test/data/baddata/bad_nodes_disk_522_elements");

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
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements_indexed_from_1");

        TS_ASSERT_EQUALS( mesh_reader.GetNumElements(), 984U);

        ElementData data = mesh_reader.GetNextElementData();
        TS_ASSERT_EQUALS(data.NodeIndices.size(), 3u);
        TS_ASSERT_EQUALS(data.NodeIndices[0], 309u);
        TS_ASSERT_EQUALS(data.NodeIndices[1], 144u);
        TS_ASSERT_EQUALS(data.NodeIndices[2], 310u);
        
        TS_ASSERT_EQUALS( mesh_reader.GetNumElementAttributes(), 0u);
        for(unsigned i=1; i<mesh_reader.GetNumElements(); i++)
        {
            ElementData data = mesh_reader.GetNextElementData();
            TS_ASSERT_EQUALS(data.AttributeValue, 0u);
        }

        TrianglesMeshReader<2,2> mesh_reader2("mesh/test/data/baddata/bad_elements_disk_522_elements");

        // Reads element 0 from file
        TS_ASSERT_THROWS_NOTHING(mesh_reader2.GetNextElementData());
        // Reads element 2 from file when expecting number 1                            
        TS_ASSERT_THROWS_ANYTHING(mesh_reader2.GetNextElementData());
    }


    /**
     * Check that the faces are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing faces) then an exception is thrown.
     */
    void TestFacesDataRead(void) throw(Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements_indexed_from_1");

        // TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 1526U); // when all faces were read
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 100U); // just boundary faces are read

        TrianglesMeshReader<2,2> mesh_reader2("mesh/test/data/baddata/bad_faces_disk_522_elements");

        // First boundary face is #20, on its way through the file there's a gap between face 1 and face 10
        TS_ASSERT_THROWS_ANYTHING(mesh_reader2.GetNextFace());
        /// \todo: the exception thrown here is thrown in the constructor as well. In that case, it's wrongly
        /// caught considering it an end of file exception and therefore setting an inconsistent number of faces.
    }



    /**
     * Checks that the reader can deal with (3-d) TetGen input files as well
     * as the previously considered (2-d) Triangles files. Checks that the
     * element output vector for a given input file is the correct length.
     */
    void Test3dDataRead(void) throw(Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/slab_138_elements");

        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 138U);
    }


    /**
     * Checks that nodes in the input data file are numbered sequentially.
     * (In the input file nodes must appear in increasing order since the node
     * number is only stored as the index of the vector in which the coordinates
     * are stored.)
     */
    void TestPermutedNodesFail(void) throw(Exception)
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/baddata/permuted_nodes_disk_522_elements");
        TS_ASSERT_THROWS_ANYTHING(for(unsigned i=0;i<reader.GetNumNodes();i++){reader.GetNextNode();})
    }


    /**
     * Checks that elements have the correct number of nodes (i.e. one more
     * node than the dimension of the mesh). If quadratic basis functions are
     * required this should be dealt with elsewhere.
     */
    void TestOrder2ElementsFail(void) throw(Exception)
    {
        TrianglesMeshReader<2,2> *p_mesh_reader;
        TS_ASSERT_THROWS_ANYTHING(
            p_mesh_reader=new READER_2D(
                            "mesh/test/data/baddata/disk_522_order_2_elements"));
    }


    /**
     * Check that GetNextNode() returns the coordinates of the correct node.
     * Compares the coordinates of the first two nodes with their known
     * values, checks that no errors are thrown for the remaining nodes and
     * that an error is thrown if we try to call the function too many times.
     */
    void TestGetNextNode(void) throw(Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");

        std::vector<double> first_node;

        first_node = mesh_reader.GetNextNode();

        TS_ASSERT_DELTA( first_node[0] ,  0.9980267283 , 1e-6 );
        TS_ASSERT_DELTA( first_node[1] , -0.0627905195 , 1e-6 )

        std::vector<double> next_node;

        next_node = mesh_reader.GetNextNode();

        TS_ASSERT_DELTA( next_node[0] , 1.0 , 1e-6 );
        TS_ASSERT_DELTA( next_node[1] , 0.0 , 1e-6 );

        for (int i = 0; i < 541; i++)
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
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");

        std::vector<unsigned> next_element;

        for (unsigned i=0; i<mesh_reader.GetNumElements(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(next_element = mesh_reader.GetNextElementData().NodeIndices);
        }

        TS_ASSERT_THROWS_ANYTHING(next_element = mesh_reader.GetNextElementData().NodeIndices);
    }


    /**
     * Check that GetNextEdge() works. Checks that no errors are thrown for
     * all of the edges and that an error is thrown if we try to call the
     * function too many times.
     */
    void TestGetNextEdge(void) throw(Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");

        std::vector<unsigned> next_edge;

        TS_ASSERT_THROWS_NOTHING(next_edge = mesh_reader.GetNextFace());
        TS_ASSERT_THROWS_NOTHING(next_edge = mesh_reader.GetNextFace());

        for (unsigned i = 2; i < mesh_reader.GetNumEdges(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(next_edge = mesh_reader.GetNextEdge());
        }

        TS_ASSERT_THROWS_ANYTHING(next_edge = mesh_reader.GetNextEdge());
    }


    /**
     * Check that the 1D data are read correctly. Check that the output vector
     * for a given input file is the correct length.
     */
    void Test1DMeshRead(void) throw(Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/trivial_1d_mesh");

        TS_ASSERT_EQUALS( mesh_reader.GetNumNodes(), 11U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumElements(), 10U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 11U);

        // Determining boundary faces is no longer done by the MeshReader
        //TS_ASSERT_EQUALS( mesh_reader.GetNumBoundaryFaces(), 2);
    }


    void Test0DMeshIn1DSpaceFails() throw(Exception)
    {
        TrianglesMeshReader<0,1> *p_mesh_reader;
        TS_ASSERT_THROWS_ANYTHING(  p_mesh_reader=new READER_0D_IN_1D(
                                                    "mesh/test/data/trivial_1d_mesh")
                                 );
    }


    void Test1DMeshIn2DSpace() throw(Exception)
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/circle_outline");
        TS_ASSERT_EQUALS( mesh_reader.GetNumNodes(), 100U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumElements(), 100U);

        //Note: don't test faces (end nodes), since they are culled later
        TrianglesMeshReader<1,2> mesh_reader2("mesh/test/data/semicircle_outline");
        TS_ASSERT_EQUALS( mesh_reader2.GetNumNodes(), 51U);
        TS_ASSERT_EQUALS( mesh_reader2.GetNumElements(), 50U);
    }


    void Test2DMeshIn3DSpace() throw(Exception)
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/slab_395_elements");

        TS_ASSERT_EQUALS( mesh_reader.GetNumNodes(), 132U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumElements(), 224U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 0U);


        TrianglesMeshReader<2,3> mesh_reader2("mesh/test/data/disk_in_3d");

        TS_ASSERT_EQUALS( mesh_reader2.GetNumNodes(), 312U);
        TS_ASSERT_EQUALS( mesh_reader2.GetNumElements(), 522U);
        // TS_ASSERT_EQUALS( mesh_reader2.GetNumFaces(), 833U); // when all faces were read
        TS_ASSERT_EQUALS( mesh_reader2.GetNumFaces(), 100U); // just boundary faces are read
    }


    void TestOtherExceptions() throw(Exception)
    {
        // these should fail because SPACE_DIM doesn't match the dimension in the file   
        TS_ASSERT_THROWS_ANYTHING( READER_1D mesh_reader("mesh/test/data/disk_984_elements") );
    }

    ////////////////////////////////////////////////////////
    // quadratic tests
    ////////////////////////////////////////////////////////
    void TestReadingQuadraticMesh1d() throw(Exception)
    {
        TS_ASSERT_THROWS_ANYTHING( READER_1D wrong_reader("mesh/test/data/1D_0_to_1_10_elements_quadratics", 1));

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_quadratic", 2);
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 21u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 10u);

		std::vector<unsigned> next_element = mesh_reader.GetNextElementData().NodeIndices;

        TS_ASSERT_EQUALS(next_element.size(), 3u);

        TS_ASSERT_EQUALS(next_element[0], 0u);   // left node
        TS_ASSERT_EQUALS(next_element[1], 1u);   // right node
        TS_ASSERT_EQUALS(next_element[2], 11u);  // middle node
    
        for(unsigned i=1; i<10; i++)
        {
        	next_element = mesh_reader.GetNextElementData().NodeIndices;
    		TS_ASSERT_EQUALS(next_element.size(), 3u);
        }
        
    }    
    
    void TestReadingQuadraticMesh2d() throw(Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic", 2);
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 17u*17u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 128u);

		std::vector<unsigned> next_element = mesh_reader.GetNextElementData().NodeIndices;

        TS_ASSERT_EQUALS(next_element.size(), 6u);

        TS_ASSERT_EQUALS(next_element[0], 53u);
        TS_ASSERT_EQUALS(next_element[1], 0u);
        TS_ASSERT_EQUALS(next_element[2], 54u);
        TS_ASSERT_EQUALS(next_element[3], 82u); // opposite to 53
        TS_ASSERT_EQUALS(next_element[4], 83u); // opposite to 0
        TS_ASSERT_EQUALS(next_element[5], 81u); // opposite to 54
    
        for(unsigned i=1; i<128; i++)
        {
        	next_element = mesh_reader.GetNextElementData().NodeIndices;
    		TS_ASSERT_EQUALS(next_element.size(), 6u);
        }
    }
    
    void TestReadingQuadraticMesh3d() throw(Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element_quadratic", 2);
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 1u);

		std::vector<unsigned> next_element = mesh_reader.GetNextElementData().NodeIndices;
    
        TS_ASSERT_EQUALS(next_element.size(), 10u);
        
        for(unsigned i=0; i<10; i++)
        {
            TS_ASSERT_EQUALS(next_element[i], i);
        }
    }
    void TestReadingElementAttributes() throw(Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_with_attributes");

        TS_ASSERT_EQUALS( mesh_reader.GetNumElements(), 10U);
        
        TS_ASSERT_EQUALS( mesh_reader.GetNumElementAttributes(), 1U);
        
        
        for(unsigned i=0; i<10; i++)
        {
            ElementData next_element_info = mesh_reader.GetNextElementData();
            std::vector<unsigned> nodes = next_element_info.NodeIndices;
            TS_ASSERT_EQUALS(nodes.size(), 2u);
            TS_ASSERT_EQUALS(next_element_info.AttributeValue, i%5+1);
        }        
    }
};

#endif //_TESTTRIANGLESMESHREADER_HPP_
