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
        TrianglesMeshReader<2,2> *p_mesh_reader;
        p_mesh_reader=new READER_2D("mesh/test/data/disk_522_elements");
        delete p_mesh_reader;
    }

    /**
     * Check that the nodes are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing nodes) then an exception is thrown.
     */
    void TestNodesDataRead(void) throw(Exception)
    {
        TrianglesMeshReader<2,2> *p_mesh_reader;
        p_mesh_reader=new TrianglesMeshReader<2,2>(
                        "mesh/test/data/disk_984_elements_indexed_from_1");

        TS_ASSERT_EQUALS( p_mesh_reader->GetNumNodes(), 543U);
        delete p_mesh_reader;
        
        READER_2D mesh_reader("mesh/test/data/baddata/bad_nodes_disk_522_elements");

        // Reads node 0 from file
        TS_ASSERT_THROWS_NOTHING(mesh_reader.GetNextNode());
        // Reads node 3 from file when expecting number 1                            
        TS_ASSERT_THROWS_ANYTHING(mesh_reader.GetNextNode());                            
    }


    /**
     * Check that the elements are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing elements) then an exception is thrown.
     */
    void TestElementsDataRead(void) throw(Exception)
    {
        TrianglesMeshReader<2,2> *p_mesh_reader;
        p_mesh_reader=new TrianglesMeshReader<2,2>(
                        "mesh/test/data/disk_984_elements_indexed_from_1");

        TS_ASSERT_EQUALS( p_mesh_reader->GetNumElements(), 984U);
        delete p_mesh_reader;

        READER_2D mesh_reader("mesh/test/data/baddata/bad_elements_disk_522_elements");

        // Reads node 0 from file
        TS_ASSERT_THROWS_NOTHING(mesh_reader.GetNextElement());
        // Reads node 2 from file when expecting number 1                            
        TS_ASSERT_THROWS_ANYTHING(mesh_reader.GetNextElement());                            
    }


    /**
     * Check that the faces are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing faces) then an exception is thrown.
     */
    void TestFacesDataRead(void) throw(Exception)
    {
        TrianglesMeshReader<2,2> *p_mesh_reader;
        p_mesh_reader=new TrianglesMeshReader<2,2>(
                        "mesh/test/data/disk_984_elements_indexed_from_1");

        // TS_ASSERT_EQUALS( p_mesh_reader->GetNumFaces(), 1526U); // when all faces were read
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumFaces(), 100U); // just boundary faces are read

        delete p_mesh_reader;

        //TS_ASSERT_THROWS_ANYTHING(
        p_mesh_reader=new READER_2D(
                            "mesh/test/data/baddata/bad_faces_disk_522_elements"/*)*/);
        delete p_mesh_reader;//\todo This is temporary since presumably the constructor ought to throw?
    }



    /**
     * Checks that the reader can deal with (3-d) TetGen input files as well
     * as the previously considered (2-d) Triangles files. Checks that the
     * element output vector for a given input file is the correct length.
     */
    void Test3dDataRead(void) throw(Exception)
    {
        TrianglesMeshReader<3,3> *p_mesh_reader;
        p_mesh_reader=new READER_3D("mesh/test/data/slab_138_elements");

        TS_ASSERT_EQUALS (p_mesh_reader->GetNumElements(), 138U);
        delete p_mesh_reader;
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
        TrianglesMeshReader<2,2> *p_mesh_reader;
        p_mesh_reader=new TrianglesMeshReader<2,2>("mesh/test/data/disk_984_elements");

        std::vector<double> FirstNode;

        FirstNode = p_mesh_reader->GetNextNode();

        TS_ASSERT_DELTA( FirstNode[0] ,  0.9980267283 , 1e-6 );
        TS_ASSERT_DELTA( FirstNode[1] , -0.0627905195 , 1e-6 )

        std::vector<double> NextNode;

        NextNode = p_mesh_reader->GetNextNode();

        TS_ASSERT_DELTA( NextNode[0] , 1.0 , 1e-6 );
        TS_ASSERT_DELTA( NextNode[1] , 0.0 , 1e-6 );

        for (int i = 0; i < 541; i++)
        {
            TS_ASSERT_THROWS_NOTHING(NextNode = p_mesh_reader->GetNextNode());
        }

        TS_ASSERT_THROWS_ANYTHING(NextNode = p_mesh_reader->GetNextNode());
        delete p_mesh_reader;
    }



    /**
     * Check that GetNextElement() works. Checks that no errors are thrown for
     * all of the elements and that an error is thrown if we try to call the
     * function too many times.
     */
    void TestGetNextElement(void) throw(Exception)
    {
        TrianglesMeshReader<2,2> *p_mesh_reader;
        p_mesh_reader=new TrianglesMeshReader<2,2>("mesh/test/data/disk_984_elements");

        std::vector<unsigned> NextElement;

        for (unsigned i = 0; i < p_mesh_reader->GetNumElements(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(NextElement = p_mesh_reader->GetNextElement());
        }

        TS_ASSERT_THROWS_ANYTHING(NextElement = p_mesh_reader->GetNextElement());
        delete p_mesh_reader;
    }


    /**
     * Check that GetNextEdge() works. Checks that no errors are thrown for
     * all of the edges and that an error is thrown if we try to call the
     * function too many times.
     */
    void TestGetNextEdge(void) throw(Exception)
    {
        TrianglesMeshReader<2,2> *p_mesh_reader;
        p_mesh_reader=new TrianglesMeshReader<2,2>("mesh/test/data/disk_984_elements");

        std::vector<unsigned> NextEdge;

        TS_ASSERT_THROWS_NOTHING(NextEdge = p_mesh_reader->GetNextFace());
        TS_ASSERT_THROWS_NOTHING(NextEdge = p_mesh_reader->GetNextFace());

        for (unsigned i = 2; i < p_mesh_reader->GetNumEdges(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(NextEdge = p_mesh_reader->GetNextEdge());
        }

        TS_ASSERT_THROWS_ANYTHING(NextEdge = p_mesh_reader->GetNextEdge());
        delete p_mesh_reader;
    }


    /**
     * Check that the 1D data are read correctly. Check that the output vector
     * for a given input file is the correct length.
     */
    void Test1DMeshRead(void) throw(Exception)
    {
        TrianglesMeshReader<1,1> *p_mesh_reader;
        p_mesh_reader=new TrianglesMeshReader<1,1>("mesh/test/data/trivial_1d_mesh");

        TS_ASSERT_EQUALS( p_mesh_reader->GetNumNodes(), 11U);
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumElements(), 10U);
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumFaces(), 11U);

        // Determining boundary faces is no longer done by the MeshReader
        //TS_ASSERT_EQUALS( p_mesh_reader->GetNumBoundaryFaces(), 2);

        delete p_mesh_reader;
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
        TrianglesMeshReader<1,2> *p_mesh_reader;

        p_mesh_reader=new TrianglesMeshReader<1,2>("mesh/test/data/circle_outline");
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumNodes(), 100U);
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumElements(), 100U);
        delete p_mesh_reader;
        //Note: don't test faces (end nodes), since they are culled later
        p_mesh_reader=new TrianglesMeshReader<1,2>("mesh/test/data/semicircle_outline");
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumNodes(), 51U);
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumElements(), 50U);
        delete p_mesh_reader;
    }


    void Test2DMeshIn3DSpace() throw(Exception)
    {
        TrianglesMeshReader<2,3> *p_mesh_reader;

        p_mesh_reader=new TrianglesMeshReader<2,3>("mesh/test/data/slab_395_elements");

        TS_ASSERT_EQUALS( p_mesh_reader->GetNumNodes(), 132U);
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumElements(), 224U);
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumFaces(), 0U);

        delete p_mesh_reader;

        p_mesh_reader=new TrianglesMeshReader<2,3>("mesh/test/data/disk_in_3d");

        TS_ASSERT_EQUALS( p_mesh_reader->GetNumNodes(), 312U);
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumElements(), 522U);
        // TS_ASSERT_EQUALS( p_mesh_reader->GetNumFaces(), 833U); // when all faces were read
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumFaces(), 100U); // just boundary faces are read

        delete p_mesh_reader;
    }


    void TestOtherExceptions() throw(Exception)
    {
        // these should fail because SPACE_DIM doesn't match the dimension in the file
        //TS_ASSERT_THROWS_ANYTHING( READER_3D mesh_reader("mesh/test/data/disk_984_elements") );   
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

		std::vector<unsigned> next_element = mesh_reader.GetNextElement();

        TS_ASSERT_EQUALS(next_element.size(), 3u);

        TS_ASSERT_EQUALS(next_element[0], 0u);   // left node
        TS_ASSERT_EQUALS(next_element[1], 1u);   // right node
        TS_ASSERT_EQUALS(next_element[2], 11u);  // middle node
    
        for(unsigned i=1; i<10; i++)
        {
        	next_element = mesh_reader.GetNextElement();
    		TS_ASSERT_EQUALS(next_element.size(), 3u);
        }
        
    }    
    
    void TestReadingQuadraticMesh2d() throw(Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic", 2);
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 17u*17u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 128u);

		std::vector<unsigned> next_element = mesh_reader.GetNextElement();

        TS_ASSERT_EQUALS(next_element.size(), 6u);

        TS_ASSERT_EQUALS(next_element[0], 53u);
        TS_ASSERT_EQUALS(next_element[1], 0u);
        TS_ASSERT_EQUALS(next_element[2], 54u);
        TS_ASSERT_EQUALS(next_element[3], 82u); // opposite to 53
        TS_ASSERT_EQUALS(next_element[4], 83u); // opposite to 0
        TS_ASSERT_EQUALS(next_element[5], 81u); // opposite to 54
    
        for(unsigned i=1; i<128; i++)
        {
        	next_element = mesh_reader.GetNextElement();
    		TS_ASSERT_EQUALS(next_element.size(), 6u);
        }
    }
    
    void TestReadingQuadraticMesh3d() throw(Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element_quadratic", 2);
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 1u);

		std::vector<unsigned> next_element = mesh_reader.GetNextElement();
    
        TS_ASSERT_EQUALS(next_element.size(), 10u);
        
        for(unsigned i=0; i<10; i++)
        {
            TS_ASSERT_EQUALS(next_element[i], i);
        }
    }
};

#endif //_TESTTRIANGLESMESHREADER_HPP_
