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


// TestTrianglesMeshReader.hpp

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


class TestTrianglesMeshReaders : public CxxTest::TestSuite
{
public:

    /**
     * Check that input files are opened correctly.
     */
    void TestFilesOpen(void)
    {
        AbstractMeshReader<2,2> *p_mesh_reader;
        p_mesh_reader=new READER_2D("mesh/test/data/disk_522_elements");
        delete p_mesh_reader;
    }

    /**
     * Check that the nodes are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing nodes) then an exception is thrown.
     */
    void TestNodesDataRead(void)
    {
        AbstractMeshReader<2,2> *p_mesh_reader;
        p_mesh_reader=new TrianglesMeshReader<2,2>(
                        "mesh/test/data/disk_984_elements_indexed_from_1");

        TS_ASSERT_EQUALS( p_mesh_reader->GetNumNodes(), 543U);
        delete p_mesh_reader;

        TS_ASSERT_THROWS_ANYTHING(
            p_mesh_reader=new READER_2D(
                            "mesh/test/data/baddata/bad_nodes_disk_522__elements_indexed_from_1"));
    }


    /**
     * Check that the elements are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing elements) then an exception is thrown.
     */
    void TestElementsDataRead(void)
    {
        AbstractMeshReader<2,2> *p_mesh_reader;
        p_mesh_reader=new TrianglesMeshReader<2,2>(
                        "mesh/test/data/disk_984_elements_indexed_from_1");

        TS_ASSERT_EQUALS( p_mesh_reader->GetNumElements(), 984U);
        delete p_mesh_reader;

        TS_ASSERT_THROWS_ANYTHING(
            p_mesh_reader=new READER_2D(
                            "mesh/test/data/baddata/bad_elements_disk_522_elements_indexed_from_1"));
    }


    /**
     * Check that the faces are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing faces) then an exception is thrown.
     */
    void TestFacesDataRead(void)
    {
        AbstractMeshReader<2,2> *p_mesh_reader;
        p_mesh_reader=new TrianglesMeshReader<2,2>(
                        "mesh/test/data/disk_984_elements_indexed_from_1");

        // TS_ASSERT_EQUALS( p_mesh_reader->GetNumFaces(), 1526U); // when all faces were read
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumFaces(), 100U); // just boundary faces are read

        delete p_mesh_reader;

        TS_ASSERT_THROWS_ANYTHING(
            p_mesh_reader=new READER_2D(
                            "mesh/test/data/baddata/bad_faces_disk_522__elements_indexed_from_1"));
    }



    /**
     * Checks that the reader can deal with (3-d) TetGen input files as well
     * as the previously considered (2-d) Triangles files. Checks that the
     * element output vector for a given input file is the correct length.
     */
    void Test3dDataRead(void)
    {
        AbstractMeshReader<3,3> *p_mesh_reader;
        p_mesh_reader=new READER_3D("mesh/test/data/slab_138_elements");

        TS_ASSERT_EQUALS (p_mesh_reader->GetNumElements(), 138U);
        delete p_mesh_reader;
    }


    /**
     * Checks that nodes are indexed from zero. Takes input file that is
     * indexed from zero and checks that the output file also is. Uses methods
     * GetMaxNodeIndex() and GetMinNodeIndex().
     */
    void TestIndexFromZero(void)
    {
        AbstractMeshReader<2,2> *p_mesh_reader;
        p_mesh_reader=new READER_2D("mesh/test/data/disk_522_elements");

        TS_ASSERT_EQUALS(p_mesh_reader->GetMaxNodeIndex(), p_mesh_reader->GetNumNodes() - 1);

        TS_ASSERT_EQUALS(p_mesh_reader->GetMinNodeIndex(), 0U);
        delete p_mesh_reader;
    }


    /**
     * Checks that nodes are indexed from zero. Takes input file that is
     * indexed from one and checks that the output file also is. Uses methods
     * GetMaxNodeIndex() and GetMinNodeIndex().
     */
    void TestIndexFromOne(void)
    {
        AbstractMeshReader<2,2> *p_mesh_reader;
        p_mesh_reader=new READER_2D(
                            "mesh/test/data/disk_522_elements_indexed_from_1");

        TS_ASSERT_EQUALS(p_mesh_reader->GetMaxNodeIndex(), p_mesh_reader->GetNumNodes() - 1);

        TS_ASSERT_EQUALS(p_mesh_reader->GetMinNodeIndex(), 0U);
        delete p_mesh_reader;
    }



    /**
     * Checks that nodes in the input data file are numbered sequentially.
     * (In the input file nodes must appear in increasing order since the node
     * number is only stored as the index of the vector in which the coordinates
     * are stored.)
     */
    void TestPermutedNodesFail(void)
    {
        AbstractMeshReader<2,2> *p_mesh_reader;
        TS_ASSERT_THROWS_ANYTHING(
            p_mesh_reader=new READER_2D(
                            "mesh/test/data/baddata/permuted_nodes_disk_522_elements"));
    }


    /**
     * Checks that elements have the correct number of nodes (i.e. one more
     * node than the dimension of the mesh). If quadratic basis functions are
     * required this should be dealt with elsewhere.
     */
    void TestOrder2ElementsFail(void)
    {
        AbstractMeshReader<2,2> *p_mesh_reader;
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
    void TestGetNextNode(void)
    {
        AbstractMeshReader<2,2> *p_mesh_reader;
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
    void TestGetNextElement(void)
    {
        AbstractMeshReader<2,2> *p_mesh_reader;
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
    void TestGetNextEdge(void)
    {
        AbstractMeshReader<2,2> *p_mesh_reader;
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
    void Test1DMeshRead(void)
    {
        AbstractMeshReader<1,1> *p_mesh_reader;
        p_mesh_reader=new TrianglesMeshReader<1,1>("mesh/test/data/trivial_1d_mesh");

        TS_ASSERT_EQUALS( p_mesh_reader->GetNumNodes(), 11U);
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumElements(), 10U);
        TS_ASSERT_EQUALS( p_mesh_reader->GetNumFaces(), 11U);

        // Determining boundary faces is no longer done by the MeshReader
        //TS_ASSERT_EQUALS( p_mesh_reader->GetNumBoundaryFaces(), 2);

        delete p_mesh_reader;
    }


    void Test0DMeshIn1DSpaceFails()
    {
        AbstractMeshReader<0,1> *p_mesh_reader;
        TS_ASSERT_THROWS_ANYTHING(  p_mesh_reader=new READER_0D_IN_1D(
                                                    "mesh/test/data/trivial_1d_mesh")
                                 );
    }


    void Test1DMeshIn2DSpace()
    {
        AbstractMeshReader<1,2> *p_mesh_reader;

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


    void Test2DMeshIn3DSpace()
    {
        AbstractMeshReader<2,3> *p_mesh_reader;

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


    void TestExceptions()
    {
        // this should fail because a 3D reader is attempting to read a 2D mesh
        TS_ASSERT_THROWS_ANYTHING( READER_3D mesh_reader("mesh/test/data/disk_984_elements") );
    }

};

#endif //_TESTTRIANGLESMESHREADER_HPP_
