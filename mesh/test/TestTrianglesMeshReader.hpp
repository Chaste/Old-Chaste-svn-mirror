// TestTrianglesMeshReader.hpp

/**
 * Test suite for the TrianglesMeshReader class.
 *
 */

#ifndef _TESTTRIANGLESMESHREADER_HPP_
#define _TESTTRIANGLESMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.cpp"

// these typedefs are just because can't have lines such as
//  TS_ASSERT_THROWS_NOTHING(pMeshReader=new TrianglesMeshReader<2,2>(name));
// because the macro thinks the comma seperates two arguments
typedef TrianglesMeshReader<1,1> READER_1D;
typedef TrianglesMeshReader<2,2> READER_2D;
typedef TrianglesMeshReader<3,3> READER_3D;
typedef TrianglesMeshReader<0,1> READER_0D_IN_1D;


class TestTrianglesMeshReaders : public CxxTest::TestSuite
{
public:

    /**
     * Check that input files are opened correctly.
     * 
     */
    
    void testFilesOpen(void)
    {
        AbstractMeshReader<2,2> *pMeshReader;
        TS_ASSERT_THROWS_NOTHING(
            pMeshReader=new READER_2D(
                            "mesh/test/data/disk_522_elements"));
        delete pMeshReader;
        
    }
    
    
    
    
    /**
     * Check that the nodes are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing nodes) then an exception is thrown.
     * 
     */
    
    void TestNodesDataRead(void)
    {
        AbstractMeshReader<2,2> *pMeshReader;
        pMeshReader=new TrianglesMeshReader<2,2>(
                        "mesh/test/data/disk_984_elements_indexed_from_1");
                        
        TS_ASSERT_EQUALS( pMeshReader->GetNumNodes(), 543);
        delete pMeshReader;
        
        TS_ASSERT_THROWS_ANYTHING(
            pMeshReader=new READER_2D(
                            "mesh/test/data/baddata/bad_nodes_disk_522__elements_indexed_from_1"));
    }
    
    
    
    /**
     * Check that the elements are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing elements) then an exception is thrown.
     * 
     */
    
    void TestElementsDataRead(void)
    {
        AbstractMeshReader<2,2> *pMeshReader;
        pMeshReader=new TrianglesMeshReader<2,2>(
                        "mesh/test/data/disk_984_elements_indexed_from_1");
                        
        TS_ASSERT_EQUALS( pMeshReader->GetNumElements(), 984);
        delete pMeshReader;
        
        TS_ASSERT_THROWS_ANYTHING(
            pMeshReader=new READER_2D(
                            "mesh/test/data/baddata/bad_elements_disk_522_elements_indexed_from_1"));
                            
    }
    
    
    
    /**
     * Check that the faces are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing faces) then an exception is thrown.
     * 
     */
    
    void TestFacesDataRead(void)
    {
        AbstractMeshReader<2,2> *pMeshReader;
        pMeshReader=new TrianglesMeshReader<2,2>(
                        "mesh/test/data/disk_984_elements_indexed_from_1");
                        
        TS_ASSERT_EQUALS( pMeshReader->GetNumFaces(), 1526);
        delete pMeshReader;
        
        TS_ASSERT_THROWS_ANYTHING(
            pMeshReader=new READER_2D(
                            "mesh/test/data/baddata/bad_faces_disk_522__elements_indexed_from_1"));
    }
    
    
    
    /**
     * Checks that the reader can deal with (3-d) TetGen input files as well
     * as the previously considered (2-d) Triangles files. Checks that the
     * element output vector for a given input file is the correct length.
     * 
     */
    
    void Test3dDataRead(void)
    {
        AbstractMeshReader<3,3> *pMeshReader;
        TS_ASSERT_THROWS_NOTHING(
            pMeshReader=new READER_3D(
                            "mesh/test/data/slab_138_elements"));
                            
        TS_ASSERT_EQUALS (pMeshReader->GetNumElements(), 138);
        delete pMeshReader;
    }
    
    
    /**
     * Checks that nodes are indexed from zero. Takes input file that is
     * indexed from zero and checks that the output file also is. Uses methods
     * GetMaxNodeIndex() and GetMinNodeIndex().
     * 
     */
    
    void TestIndexFromZero(void)
    {
        AbstractMeshReader<2,2> *pMeshReader;
        TS_ASSERT_THROWS_NOTHING(
            pMeshReader=new READER_2D(
                            "mesh/test/data/disk_522_elements"));
                            
        TS_ASSERT_EQUALS(pMeshReader->GetMaxNodeIndex(), pMeshReader->GetNumNodes() - 1);
        
        TS_ASSERT_EQUALS(pMeshReader->GetMinNodeIndex(), 0);
        delete pMeshReader;
    }
    
    
    
    /**
     * Checks that nodes are indexed from zero. Takes input file that is
     * indexed from one and checks that the output file also is. Uses methods
     * GetMaxNodeIndex() and GetMinNodeIndex().
     * 
     */
    
    void TestIndexFromOne(void)
    {
        AbstractMeshReader<2,2> *pMeshReader;
        TS_ASSERT_THROWS_NOTHING(
            pMeshReader=new READER_2D(
                            "mesh/test/data/disk_522_elements_indexed_from_1"));
                            
        TS_ASSERT_EQUALS(pMeshReader->GetMaxNodeIndex(), pMeshReader->GetNumNodes() - 1);
        
        TS_ASSERT_EQUALS(pMeshReader->GetMinNodeIndex(), 0);
        delete pMeshReader;
    }
    
    
    
    /**
     * Checks that nodes in the input data file are numbered sequentially.
     * (In the input file nodes must appear in increasing order since the node
     * number is only stored as the index of the vector in which the coordinates
     * are stored.)
     * 
     */
    
    void TestPermutedNodesFail(void)
    {
        AbstractMeshReader<2,2> *pMeshReader;
        TS_ASSERT_THROWS_ANYTHING(
            pMeshReader=new READER_2D(
                            "mesh/test/data/baddata/permuted_nodes_disk_522_elements"));
                            
    }
    
    
    /**
     * Checks that elements have the correct number of nodes (i.e. one more
     * node than the dimension of the mesh). If quadratic basis functions are
     * required this should be dealt with elsewhere. 
     * 
     */
    
    void TestOrder2ElementsFail(void)
    {
        AbstractMeshReader<2,2> *pMeshReader;
        TS_ASSERT_THROWS_ANYTHING(
            pMeshReader=new READER_2D(
                            "mesh/test/data/baddata/disk_522_order_2_elements"));
    }
    
    
    /**
     * Check that GetNextNode() returns the coordinates of the correct node.
     * Compares the coordinates of the first two nodes with their known
     * values, checks that no errors are thrown for the remaining nodes and
     * that an error is thrown if we try to call the function too many times.
     * 
     */
    
    void TestGetNextNode(void)
    {
        AbstractMeshReader<2,2> *pMeshReader;
        pMeshReader=new TrianglesMeshReader<2,2>(
                        "mesh/test/data/disk_984_elements");
                        
        std::vector<double> FirstNode;
        
        FirstNode = pMeshReader->GetNextNode();
        
        TS_ASSERT_DELTA( FirstNode[0] ,  0.9980267283 , 1e-6 );
        TS_ASSERT_DELTA( FirstNode[1] , -0.0627905195 , 1e-6 )
        
        std::vector<double> NextNode;
        
        NextNode = pMeshReader->GetNextNode();
        
        TS_ASSERT_DELTA( NextNode[0] , 1.0 , 1e-6 );
        TS_ASSERT_DELTA( NextNode[1] , 0.0 , 1e-6 );
        
        for (int i = 0; i < 541; i++)
        {
            TS_ASSERT_THROWS_NOTHING(NextNode = pMeshReader->GetNextNode());
        }
        
        TS_ASSERT_THROWS_ANYTHING(NextNode = pMeshReader->GetNextNode());
        delete pMeshReader;
    }
    
    
    
    /**
     * Check that GetNextElement() works. Checks that no errors are thrown for 
     * all of the elements and that an error is thrown if we try to call the 
     * function too many times.
     * 
     */
    
    void TestGetNextElement(void)
    {
        AbstractMeshReader<2,2> *pMeshReader;
        pMeshReader=new TrianglesMeshReader<2,2>(
                        "mesh/test/data/disk_984_elements");
                        
        std::vector<int> NextElement;
        
        for (int i = 0; i < pMeshReader->GetNumElements(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(NextElement = pMeshReader->GetNextElement());
        }
        
        TS_ASSERT_THROWS_ANYTHING(NextElement = pMeshReader->GetNextElement());
        delete pMeshReader;
    }
    
    
    /**
     * Check that GetNextEdge() works. Checks that no errors are thrown for 
     * all of the edges and that an error is thrown if we try to call the 
     * function too many times.
     * 
     */
    
    void TestGetNextEdge(void)
    {
        AbstractMeshReader<2,2> *pMeshReader;
        pMeshReader=new TrianglesMeshReader<2,2>(
                        "mesh/test/data/disk_984_elements");
                        
        std::vector<int> NextEdge;
        
        TS_ASSERT_THROWS_NOTHING(NextEdge = pMeshReader->GetNextFace());
        TS_ASSERT_THROWS_NOTHING(NextEdge = pMeshReader->GetNextFace());
        
        for (int i = 2; i < pMeshReader->GetNumEdges(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(NextEdge = pMeshReader->GetNextEdge());
        }
        
        TS_ASSERT_THROWS_ANYTHING(NextEdge = pMeshReader->GetNextEdge());
        delete pMeshReader;
    }
    
    /**
     * Check that the 1D data are read correctly. Check that the output vector
     * for a given input file is the correct length.
     */
    
    void Test1DMeshRead(void)
    {
        AbstractMeshReader<1,1> *pMeshReader;
        pMeshReader=new TrianglesMeshReader<1,1>(
                        "mesh/test/data/trivial_1d_mesh");
                        
        TS_ASSERT_EQUALS( pMeshReader->GetNumNodes(), 11);
        
        TS_ASSERT_EQUALS( pMeshReader->GetNumElements(), 10);
        
        TS_ASSERT_EQUALS( pMeshReader->GetNumFaces(), 11);
        
        // Determining boundary faces is no longer done by the MeshReader
        //TS_ASSERT_EQUALS( pMeshReader->GetNumBoundaryFaces(), 2);
        
        delete pMeshReader;
        
    }
    
    
    void Test0DMeshIn1DSpaceFails()
    {
    
        AbstractMeshReader<0,1> *pMeshReader;
        
        TS_ASSERT_THROWS_ANYTHING(  pMeshReader=new READER_0D_IN_1D(
                                                    "mesh/test/data/trivial_1d_mesh")
                                 );
    }
    
    void Test1DMeshIn2DSpace()
    {
    
        AbstractMeshReader<1,2> *pMeshReader;
        
        pMeshReader=new TrianglesMeshReader<1,2>(
                        "mesh/test/data/circle_outline");
        TS_ASSERT_EQUALS( pMeshReader->GetNumNodes(), 100);
        TS_ASSERT_EQUALS( pMeshReader->GetNumElements(), 100);
        delete pMeshReader;
        //Note: don't test faces (end nodes), since they are culled later
        pMeshReader=new TrianglesMeshReader<1,2>(
                        "mesh/test/data/semicircle_outline");
        TS_ASSERT_EQUALS( pMeshReader->GetNumNodes(), 51);
        TS_ASSERT_EQUALS( pMeshReader->GetNumElements(), 50);
        delete pMeshReader;
    }
    
    void Test2DMeshIn3DSpace()
    {
    
        AbstractMeshReader<2,3> *pMeshReader;
        
        pMeshReader=new TrianglesMeshReader<2,3>(
                        "mesh/test/data/slab_395_elements");
                        
                        
        TS_ASSERT_EQUALS( pMeshReader->GetNumNodes(), 132);
        
        TS_ASSERT_EQUALS( pMeshReader->GetNumElements(), 224);
        
        TS_ASSERT_EQUALS( pMeshReader->GetNumFaces(), 0);
        
        delete pMeshReader;
        
        pMeshReader=new TrianglesMeshReader<2,3>(
                        "mesh/test/data/disk_in_3d");
                        
                        
        TS_ASSERT_EQUALS( pMeshReader->GetNumNodes(), 312);
        
        TS_ASSERT_EQUALS( pMeshReader->GetNumElements(), 522);
        
        TS_ASSERT_EQUALS( pMeshReader->GetNumFaces(), 833);
        
        delete pMeshReader;
        
    }
    
    void TestExceptions()
    {
        // this should fail because a 3D reader is attempting to read a 2D mesh
        TS_ASSERT_THROWS_ANYTHING( READER_3D mesh_reader("mesh/test/data/disk_984_elements") );
    }
    
};

#endif //_TESTTRIANGLESMESHREADER_HPP_
