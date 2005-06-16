// TestTrianglesMeshReader.hpp

/**
 * Test suite for the TrianglesMeshReader class.
 * 
 */

#ifndef _TESTTRIANGLESMESHREADER_HPP_
#define _TESTTRIANGLESMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.hpp"


class TestTrianglesMeshReaders : public CxxTest::TestSuite
{
    private:
          AbstractMeshReader *spMeshReader;
	public:
	
	/**
	 * Check that input files are opened correctly.
	 * 
	 */
	
	void testFilesOpen(void)
	{
        AbstractMeshReader *pMeshReader;
		TS_ASSERT_THROWS_NOTHING(
		                  pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements"));
        delete pMeshReader;                          
		
	}
	
	
	
	/**
	 * Check that large data files can be opened without problems. This
	 * currently causes the whole system to slow down.
	 * 
	 */
	
	//void testTulaneFilesOpen(void)
	//{
        //AbstractMeshReader *pMeshReader;
		//TS_ASSERT_THROWS_NOTHING(
		//                  pMeshReader=new TrianglesMeshReader(
		//                 "pdes/tests/meshdata/tulane_data_about_400k_elements"));
        //delete pMeshReader;
		
	//}
	
	
	/**
	 * Check that the nodes are read correctly. Checks that the output vector
	 * for a given input file is the correct length and that if the input file
	 * is corrupted (missing nodes) then an exception is thrown.
	 * 
	 */
	
	void TestNodesDataRead(void)
	{
        AbstractMeshReader *pMeshReader;
		pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements_indexed_from_1");
		
		TS_ASSERT( pMeshReader->GetNumNodes() == 543); 
        delete pMeshReader;
		
		TS_ASSERT_THROWS_ANYTHING(
		                  pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/baddata/bad_nodes_disk_522__elements_indexed_from_1"));		
	}
	
	
	
	/**
	 * Check that the elements are read correctly. Checks that the output vector
	 * for a given input file is the correct length and that if the input file
	 * is corrupted (missing elements) then an exception is thrown.
	 * 
	 */
	
	void TestElementsDataRead(void)
	{
        AbstractMeshReader *pMeshReader;
		pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements_indexed_from_1");
		
		TS_ASSERT( pMeshReader->GetNumElements() == 984); 
        delete pMeshReader;
		
		TS_ASSERT_THROWS_ANYTHING(
		                  pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/baddata/bad_elements_disk_522_elements_indexed_from_1"));

	}
	
	
	
	/**
	 * Check that the faces are read correctly. Checks that the output vector
	 * for a given input file is the correct length and that if the input file
	 * is corrupted (missing faces) then an exception is thrown.
	 * 
	 */
	
	void TestFacesDataRead(void)
	{
        AbstractMeshReader *pMeshReader;
		pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements_indexed_from_1");
		
		TS_ASSERT( pMeshReader->GetNumFaces() == 1526); 
		delete pMeshReader;
		
		TS_ASSERT_THROWS_ANYTHING(
		                  pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/baddata/bad_faces_disk_522__elements_indexed_from_1"));		
	}
	
	
	
	/**
	 * Checks that the reader can deal with (3-d) TetGen input files as well
	 * as the previously considered (2-d) Triangles files. Checks that the
	 * element output vector for a given input file is the correct length.
	 * 
	 */
	
	void Test3dDataRead(void)
	{			
        AbstractMeshReader *pMeshReader;
		TS_ASSERT_THROWS_NOTHING(
		                  pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/slab_138_elements"));
			
		TS_ASSERT (pMeshReader->GetNumElements() == 138);
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
        AbstractMeshReader *pMeshReader;
		TS_ASSERT_THROWS_NOTHING(
		                  pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements"));
		
		TS_ASSERT(pMeshReader->GetMaxNodeIndex() == pMeshReader->GetNumNodes() - 1);
		
		TS_ASSERT(pMeshReader->GetMinNodeIndex() == 0);
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
        AbstractMeshReader *pMeshReader;
		TS_ASSERT_THROWS_NOTHING(
		                  pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements_indexed_from_1"));
		
		TS_ASSERT(pMeshReader->GetMaxNodeIndex() == pMeshReader->GetNumNodes() - 1);
		
		TS_ASSERT(pMeshReader->GetMinNodeIndex() == 0);
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
        AbstractMeshReader *pMeshReader;
		TS_ASSERT_THROWS_ANYTHING(
		                  pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/baddata/permuted_nodes_disk_522_elements"));	

	}
	
	
	/**
	 * Checks that elements have the correct number of nodes (i.e. one more
	 * node than the dimension of the mesh). If quadratic basis functions are
	 * required this should be dealt with elsewhere. 
	 * 
	 */
	
	void TestOrder2ElementsFail(void)
	{
        AbstractMeshReader *pMeshReader;
		TS_ASSERT_THROWS_ANYTHING(
		                  pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/baddata/disk_522_order_2_elements"));	
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
        AbstractMeshReader *pMeshReader;
		pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements");
		
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
        AbstractMeshReader *pMeshReader;
		pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements");
		
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
        AbstractMeshReader *pMeshReader;
		pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements");
		
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
	
	void TestGetNextBoundaryEdge(void)
	{
        AbstractMeshReader *pMeshReader;
		pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements");
		
		std::vector<int> NextEdge;
		
		TS_ASSERT(pMeshReader->GetNumBoundaryEdges() == 100);
		
		TS_ASSERT_THROWS_NOTHING(NextEdge = pMeshReader->GetNextBoundaryFace());
		TS_ASSERT_THROWS_NOTHING(NextEdge = pMeshReader->GetNextBoundaryFace());
		    		
		for (int i = 2; i < pMeshReader->GetNumBoundaryEdges(); i++)
		{
			TS_ASSERT_THROWS_NOTHING(NextEdge = pMeshReader->GetNextBoundaryEdge());
		}
		
		TS_ASSERT_THROWS_ANYTHING(NextEdge = pMeshReader->GetNextBoundaryEdge());
		delete pMeshReader;
	}	

	/**
	 * Check that the 1D data are read correctly. Check that the output vector
	 * for a given input file is the correct length. Check that the number of 
	 * boundary faces is 2.
	 */
	
	void Test1DMeshRead(void)
	{
        AbstractMeshReader *pMeshReader;
		pMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/trivial_1d_mesh");
		
		TS_ASSERT( pMeshReader->GetNumNodes() == 11); 
		
		TS_ASSERT( pMeshReader->GetNumElements() == 10); 
		
		TS_ASSERT( pMeshReader->GetNumFaces() == 11); 		
		
		TS_ASSERT( pMeshReader->GetNumBoundaryFaces() == 2); 		
	    
        delete pMeshReader;
       
	}

};

#endif //_TESTTRIANGLESMESHREADER_HPP_
