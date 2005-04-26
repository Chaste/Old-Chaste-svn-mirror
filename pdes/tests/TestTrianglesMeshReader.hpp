// TestTrianglesMeshReader.hpp

/**
 * Test suite for the TrianglesMeshReader class.
 * 
 */

#ifndef _TESTTRIANGLESMESHREADER_HPP_
#define _TESTTRIANGLESMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "../TrianglesMeshReader.hpp"

static		AbstractMeshReader *spMeshReader;
class TestTrianglesMeshReaders : public CxxTest::TestSuite
{
	public:
	
	/**
	 * Check that input files are opened correctly.
	 * 
	 */
	
	void testFilesOpen(void)
	{
		TS_ASSERT_THROWS_NOTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements"));
		
	}
	
	
	
	/**
	 * Check that large data files can be opened without problems. This
	 * currently causes the whole system to slow down.
	 * 
	 */
	
	//void testTulaneFilesOpen(void)
	//{
		//TS_ASSERT_THROWS_NOTHING(
		//                  spMeshReader=new TrianglesMeshReader(
		//                 "pdes/tests/meshdata/tulane_data_about_400k_elements"));
		
	//}
	
	
	/**
	 * Check that the nodes are read correctly. Checks that the output vector
	 * for a given input file is the correct length and that if the input file
	 * is corrupted (missing nodes) then an exception is thrown.
	 * 
	 */
	
	void testNodesDataRead(void)
	{
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements_indexed_from_1");
		
		TS_ASSERT( spMeshReader->GetNumNodes() == 543); 
		
		
		TS_ASSERT_THROWS_ANYTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/bad_nodes_disk_522__elements_indexed_from_1"));		
		
	}
	
	
	
	/**
	 * Check that the elements are read correctly. Checks that the output vector
	 * for a given input file is the correct length and that if the input file
	 * is corrupted (missing elements) then an exception is thrown.
	 * 
	 */
	
	void testElementsDataRead(void)
	{
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements_indexed_from_1");
		
		TS_ASSERT( spMeshReader->GetNumElements() == 984); 
		
		
		TS_ASSERT_THROWS_ANYTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/bad_elements_disk_522_elements_indexed_from_1"));
	
	
	}
	
	
	
	/**
	 * Check that the faces are read correctly. Checks that the output vector
	 * for a given input file is the correct length and that if the input file
	 * is corrupted (missing faces) then an exception is thrown.
	 * 
	 */
	
	void testFacesDataRead(void)
	{
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements_indexed_from_1");
		
		TS_ASSERT( spMeshReader->GetNumFaces() == 1526); 
		
		
		TS_ASSERT_THROWS_ANYTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/bad_faces_disk_522__elements_indexed_from_1"));		
		
	}
	
	
	
	/**
	 * Checks that the reader can deal with (3-d) TetGen input files as well
	 * as the previously considered (2-d) Triangles files. Checks that the
	 * element output vector for a given input file is the correct length.
	 * 
	 */
	
	void test3dDataRead(void)
	{			
		TS_ASSERT_THROWS_NOTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/slab_138_elements"));
			
		TS_ASSERT (spMeshReader->GetNumElements() == 138);
		
	}
	
	
	/**
	 * Checks that nodes are indexed from zero. Takes input file that is
	 * indexed from zero and checks that the output file also is. Uses methods
	 * GetMaxNodeIndex() and GetMinNodeIndex().
	 * 
	 */
	
	void testIndexFromZero(void)
	{		
		TS_ASSERT_THROWS_NOTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements"));
		
		TS_ASSERT(spMeshReader->GetMaxNodeIndex() == spMeshReader->GetNumNodes() - 1);
		
		TS_ASSERT(spMeshReader->GetMinNodeIndex() == 0);
		
	}
	
	
	
	/**
	 * Checks that nodes are indexed from zero. Takes input file that is
	 * indexed from one and checks that the output file also is. Uses methods
	 * GetMaxNodeIndex() and GetMinNodeIndex().
	 * 
	 */
	
	void testIndexFromOne(void)
	{
		TS_ASSERT_THROWS_NOTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements_indexed_from_1"));
		
		TS_ASSERT(spMeshReader->GetMaxNodeIndex() == spMeshReader->GetNumNodes() - 1);
		
		TS_ASSERT(spMeshReader->GetMinNodeIndex() == 0);
		
	}
	
	
	
	/**
	 * Checks that nodes in the input data file are numbered sequentially.
	 * (In the input file nodes must appear in increasing order since the node
	 * number is only stored as the index of the vector in which the coordinates
	 * are stored.)
	 * 
	 */
	
	void testPermutedNodesFail(void)
	{
		TS_ASSERT_THROWS_ANYTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/permuted_nodes_disk_522_elements"));	
		
	}
	
	
	/**
	 * Checks that elements have the correct number of nodes (i.e. one more
	 * node than the dimension of the mesh). If quadratic basis functions are
	 * required this should be dealt with elsewhere. 
	 * 
	 */
	
	void testOrder2ElementsFail(void)
	{
		TS_ASSERT_THROWS_ANYTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_order_2_elements"));	
		
	}
	
	
	/**
	 * Check that GetNextNode() returns the coordinates of the correct node.
	 * Compares the coordinates of the first two nodes with their known
	 * values, checks that no errors are thrown for the remaining nodes and
	 * that an error is thrown if we try to call the function too many times.
	 * 
	 */
	
	void testGetNextNode(void)
	{
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements");
		
		std::vector<double> FirstNode;
		                  
		FirstNode = spMeshReader->GetNextNode();
		
		TS_ASSERT_DELTA( FirstNode[0] ,  0.9980267283 , 1e-6 );
		TS_ASSERT_DELTA( FirstNode[1] , -0.0627905195 , 1e-6 )
		
		std::vector<double> NextNode;
		                  
		NextNode = spMeshReader->GetNextNode();
		
		TS_ASSERT_DELTA( NextNode[0] , 1.0 , 1e-6 );
		TS_ASSERT_DELTA( NextNode[1] , 0.0 , 1e-6 )
			    		
		for (int i = 0; i < 541; i++)
		{
			TS_ASSERT_THROWS_NOTHING(NextNode = spMeshReader->GetNextNode());
		}
		
		TS_ASSERT_THROWS_ANYTHING(NextNode = spMeshReader->GetNextNode());
		
	}
	
	
	
	/**
	 * Check that GetNextElement() works. Checks that no errors are thrown for 
	 * all of the elements and that an error is thrown if we try to call the 
	 * function too many times.
	 * 
	 */	
	
	void testGetNextElement(void)
	{
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements");
		
		std::vector<int> NextElement;
		    		
		for (int i = 0; i < spMeshReader->GetNumElements(); i++)
		{
			TS_ASSERT_THROWS_NOTHING(NextElement = spMeshReader->GetNextElement());
		}
		
		TS_ASSERT_THROWS_ANYTHING(NextElement = spMeshReader->GetNextElement());
		
	}
	
	
	/**
	 * Check that GetNextEdge() works. Checks that no errors are thrown for 
	 * all of the edges and that an error is thrown if we try to call the 
	 * function too many times.
	 * 
	 */	
	
	void testGetNextEdge(void)
	{
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements");
		
		std::vector<int> NextEdge;

		TS_ASSERT_THROWS_NOTHING(NextEdge = spMeshReader->GetNextFace());
		TS_ASSERT_THROWS_NOTHING(NextEdge = spMeshReader->GetNextFace());
		    		
		for (int i = 2; i < spMeshReader->GetNumEdges(); i++)
		{
			TS_ASSERT_THROWS_NOTHING(NextEdge = spMeshReader->GetNextEdge());
		}
		
		TS_ASSERT_THROWS_ANYTHING(NextEdge = spMeshReader->GetNextEdge());
		
	}
	
	void testGetNextBoundaryEdge(void)
	{
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements");
		
		std::vector<int> NextEdge;
		
		TS_ASSERT(spMeshReader->GetNumBoundaryEdges() == 100);
		
		TS_ASSERT_THROWS_NOTHING(NextEdge = spMeshReader->GetNextBoundaryFace());
		TS_ASSERT_THROWS_NOTHING(NextEdge = spMeshReader->GetNextBoundaryFace());
		    		
		for (int i = 2; i < spMeshReader->GetNumBoundaryEdges(); i++)
		{
			TS_ASSERT_THROWS_NOTHING(NextEdge = spMeshReader->GetNextBoundaryEdge());
		}
		
		TS_ASSERT_THROWS_ANYTHING(NextEdge = spMeshReader->GetNextBoundaryEdge());
		
	}	
};

#endif //_TESTTRIANGLESMESHREADER_HPP_
