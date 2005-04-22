#ifndef _TESTTRIANGLESMESHREADER_HPP_
#define _TESTTRIANGLESMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "../TrianglesMeshReader.hpp"

static		AbstractMeshReader *spMeshReader;
class TestTrianglesMeshReaders : public CxxTest::TestSuite
{
	public:
	void testFilesOpen(void)
	{
		

		
		TS_ASSERT_THROWS_NOTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements"));
		
	
	}
	
	//void testTulaneFilesOpen(void)
	//{
		
		//std::cout<<"\nDoing a long test\n";
		
		//TS_ASSERT_THROWS_NOTHING(
		//                  spMeshReader=new TrianglesMeshReader(
		 //                 "pdes/tests/meshdata/tulane_data_about_400k_elements"));
		
		//std::cout<<"Long test finished\n";
	
	//ls
	//}
	
	void testNodesDataRead(void)
	{
		
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements_indexed_from_1");
		
		TS_ASSERT( spMeshReader->GetNumNodes() == 543); 
		
		
		TS_ASSERT_THROWS_ANYTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/bad_nodes_disk_522__elements_indexed_from_1"));		
		
	}
	
	void testElementsDataRead(void)
	{
		
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements_indexed_from_1");
		
		TS_ASSERT( spMeshReader->GetNumElements() == 984); 
		
		
		TS_ASSERT_THROWS_ANYTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/bad_elements_disk_522_elements_indexed_from_1"));
	
	
	}
	
	void testFacesDataRead(void)
	{
		
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements_indexed_from_1");
		
		TS_ASSERT( spMeshReader->GetNumFaces() == 1526); 
		
		
		TS_ASSERT_THROWS_ANYTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/bad_faces_disk_522__elements_indexed_from_1"));		
		
	}
	
	void test3dDataRead(void)
	{
		
			
		TS_ASSERT_THROWS_NOTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/slab_138_elements"));
		
		
		TS_ASSERT (spMeshReader->GetNumElements() == 138);
		
	}
	
	void testIndexFromZero(void)
	{
		
		TS_ASSERT_THROWS_NOTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements"));
		//spMeshReader=new TrianglesMeshReader("pdes/tests/meshdata/disk_522__elements");
		
		TS_ASSERT(spMeshReader->GetMaxNodeIndex() == spMeshReader->GetNumNodes() - 1);
		TS_ASSERT(spMeshReader->GetMinNodeIndex() == 0);
		
	}
	
	
	void testIndexFromOne(void)
	{
		TS_ASSERT_THROWS_NOTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements_indexed_from_1"));
		
		//spMeshReader=new TrianglesMeshReader("pdes/tests/meshdata/disk_522__elements_indexed_from_1");
		
		TS_ASSERT(spMeshReader->GetMaxNodeIndex() == spMeshReader->GetNumNodes() - 1);
		TS_ASSERT(spMeshReader->GetMinNodeIndex() == 0);
		
	}
	
	void testPermutedNodesFail(void)
	{
		
		TS_ASSERT_THROWS_ANYTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/permuted_nodes_disk_522_elements"));	
		
	}
	
	void testOrder2ElementsFail(void)
	{
		
		TS_ASSERT_THROWS_ANYTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_order_2_elements"));	
		
	}
	
	void testGetFirstNode(void)
	{
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements");
		
		std::vector<double> FirstNode;
		                  
		FirstNode = spMeshReader->GetNextNode();
		
		TS_ASSERT_DELTA( FirstNode[0] ,  0.9980267283 , 1e-6 );
		TS_ASSERT_DELTA( FirstNode[1] , -0.0627905195 , 1e-6 )
		
		
	}
	
	void testGetNextNode(void)
	{
		std::vector<double> NextNode;
		                  
		NextNode = spMeshReader->GetNextNode();
		
		TS_ASSERT_DELTA( NextNode[0] , 1.0 , 1e-6 );
		TS_ASSERT_DELTA( NextNode[1] , 0.0 , 1e-6 )
		
		
	}
	
	void testGetTooManyNodesFails(void)
	{
		std::vector<double> NextNode;
		    		
		for (int i = 0; i < 541; i++)
		{
			TS_ASSERT_THROWS_NOTHING(NextNode = spMeshReader->GetNextNode());
		}
		
		TS_ASSERT_THROWS_ANYTHING(NextNode = spMeshReader->GetNextNode());
		
	}
	
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
	
	
};

#endif //_TESTTRIANGLESMESHREADER_HPP_
