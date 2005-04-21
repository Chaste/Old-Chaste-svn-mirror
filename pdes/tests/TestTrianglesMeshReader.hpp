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
		                  "pdes/tests/meshdata/disk_522_elements_indexed_from_1"));
		
	
	}
	
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
	
	
};

#endif //_TESTTRIANGLESMESHREADER_HPP_
