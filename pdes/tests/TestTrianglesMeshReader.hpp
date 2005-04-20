#ifndef _TESTTRIANGLESMESHREADER_HPP_
#define _TESTTRIANGLESMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "../TrianglesMeshReader.hpp"

class TestTrianglesMeshReaders : public CxxTest::TestSuite
{
	public:
	void testFilesOpen(void)
	{
		
		//TriangleMeshReader meshReader("disk_522_elements");
		AbstractMeshReader *pmeshReader;
		
		
		TS_ASSERT_THROWS_NOTHING(
		                  pmeshReader=new TrianglesMeshReader("pdes/tests/meshdata/disk_522_elements"));
		
	
	}
};

#endif //_TESTTRIANGLESMESHREADER_HPP_
