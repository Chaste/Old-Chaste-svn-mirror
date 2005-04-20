#ifndef _TESTTRIANGLESMESHREADER_HPP_
#define _TESTTRIANGLESMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "../AbstractMeshReader.hpp"

class TestTrianglesMeshReaders : public CxxTest::TestSuite
{
	public:
	void testFilesOpen(void)
	{
		bool success;
		
		//TriangleMeshReader meshReader("disk_522_elements");
		AbstractMeshReader *pmeshReader=new AbstractMeshReader("disk_522_elements");
		
		//success=meshReader.IsReaderSuccess();
		//TS_ASSERT(success);
		
		success=pmeshReader->IsReaderSuccess();
		
		TS_ASSERT(success);
	}
};

#endif //_TESTTRIANGLESMESHREADER_HPP_
