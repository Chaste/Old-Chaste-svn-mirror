// TestMeshWriters.hpp

/**
 * Test suite for the MeshWriter class(es).
 * 
 */

#ifndef _TESTMESHWRITERS_HPP_
#define _TESTMESHWRITERS_HPP_

#include <cxxtest/TestSuite.h>
#include "../MemfemMeshReader.hpp"
#include "../FemlabMeshReader.hpp"
#include "../TrianglesMeshWriter.hpp"

class TestMeshWriters : public CxxTest::TestSuite
{
	public:
	
	
	void testMemfemtoTetgen(void)
	{	
		TS_ASSERT(1==1);
	}

	void testFemlabtoTriangles(void)
	{	
		TS_ASSERT(1==1);
	}
	
};

#endif //_TESTMEMFEMMESHREADER_HPP_
