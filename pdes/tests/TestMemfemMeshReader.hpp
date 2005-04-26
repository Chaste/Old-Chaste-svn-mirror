// TestMemfemMeshReader.hpp

/**
 * Test suite for the MemfemMeshReader class.
 * 
 */

#ifndef _TESTMEMFEMMESHREADER_HPP_
#define _TESTMEMFEMMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "../MemfemMeshReader.hpp"

static		AbstractMeshReader *spAbMeshReader;
class TestMemfemMeshReaders : public CxxTest::TestSuite
{
	public:
	
	/**
	 * Check that input files are opened correctly.
	 * 
	 */
	
	void testFilesOpen(void)
	{
		TS_ASSERT_THROWS_NOTHING(
		                  spAbMeshReader=new MemfemMeshReader(
		                  "pdes/tests/meshdata/Memfem_slab"));
		                  
		             
		TS_ASSERT(spAbMeshReader->GetNumNodes() == 381);
		TS_ASSERT(spAbMeshReader->GetNumElements() == 1030);
		TS_ASSERT(spAbMeshReader->GetNumBoundaryFaces() == 758);
		
		std::vector<int> NextBoundaryFace;
		                  
		NextBoundaryFace = spAbMeshReader->GetNextBoundaryFace();
		
		TS_ASSERT( NextBoundaryFace[0] == 338  );
		TS_ASSERT( NextBoundaryFace[1] == 23 );
		TS_ASSERT( NextBoundaryFace[2] == 374 );
		
		TS_ASSERT(spAbMeshReader->GetMaxNodeIndex() == spAbMeshReader->GetNumNodes() - 1);
		
		TS_ASSERT(spAbMeshReader->GetMinNodeIndex() == 0);
		
			    		
	}
	
	
	};

#endif //_TESTMEMFEMMESHREADER_HPP_
