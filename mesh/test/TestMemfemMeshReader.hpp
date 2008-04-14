/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

// TestMemfemMeshReader.hpp

/**
 * Test suite for the MemfemMeshReader class.
 *
 */

#ifndef _TESTMEMFEMMESHREADER_HPP_
#define _TESTMEMFEMMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "MemfemMeshReader.hpp"

typedef MemfemMeshReader<3,3> READER_3D;

class TestMemfemMeshReaders : public CxxTest::TestSuite
{
public:

    /**
     * Check that input files are opened correctly.
     * 
     */
    void TestFilesOpen(void)
    {
        MemfemMeshReader<3,3> *pMeshReader;
        pMeshReader = new READER_3D(
                              "mesh/test/data/Memfem_slab");
                              
        TS_ASSERT(pMeshReader->GetNumNodes() == 381U);
        TS_ASSERT(pMeshReader->GetNumElements() == 1030U);
        TS_ASSERT(pMeshReader->GetNumFaces() == 758U);
        
        std::vector<unsigned> NextFace;
        
        NextFace = pMeshReader->GetNextFace();
        
        TS_ASSERT( NextFace[0] == 338  );
        TS_ASSERT( NextFace[1] == 23 );
        TS_ASSERT( NextFace[2] == 374 );
        
        TS_ASSERT(pMeshReader->GetMaxNodeIndex() == pMeshReader->GetNumNodes() - 1);
        
        TS_ASSERT(pMeshReader->GetMinNodeIndex() == 0);
        
        delete pMeshReader;
    }
    
};

#endif //_TESTMEMFEMMESHREADER_HPP_
