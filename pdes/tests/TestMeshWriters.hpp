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
#include "../TrianglesMeshReader.hpp"
#include "../TrianglesMeshWriter.hpp"

static AbstractMeshReader *spImportMeshReader;
static AbstractMeshReader *spNewMeshReader;
static TrianglesMeshWriter *spMeshWriter;

class TestMeshWriters : public CxxTest::TestSuite
{
	public:
	
	
	void testMemfemtoTetgen(void)
	{	
		spImportMeshReader=new MemfemMeshReader(
							"pdes/tests/meshdata/Memfem_slab");
		spMeshWriter=new TrianglesMeshWriter(
							"/tmp/MeshFromMemfem",3);
	
		int i;
		for (i=0; i<spImportMeshReader->GetNumNodes();i++)
		{
			spMeshWriter->SetNextNode(spImportMeshReader->GetNextNode());
		}
		for (i=0; i<spImportMeshReader->GetNumElements();i++)
		{
			spMeshWriter->SetNextElement(spImportMeshReader->GetNextElement());
		}
		for (i=0; i<spImportMeshReader->GetNumBoundaryFaces();i++)
		{
			spMeshWriter->SetNextBoundaryFace(spImportMeshReader->GetNextBoundaryFace());
		}
		
		spMeshWriter->WriteFiles();
		
		TS_ASSERT_THROWS_NOTHING(spNewMeshReader = 
								new TrianglesMeshReader("/tmp/MeshFromMemfem"));
	
	}

	void testFemlabtoTriangles(void)
	{	
		spImportMeshReader=new FemlabMeshReader(
							"pdes/tests/meshdata/",
		                  	"femlab_lshape_nodes.dat",
		                  	"femlab_lshape_elements.dat",
		                  	"femlab_lshape_edges.dat");
		spMeshWriter=new TrianglesMeshWriter(
							"/tmp/MeshFromFemlab",2);
		int i;
		for (i=0; i<spImportMeshReader->GetNumNodes();i++)
		{
			spMeshWriter->SetNextNode(spImportMeshReader->GetNextNode());
		}
		for (i=0; i<spImportMeshReader->GetNumElements();i++)
		{
			spMeshWriter->SetNextElement(spImportMeshReader->GetNextElement());
		}
		for (i=0; i<spImportMeshReader->GetNumBoundaryFaces();i++)
		{
			spMeshWriter->SetNextBoundaryFace(spImportMeshReader->GetNextBoundaryFace());
		}
		
		spMeshWriter->WriteFiles();


		TS_ASSERT_THROWS_NOTHING(spNewMeshReader = 
							new TrianglesMeshReader("/tmp/MeshFromFemlab"));
	}
	
};

#endif //_TESTMEMFEMMESHREADER_HPP_
