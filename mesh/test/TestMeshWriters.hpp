#ifndef _TESTMESHWRITERS_HPP_
#define _TESTMESHWRITERS_HPP_
/**
 * Test suite for the MeshWriter class(es).
 * 
 */


#include <cxxtest/TestSuite.h>
#include "MemfemMeshReader.hpp"
#include "FemlabMeshReader.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include <math.h>

AbstractMeshReader *pImportMeshReader;
AbstractMeshReader *pNewMeshReader;
AbstractMeshWriter *pMeshWriter;

class TestMeshWriters : public CxxTest::TestSuite
{
	public:
		
	void TestMemfemtoTetgen(void)
	{	
		pImportMeshReader=new MemfemMeshReader(
							"mesh/test/data/Memfem_slab");
		pMeshWriter=new TrianglesMeshWriter(
							"testoutput/MeshFromMemfem",3);
	
		int i;
		for (i=0; i<pImportMeshReader->GetNumNodes();i++)
		{
			pMeshWriter->SetNextNode(pImportMeshReader->GetNextNode());
		}
		for (i=0; i<pImportMeshReader->GetNumElements();i++)
		{
			pMeshWriter->SetNextElement(pImportMeshReader->GetNextElement());
		}
        // Note: the results of this may not be as expected!
		for (i=0; i<pImportMeshReader->GetNumFaces();i++)
		{
			pMeshWriter->SetNextBoundaryFace(pImportMeshReader->GetNextFace());
		}
		
		pMeshWriter->WriteFiles();
		
		TS_ASSERT_THROWS_NOTHING(pNewMeshReader = 
								new TrianglesMeshReader("testoutput/MeshFromMemfem"));
								
		delete pImportMeshReader;
		delete pMeshWriter;
		delete pNewMeshReader;
	
	}

	void TestFemlabtoTriangles(void)
	{	
		pImportMeshReader=new FemlabMeshReader(
							"mesh/test/data/",
		                  	"femlab_lshape_nodes.dat",
		                  	"femlab_lshape_elements.dat",
		                  	"femlab_lshape_edges.dat");
		pMeshWriter=new TrianglesMeshWriter(
							"testoutput/MeshFromFemlab",2);
		int i;
		for (i=0; i<pImportMeshReader->GetNumNodes();i++)
		{
			pMeshWriter->SetNextNode(pImportMeshReader->GetNextNode());
		}
		for (i=0; i<pImportMeshReader->GetNumElements();i++)
		{
			pMeshWriter->SetNextElement(pImportMeshReader->GetNextElement());
		}
		for (i=0; i<pImportMeshReader->GetNumFaces();i++)
		{
			pMeshWriter->SetNextBoundaryFace(pImportMeshReader->GetNextFace());
		}
		
		pMeshWriter->WriteFiles();

		TS_ASSERT_THROWS_NOTHING(pNewMeshReader = 
							new TrianglesMeshReader("testoutput/MeshFromFemlab"));
							
		delete pImportMeshReader;
		delete pMeshWriter;
		delete pNewMeshReader;							
	}
	
	void TestTrianglesToMeshalyzer(void)
	{	
		pImportMeshReader=new TrianglesMeshReader(
							"mesh/test/data/slab_138_elements");
		pMeshWriter=new MeshalyzerMeshWriter(
							"testoutput/MeshFromTetgen");
	
		int i;
		for (i=0; i<pImportMeshReader->GetNumNodes();i++)
		{
			pMeshWriter->SetNextNode(pImportMeshReader->GetNextNode());
		}
		for (i=0; i<pImportMeshReader->GetNumElements();i++)
		{
			pMeshWriter->SetNextElement(pImportMeshReader->GetNextElement());
		}
		for (i=0; i<pImportMeshReader->GetNumFaces();i++)
		{
			pMeshWriter->SetNextBoundaryFace(pImportMeshReader->GetNextFace());
		}
				
		TS_ASSERT_THROWS_NOTHING(pMeshWriter->WriteFiles());
		
		int num_tsteps=500;
		int num_nodes = pImportMeshReader->GetNumNodes();
		char fake_data_name[40];
		sprintf(fake_data_name, "testoutput/MeshFromTetgen.tdat");
		std::ofstream fake_data(fake_data_name) ;
		for(int t= 0; t<num_tsteps ;t++)
		{
			for(int n=0 ; n<num_nodes ; n++)
			{
				fake_data<<sin((t*M_PI)/num_tsteps) <<"\n";
			}
		} 
		fake_data.close() ;		
		
		delete pImportMeshReader;
		delete pMeshWriter;

	}

	void TestTrianglesToCoolGraphics(void)
	{	
		pImportMeshReader=new TrianglesMeshReader(
							"mesh/test/data/slab_138_elements");
		bool set_CG_format=true;
		pMeshWriter=new MeshalyzerMeshWriter(
							"testoutput/CGFromTetgen", set_CG_format);
	
		int i;
		for (i=0; i<pImportMeshReader->GetNumNodes();i++)
		{
			pMeshWriter->SetNextNode(pImportMeshReader->GetNextNode());
		}
		for (i=0; i<pImportMeshReader->GetNumElements();i++)
		{
			pMeshWriter->SetNextElement(pImportMeshReader->GetNextElement());
		}
		for (i=0; i<pImportMeshReader->GetNumFaces();i++)
		{
			pMeshWriter->SetNextBoundaryFace(pImportMeshReader->GetNextFace());
		}				
		
		TS_ASSERT_THROWS_NOTHING(pMeshWriter->WriteFiles());
				
		int num_tsteps=500;
		int num_nodes = pImportMeshReader->GetNumNodes();
		for(int t= 0; t<num_tsteps ;t++)
		{
			char fake_data_name[40];
			sprintf(fake_data_name, "testoutput/CGFromTetgen.t%i", t);
			std::ofstream fake_data(fake_data_name) ;
			fake_data << "t = "<<t<<"\n";
			for(int n=0 ; n<num_nodes ; n++)
			{
				fake_data<<0.0<<"\t"<<0.0<<"\t"<<sin((t*M_PI)/num_tsteps) <<"\n";
			}
		 
			fake_data.close() ;
		}
		
		delete pImportMeshReader;
		delete pMeshWriter;

	}	
};

#endif //_TESTMEMFEMMESHREADER_HPP_
