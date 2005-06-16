// TestMeshWriters.hpp

/**
 * Test suite for the MeshWriter class(es).
 * 
 */

#ifndef _TESTMESHWRITERS_HPP_
#define _TESTMESHWRITERS_HPP_

#include <cxxtest/TestSuite.h>
#include "MemfemMeshReader.hpp"
#include "FemlabMeshReader.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include <math.h>

static AbstractMeshReader *spImportMeshReader;
static AbstractMeshReader *spNewMeshReader;
static AbstractMeshWriter *spMeshWriter;

class TestMeshWriters : public CxxTest::TestSuite
{
	public:
		
	void TestMemfemtoTetgen(void)
	{	
		spImportMeshReader=new MemfemMeshReader(
							"mesh/test/data/Memfem_slab");
		spMeshWriter=new TrianglesMeshWriter(
							"testoutput/MeshFromMemfem",3);
	
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
								new TrianglesMeshReader("testoutput/MeshFromMemfem"));
	
	}

	void TestFemlabtoTriangles(void)
	{	
		spImportMeshReader=new FemlabMeshReader(
							"mesh/test/data/",
		                  	"femlab_lshape_nodes.dat",
		                  	"femlab_lshape_elements.dat",
		                  	"femlab_lshape_edges.dat");
		spMeshWriter=new TrianglesMeshWriter(
							"testoutput/MeshFromFemlab",2);
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
							new TrianglesMeshReader("testoutput/MeshFromFemlab"));
	}
	
	void TestTrianglesToMeshalyzer(void)
	{	
		spImportMeshReader=new TrianglesMeshReader(
							"mesh/test/data/slab_138_elements");
		spMeshWriter=new MeshalyzerMeshWriter(
							"testoutput/MeshFromTetgen");
	
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
				
		TS_ASSERT_THROWS_NOTHING(spMeshWriter->WriteFiles());
		
		int num_tsteps=500;
		int num_nodes = spImportMeshReader->GetNumNodes();
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
	}

	void TestTrianglesToCoolGraphics(void)
	{	
		spImportMeshReader=new TrianglesMeshReader(
							"mesh/test/data/slab_138_elements");
		bool set_CG_format=true;
		spMeshWriter=new MeshalyzerMeshWriter(
							"testoutput/CGFromTetgen", set_CG_format);
	
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
		
		TS_ASSERT_THROWS_NOTHING(spMeshWriter->WriteFiles());
				
		int num_tsteps=500;
		int num_nodes = spImportMeshReader->GetNumNodes();
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
	}	
};

#endif //_TESTMEMFEMMESHREADER_HPP_
