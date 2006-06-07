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
#include "OutputFileHandler.hpp"
#include <cmath>
#include <sys/stat.h>

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
                            "", "MeshFromMemfem",3);
//                            "/tmp/testoutput/MeshFromMemfem",3);
	
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
								new TrianglesMeshReader("/tmp/testoutput/MeshFromMemfem"));
								
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
							"","MeshFromFemlab",2);
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
							new TrianglesMeshReader("/tmp/testoutput/MeshFromFemlab"));
							
		delete pImportMeshReader;
		delete pMeshWriter;
		delete pNewMeshReader;							
	}
	
	void TestTrianglesToMeshalyzer(void)
	{	
		pImportMeshReader=new TrianglesMeshReader(
							"mesh/test/data/slab_138_elements");
		pMeshWriter=new MeshalyzerMeshWriter(
							"", "MeshFromTetgen");
	
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
		char *fake_data_name = "/tmp/testoutput/MeshFromTetgen.tdat";
		std::ofstream fake_data(fake_data_name);
        TS_ASSERT(fake_data.is_open());
        chmod(fake_data_name, 0666);
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

//        std::string output_directory = "/tmp/testoutput/CGFromTetgen/";
//        system(("mkdir -p "+output_directory).c_str());
//        chmod(output_directory.c_str(), 0777);
            
		pMeshWriter=new MeshalyzerMeshWriter("CGFromTetgen", "CGFromTetgen", set_CG_format);
	
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

        OutputFileHandler output_file_handler("CGFromTetgen");
		for(int t= 0; t<num_tsteps ;t++)
		{
			char fake_data_name[100];
			sprintf(fake_data_name, std::string("CGFromTetgen.t%i").c_str(), t);
            std::ofstream *p_fake_data = output_file_handler.OpenOutputFile(fake_data_name);
//			std::ofstream fake_data(fake_data_name) ;
//            TS_ASSERT(fake_data->is_open());
//            chmod(fake_data_name, 0666);
			*p_fake_data << "t = "<<t<<"\n";
			for(int n=0 ; n<num_nodes ; n++)
			{
				*p_fake_data<<0.0<<"\t"<<0.0<<"\t"<<sin((t*M_PI)/num_tsteps) <<"\n";
			}
		 
			p_fake_data->close() ;
            delete p_fake_data;
		}
		
		delete pImportMeshReader;
		delete pMeshWriter;

	}	
};

#endif //_TESTMEMFEMMESHREADER_HPP_
