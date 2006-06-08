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
//#include <iostream>

AbstractMeshReader *pImportMeshReader;
AbstractMeshReader *pNewMeshReader;
AbstractMeshWriter *pMeshWriter;

class TestMeshWriters : public CxxTest::TestSuite
{
    public:
        
    void TestMemfemtoTetgen(void)
    {
        pImportMeshReader=new MemfemMeshReader("mesh/test/data/Memfem_slab");
        pMeshWriter = new TrianglesMeshWriter("", "MeshFromMemfem",3);
    
        for (int i=0; i<pImportMeshReader->GetNumNodes();i++)
        {
            pMeshWriter->SetNextNode(pImportMeshReader->GetNextNode());
        }
        for (int i=0; i<pImportMeshReader->GetNumElements();i++)
        {
            pMeshWriter->SetNextElement(pImportMeshReader->GetNextElement());
        }
        // Note: the results of this may not be as expected!
        for (int i=0; i<pImportMeshReader->GetNumFaces();i++)
        {
            pMeshWriter->SetNextBoundaryFace(pImportMeshReader->GetNextFace());
        }
        
        pMeshWriter->WriteFiles();
        std::string output_dir = pMeshWriter->GetOutputDirectory();
        
        TS_ASSERT_THROWS_NOTHING(pNewMeshReader = 
            new TrianglesMeshReader(output_dir + "MeshFromMemfem"));
                                
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
        std::string output_dir = pMeshWriter->GetOutputDirectory();

        TS_ASSERT_THROWS_NOTHING(pNewMeshReader = 
            new TrianglesMeshReader(output_dir + "MeshFromFemlab"));
                            
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
        
        // Fake data on the above mesh
        // For use in testing by eye
//        int num_tsteps=500;
//        int num_nodes = pImportMeshReader->GetNumNodes();
//        OutputFileHandler output_file_handler("");
//        out_stream p_fake_data = output_file_handler.OpenOutputFile(
//            "MeshFromTetgen.tdat");
//        for (int t= 0; t<num_tsteps ;t++)
//        {
//            for (int n=0 ; n<num_nodes ; n++)
//            {
//                *p_fake_data<<sin((t*M_PI)/num_tsteps) <<"\n";
//            }
//        } 
//        p_fake_data->close();
        
        delete pImportMeshReader;
        delete pMeshWriter;
    }

    void TestTrianglesToCoolGraphics(void)
    {    
        pImportMeshReader=new TrianglesMeshReader(
                            "mesh/test/data/slab_138_elements");
        bool set_CG_format=true;
            
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
                
        // Fake data on the above mesh
        // For use in testing by eye
//        int num_tsteps=500;
//        int num_nodes = pImportMeshReader->GetNumNodes();
//        OutputFileHandler output_file_handler("CGFromTetgen");
//        for(int t= 0; t<num_tsteps ;t++)
//        {
//            char fake_data_name[100];
//            sprintf(fake_data_name, std::string("CGFromTetgen.t%i").c_str(), t);
//            out_stream p_fake_data = output_file_handler.OpenOutputFile(fake_data_name);
//            *p_fake_data << "t = "<<t<<"\n";
//            for(int n=0 ; n<num_nodes ; n++)
//            {
//                *p_fake_data<<0.0<<"\t"<<0.0<<"\t"<<sin((t*M_PI)/num_tsteps) <<"\n";
//            }
//         
//            p_fake_data->close() ;
//        }
        
        delete pImportMeshReader;
        delete pMeshWriter;

    }
};

#endif //_TESTMEMFEMMESHREADER_HPP_
