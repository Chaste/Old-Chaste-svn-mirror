#ifndef _TESTMESHWRITERS_HPP_
#define _TESTMESHWRITERS_HPP_
/**
 * Test suite for the MeshWriter class(es).
 * 
 */


#include <cxxtest/TestSuite.h>
#include "MemfemMeshReader.cpp"
#include "FemlabMeshReader.cpp"
#include "TrianglesMeshReader.cpp"
#include "TrianglesMeshWriter.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "OutputFileHandler.hpp"
#include <cmath>
//#include <iostream>

AbstractMeshWriter *pMeshWriter;
typedef TrianglesMeshReader<3,3> TRI_READ_3;
typedef TrianglesMeshReader<2,2> TRI_READ_2;

class TestMeshWriters : public CxxTest::TestSuite
{
    public:
        
    void TestMemfemtoTetgen(void)
    {
        MemfemMeshReader<3,3>* p_import_mesh_reader=new MemfemMeshReader<3,3>("mesh/test/data/Memfem_slab");
        pMeshWriter = new TrianglesMeshWriter("", "MeshFromMemfem",3);
    
        for (int i=0; i<p_import_mesh_reader->GetNumNodes();i++)
        {
            pMeshWriter->SetNextNode(p_import_mesh_reader->GetNextNode());
        }
        for (int i=0; i<p_import_mesh_reader->GetNumElements();i++)
        {
            pMeshWriter->SetNextElement(p_import_mesh_reader->GetNextElement());
        }
        // Note: the results of this may not be as expected!
        for (int i=0; i<p_import_mesh_reader->GetNumFaces();i++)
        {
            pMeshWriter->SetNextBoundaryFace(p_import_mesh_reader->GetNextFace());
        }
        
        pMeshWriter->WriteFiles();
        std::string output_dir = pMeshWriter->GetOutputDirectory();
        TRI_READ_3 *p_new_mesh_reader;
        TS_ASSERT_THROWS_NOTHING(p_new_mesh_reader = 
            new TRI_READ_3(output_dir + "MeshFromMemfem"));
                                
        delete p_import_mesh_reader;
        delete pMeshWriter;
        delete p_new_mesh_reader;
    }

    void TestFemlabtoTriangles(void)
    {    
        FemlabMeshReader<2,2>* p_import_mesh_reader=new FemlabMeshReader<2,2>(
                            "mesh/test/data/",
                              "femlab_lshape_nodes.dat",
                              "femlab_lshape_elements.dat",
                              "femlab_lshape_edges.dat");
        pMeshWriter=new TrianglesMeshWriter(
                            "","MeshFromFemlab",2);
        int i;
        for (i=0; i<p_import_mesh_reader->GetNumNodes();i++)
        {
            pMeshWriter->SetNextNode(p_import_mesh_reader->GetNextNode());
        }
        for (i=0; i<p_import_mesh_reader->GetNumElements();i++)
        {
            pMeshWriter->SetNextElement(p_import_mesh_reader->GetNextElement());
        }
        for (i=0; i<p_import_mesh_reader->GetNumEdges();i++)
        {
            pMeshWriter->SetNextBoundaryEdge(p_import_mesh_reader->GetNextEdge());
        }
        
        pMeshWriter->WriteFiles();
        std::string output_dir = pMeshWriter->GetOutputDirectory();

        TRI_READ_2 *p_new_mesh_reader;
        TS_ASSERT_THROWS_NOTHING(p_new_mesh_reader = 
            new TRI_READ_2(output_dir + "MeshFromFemlab"));
                            
        delete p_import_mesh_reader;
        delete pMeshWriter;
        delete p_new_mesh_reader;                            
    }
    
    void TestTrianglesToMeshalyzer(void)
    {    
        TrianglesMeshReader<3,3>* p_import_mesh_reader=new TrianglesMeshReader<3,3>(
                            "mesh/test/data/slab_138_elements");
        pMeshWriter=new MeshalyzerMeshWriter(
                            "", "MeshFromTetgen");
    
        int i;
        for (i=0; i<p_import_mesh_reader->GetNumNodes();i++)
        {
            pMeshWriter->SetNextNode(p_import_mesh_reader->GetNextNode());
        }
        for (i=0; i<p_import_mesh_reader->GetNumElements();i++)
        {
            pMeshWriter->SetNextElement(p_import_mesh_reader->GetNextElement());
        }
        for (i=0; i<p_import_mesh_reader->GetNumFaces();i++)
        {
            pMeshWriter->SetNextBoundaryFace(p_import_mesh_reader->GetNextFace());
        }
                
        TS_ASSERT_THROWS_NOTHING(pMeshWriter->WriteFiles());
        
        // Fake data on the above mesh
        // For use in testing by eye
//        int num_tsteps=500;
//        int num_nodes = p_import_mesh_reader->GetNumNodes();
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
        
        delete p_import_mesh_reader;
        delete pMeshWriter;
    }

    void TestTrianglesToCoolGraphics(void)
    {    
        TrianglesMeshReader<3,3>* p_import_mesh_reader=new TrianglesMeshReader<3,3>(
                            "mesh/test/data/slab_138_elements");
        bool set_CG_format=true;
            
        pMeshWriter=new MeshalyzerMeshWriter("CGFromTetgen", "CGFromTetgen", set_CG_format);
    
        int i;
        for (i=0; i<p_import_mesh_reader->GetNumNodes();i++)
        {
            pMeshWriter->SetNextNode(p_import_mesh_reader->GetNextNode());
        }
        for (i=0; i<p_import_mesh_reader->GetNumElements();i++)
        {
            pMeshWriter->SetNextElement(p_import_mesh_reader->GetNextElement());
        }
        for (i=0; i<p_import_mesh_reader->GetNumFaces();i++)
        {
            pMeshWriter->SetNextBoundaryFace(p_import_mesh_reader->GetNextFace());
        }                
        
        TS_ASSERT_THROWS_NOTHING(pMeshWriter->WriteFiles());
                
        // Fake data on the above mesh
        // For use in testing by eye
//        int num_tsteps=500;
//        int num_nodes = p_import_mesh_reader->GetNumNodes();
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
        
        delete p_import_mesh_reader;
        delete pMeshWriter;

    }
};

#endif //_TESTMEMFEMMESHREADER_HPP_
