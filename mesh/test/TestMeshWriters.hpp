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
#include "TrianglesMeshWriter.cpp"
#include "MeshalyzerMeshWriter.cpp"
#include "OutputFileHandler.hpp"
#include <cmath>
//#include <iostream>

typedef TrianglesMeshReader<3,3> TRI_READ_3;
typedef TrianglesMeshReader<2,2> TRI_READ_2;

class TestMeshWriters : public CxxTest::TestSuite
{
    public:
        
    void TestMemfemtoTetgen(void)
    {
        TrianglesMeshWriter<3,3> mesh_writer("", "MeshFromMemfem",3);
        MemfemMeshReader<3,3> import_mesh_reader("mesh/test/data/Memfem_slab");
    
        for (int i=0; i<import_mesh_reader.GetNumNodes();i++)
        {
            mesh_writer.SetNextNode(import_mesh_reader.GetNextNode());
        }
        for (int i=0; i<import_mesh_reader.GetNumElements();i++)
        {
            mesh_writer.SetNextElement(import_mesh_reader.GetNextElement());
        }
        // Note: the results of this may not be as expected!
        for (int i=0; i<import_mesh_reader.GetNumFaces();i++)
        {
            mesh_writer.SetNextBoundaryFace(import_mesh_reader.GetNextFace());
        }
        
        mesh_writer.WriteFiles();
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TRI_READ_3 *p_new_mesh_reader;
        TS_ASSERT_THROWS_NOTHING(p_new_mesh_reader = 
            new TRI_READ_3(output_dir + "MeshFromMemfem"));
                                
        delete p_new_mesh_reader;
    }

    void TestFemlabtoTriangles(void)
    {   
        TrianglesMeshWriter<2,2> mesh_writer("","MeshFromFemlab",2);

        FemlabMeshReader<2,2> import_mesh_reader(
                              "mesh/test/data/",
                              "femlab_lshape_nodes.dat",
                              "femlab_lshape_elements.dat",
                              "femlab_lshape_edges.dat");
        int i;
        for (i=0; i<import_mesh_reader.GetNumNodes();i++)
        {
            mesh_writer.SetNextNode(import_mesh_reader.GetNextNode());
        }
        for (i=0; i<import_mesh_reader.GetNumElements();i++)
        {
            mesh_writer.SetNextElement(import_mesh_reader.GetNextElement());
        }
        for (i=0; i<import_mesh_reader.GetNumEdges();i++)
        {
            mesh_writer.SetNextBoundaryEdge(import_mesh_reader.GetNextEdge());
        }
        
        mesh_writer.WriteFiles();
        std::string output_dir = mesh_writer.GetOutputDirectory();

        TRI_READ_2 *p_new_mesh_reader;
        TS_ASSERT_THROWS_NOTHING(p_new_mesh_reader = 
            new TRI_READ_2(output_dir + "MeshFromFemlab"));
                            
        delete p_new_mesh_reader;                            
    }
    
    void TestTrianglesToMeshalyzer(void)
    {    
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        MeshalyzerMeshWriter<3,3> mesh_writer("", "MeshFromTetgen");
    
        int i;
        for (i=0; i<import_mesh_reader.GetNumNodes();i++)
        {
            mesh_writer.SetNextNode(import_mesh_reader.GetNextNode());
        }
        for (i=0; i<import_mesh_reader.GetNumElements();i++)
        {
            mesh_writer.SetNextElement(import_mesh_reader.GetNextElement());
        }
        for (i=0; i<import_mesh_reader.GetNumFaces();i++)
        {
            mesh_writer.SetNextBoundaryFace(import_mesh_reader.GetNextFace());
        }
                
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFiles());
    }

    void TestTrianglesToCoolGraphics(void)
    {    
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        bool set_CG_format=true;
            
        MeshalyzerMeshWriter<3,3> mesh_writer("CGFromTetgen", "CGFromTetgen", set_CG_format);
    
        int i;
        for (i=0; i<import_mesh_reader.GetNumNodes();i++)
        {
            mesh_writer.SetNextNode(import_mesh_reader.GetNextNode());
        }
        for (i=0; i<import_mesh_reader.GetNumElements();i++)
        {
            mesh_writer.SetNextElement(import_mesh_reader.GetNextElement());
        }
        for (i=0; i<import_mesh_reader.GetNumFaces();i++)
        {
            mesh_writer.SetNextBoundaryFace(import_mesh_reader.GetNextFace());
        }                
        
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFiles());
    }
};

#endif //_TESTMEMFEMMESHREADER_HPP_
