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
        TrianglesMeshWriter<3,3> mesh_writer("", "MeshFromMemfem");
        MemfemMeshReader<3,3> import_mesh_reader("mesh/test/data/Memfem_slab");
            
        mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader);
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TRI_READ_3 *p_new_mesh_reader;
        TS_ASSERT_THROWS_NOTHING(p_new_mesh_reader = 
            new TRI_READ_3(output_dir + "MeshFromMemfem"));
                                
        delete p_new_mesh_reader;
    }

    void TestFemlabtoTriangles(void)
    {   
        TrianglesMeshWriter<2,2> mesh_writer("","MeshFromFemlab");

        FemlabMeshReader<2,2> import_mesh_reader(
                              "mesh/test/data/",
                              "femlab_lshape_nodes.dat",
                              "femlab_lshape_elements.dat",
                              "femlab_lshape_edges.dat");
        
        mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader);
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
    
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }

    void TestTrianglesToCoolGraphics(void)
    {    
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        bool set_CG_format=true;
            
        MeshalyzerMeshWriter<3,3> mesh_writer("CGFromTetgen", "CGFromTetgen", set_CG_format);
    
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }
    
    
    void TestFemlabtoTrianglesViaMesh(void)
    {   
        TrianglesMeshWriter<2,2> mesh_writer("","MeshFromFemlabViaMesh");

        FemlabMeshReader<2,2> import_mesh_reader(
                              "mesh/test/data/",
                              "femlab_lshape_nodes.dat",
                              "femlab_lshape_elements.dat",
                              "femlab_lshape_edges.dat");
        
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);

        mesh_writer.WriteFilesUsingMesh(mesh);
        std::string output_dir = mesh_writer.GetOutputDirectory();

        TRI_READ_2 *p_new_mesh_reader;
        TS_ASSERT_THROWS_NOTHING(p_new_mesh_reader = 
            new TRI_READ_2(output_dir + "MeshFromFemlabViaMesh"));
                            
        delete p_new_mesh_reader;                            
    }
    
    void TestTrianglesToMeshalyzerViaMesh(void)
    {    
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        MeshalyzerMeshWriter<3,3> mesh_writer("", "MeshFromTetgenViaMesh");
    
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);
        
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
    }

    void TestTrianglesToCoolGraphicsViaMesh(void)
    {    
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        bool set_CG_format=true;            
        MeshalyzerMeshWriter<3,3> mesh_writer("CGFromTetgenViaMesh", "CGFromTetgenViaMesh", set_CG_format);

        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);
    
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
    }
    
    
    void TestTriangles1DMeshIn2DSpace()
    {        
        TrianglesMeshReader<1,2> mesh_reader(
                          "mesh/test/data/circle_outline");
        ConformingTetrahedralMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);         

        TrianglesMeshWriter<1,2> mesh_writer("","1dMeshIn2dSpace");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();
        
        TrianglesMeshReader<1,2> mesh_reader2(
                          output_dir+"1dMeshIn2dSpace");
        TS_ASSERT_EQUALS( mesh_reader2.GetNumNodes(), 100); 
        
        TS_ASSERT_EQUALS( mesh_reader2.GetNumElements(),100); 
        
        TS_ASSERT_EQUALS( mesh_reader2.GetNumFaces(), 0);
                          
        
    }         
 
    void Test2DMeshIn3DSpace()
    {
        
        TrianglesMeshReader<2,3> mesh_reader(
                          "mesh/test/data/slab_395_elements");   
                          
        ConformingTetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<2,3> mesh_writer("","2dMeshIn3dSpace");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,3> mesh_reader2(
                          output_dir+"2dMeshIn3dSpace");
        TS_ASSERT_EQUALS( mesh_reader2.GetNumNodes(), 132); 
        
        TS_ASSERT_EQUALS( mesh_reader2.GetNumElements(), 224); 
        
        TS_ASSERT_EQUALS( mesh_reader2.GetNumFaces(), 0);


    }             
   
};

#endif //_TESTMEMFEMMESHREADER_HPP_
