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
    
    
    
    
    void TestTrianglesToMeshalyzer1d(void)
    {
        TrianglesMeshReader<1,1> import_mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MeshalyzerMeshWriter<1,1> mesh_writer("", "MeshFromTetgen");
        
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }
    
    void TestTrianglesToMeshalyzer2d(void)
    {
        TrianglesMeshReader<2,2> import_mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        MeshalyzerMeshWriter<2,2> mesh_writer("", "MeshFromTetgen");
        
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }
    
    void TestTrianglesToMeshalyzer3d(void)
    {
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        MeshalyzerMeshWriter<3,3> mesh_writer("", "MeshFromTetgen");
        
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }
    
    void TestTrianglesToCoolGraphics(void)
    {
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        bool set_CG_format=true;
        
        MeshalyzerMeshWriter<3,3> mesh_writer("CGFromTetgen", "CGFromTetgen", true, set_CG_format);
        
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
    
    
    void TestTrianglesToMeshalyzerViaMesh1d(void)
    {
        TrianglesMeshReader<1,1> import_mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MeshalyzerMeshWriter<1,1> mesh_writer("", "MeshFromTetgenViaMesh");
        
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);
        
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
    }
    
    void TestTrianglesToMeshalyzerViaMesh2d(void)
    {
        TrianglesMeshReader<2,2> import_mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        MeshalyzerMeshWriter<2,2> mesh_writer("", "MeshFromTetgenViaMesh");
        
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);
        
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
    }
    
    void TestTrianglesToMeshalyzerViaMesh3d(void)
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
        MeshalyzerMeshWriter<3,3> mesh_writer("CGFromTetgenViaMesh", "CGFromTetgenViaMesh", true, set_CG_format);
        
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);
        
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
    }
    
    
    void TestTriangles1DClosedMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader(
            "mesh/test/data/circle_outline");
        ConformingTetrahedralMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TrianglesMeshWriter<1,2> mesh_writer("","1dClosedMeshIn2dSpace");
        
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
        
        std::string output_dir = mesh_writer.GetOutputDirectory();
        
        TrianglesMeshReader<1,2> mesh_reader2(
            output_dir+"1dClosedMeshIn2dSpace");
        TS_ASSERT_EQUALS( mesh_reader2.GetNumNodes(), 100U);
        
        TS_ASSERT_EQUALS( mesh_reader2.GetNumElements(),100U);
        
        //By default all nodes are read as potential "faces"
        TS_ASSERT_EQUALS( mesh_reader2.GetNumFaces(), 100U);
    }
    
    void TestTriangles1DMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader(
            "mesh/test/data/semicircle_outline");
        ConformingTetrahedralMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TrianglesMeshWriter<1,2> mesh_writer("","1dMeshIn2dSpace");
        
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
        
        std::string output_dir = mesh_writer.GetOutputDirectory();
        
        TrianglesMeshReader<1,2> mesh_reader2(
            output_dir+"1dMeshIn2dSpace");
        TS_ASSERT_EQUALS( mesh_reader2.GetNumNodes(), 51U);
        
        TS_ASSERT_EQUALS( mesh_reader2.GetNumElements(), 50U);
        
        TS_ASSERT_EQUALS( mesh_reader2.GetNumFaces(), 51U);
        
        
    }
    
    void Test2DClosedMeshIn3DSpace()
    {
    
        TrianglesMeshReader<2,3> mesh_reader(
            "mesh/test/data/slab_395_elements");
            
        ConformingTetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TrianglesMeshWriter<2,3> mesh_writer("","2dClosedMeshIn3dSpace");
        
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
        
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,3> mesh_reader2(
            output_dir+"2dClosedMeshIn3dSpace");
        TS_ASSERT_EQUALS( mesh_reader2.GetNumNodes(), 132U);
        
        TS_ASSERT_EQUALS( mesh_reader2.GetNumElements(), 224U);
        
        TS_ASSERT_EQUALS( mesh_reader2.GetNumFaces(), 0U);
        
        
    }
    
    void Test2DMeshIn3DSpace()
    {
    
        TrianglesMeshReader<2,3> mesh_reader(
            "mesh/test/data/disk_in_3d");
            
        ConformingTetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TrianglesMeshWriter<2,3> mesh_writer("","2dMeshIn3dSpace");
        
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
        
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,3> mesh_reader2(
            output_dir+"2dMeshIn3dSpace");
        TS_ASSERT_EQUALS( mesh_reader2.GetNumNodes(), 312U);
        TS_ASSERT_EQUALS( mesh_reader2.GetNumElements(), 522U);
        TS_ASSERT_EQUALS( mesh_reader2.GetNumNodes(), mesh_reader.GetNumNodes());
        TS_ASSERT_EQUALS( mesh_reader2.GetNumElements(), mesh_reader.GetNumElements());
        
        //Note that this not a straight conversion, since we have culled internal data
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 833U);
        TS_ASSERT_EQUALS( mesh_reader2.GetNumFaces(), 100U);
        
    }
    
};

#endif //_TESTMEMFEMMESHREADER_HPP_
