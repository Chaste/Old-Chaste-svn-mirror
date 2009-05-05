/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef _TESTMESHWRITERS_HPP_
#define _TESTMESHWRITERS_HPP_

#include <cxxtest/TestSuite.h>
#include "MemfemMeshReader.hpp"
#include "FemlabMeshReader.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "OutputFileHandler.hpp"
#include "MutableMesh.hpp"
#include "CmguiWriter.hpp"
#include "VtkWriter.hpp"


class TestMeshWriters : public CxxTest::TestSuite
{
public:

    void TestMemfemtoTetgen()
    {
        TrianglesMeshWriter<3,3> mesh_writer("", "MeshFromMemfem");
        MemfemMeshReader<3,3> import_mesh_reader("mesh/test/data/Memfem_slab");

        mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader);
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<3,3>* p_new_mesh_reader;
        p_new_mesh_reader = new TrianglesMeshReader<3,3>(output_dir + "MeshFromMemfem");

        delete p_new_mesh_reader;
    }

    void TestFemlabtoTriangles()
    {
        TrianglesMeshWriter<2,2> mesh_writer("", "MeshFromFemlab");

        FemlabMeshReader<2,2> import_mesh_reader("mesh/test/data/",
                                                 "femlab_lshape_nodes.dat",
                                                 "femlab_lshape_elements.dat",
                                                 "femlab_lshape_edges.dat");

        mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader);
        std::string output_dir = mesh_writer.GetOutputDirectory();

        TrianglesMeshReader<2,2>* p_new_mesh_reader;
        p_new_mesh_reader = new TrianglesMeshReader<2,2>(output_dir + "MeshFromFemlab");

        delete p_new_mesh_reader;
    }

    void TestTrianglesToMeshalyzer1d()
    {
        TrianglesMeshReader<1,1> import_mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MeshalyzerMeshWriter<1,1> mesh_writer("", "MeshFromTetgen");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }

    void TestTrianglesToMeshalyzer2d()
    {
        TrianglesMeshReader<2,2> import_mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        MeshalyzerMeshWriter<2,2> mesh_writer("", "MeshFromTetgen");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }

    void TestTrianglesToMeshalyzer3d()
    {
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        MeshalyzerMeshWriter<3,3> mesh_writer("", "MeshFromTetgen");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }
    
    void TestTrianglesToMeshalyzer1din3d()
    {
        TrianglesMeshReader<1,3> import_mesh_reader("mesh/test/data/trivial_1d_in_3d_mesh");
        MeshalyzerMeshWriter<1,3> mesh_writer("", "Mesh");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }    

    void TestTrianglesToCoolGraphics()
    {
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        bool set_CG_format = true;

        MeshalyzerMeshWriter<3,3> mesh_writer("CGFromTetgen", "CGFromTetgen", true, set_CG_format);

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader));
    }

    void TestFemlabtoTrianglesViaMesh()
    {
        TrianglesMeshWriter<2,2> mesh_writer("", "MeshFromFemlabViaMesh");

        FemlabMeshReader<2,2> import_mesh_reader("mesh/test/data/",
                                                 "femlab_lshape_nodes.dat",
                                                 "femlab_lshape_elements.dat",
                                                 "femlab_lshape_edges.dat");

        TetrahedralMesh<2,2> mesh;
        bool cull_internal_faces = true;
        mesh.ConstructFromMeshReader(import_mesh_reader, cull_internal_faces);

        mesh_writer.WriteFilesUsingMesh(mesh);
        std::string output_dir = mesh_writer.GetOutputDirectory();

        TrianglesMeshReader<2,2>* p_new_mesh_reader;
        p_new_mesh_reader = new TrianglesMeshReader<2,2>(output_dir + "MeshFromFemlabViaMesh");

        delete p_new_mesh_reader;
    }

    void TestTrianglesToMeshalyzerViaMesh1d()
    {
        TrianglesMeshReader<1,1> import_mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MeshalyzerMeshWriter<1,1> mesh_writer("", "MeshFromTetgenViaMesh");

        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
    }

    void TestTrianglesToMeshalyzerViaMesh2d()
    {
        TrianglesMeshReader<2,2> import_mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        MeshalyzerMeshWriter<2,2> mesh_writer("", "MeshFromTetgenViaMesh");

        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
    }

    void TestTrianglesToMeshalyzerViaMesh3d()
    {
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        MeshalyzerMeshWriter<3,3> mesh_writer("", "MeshFromTetgenViaMesh");

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
    }

    void TestTrianglesToCoolGraphicsViaMesh()
    {
        TrianglesMeshReader<3,3> import_mesh_reader("mesh/test/data/slab_138_elements");
        bool set_CG_format = true;
        MeshalyzerMeshWriter<3,3> mesh_writer("CGFromTetgenViaMesh", "CGFromTetgenViaMesh", true, set_CG_format);

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(import_mesh_reader);

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));
    }

    void TestTriangles1DClosedMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/circle_outline");
        TetrahedralMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<1,2> mesh_writer("", "1dClosedMeshIn2dSpace");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();

        std::cout << output_dir + "1dClosedMeshIn2dSpace" << std::endl;

        TrianglesMeshReader<1,2> mesh_reader2(output_dir + "1dClosedMeshIn2dSpace");

        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 100u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 100u);

        // By default all nodes are read as potential "faces"
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 100u);
    }

    void TestTriangles1DMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/semicircle_outline");
        TetrahedralMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<1,2> mesh_writer("", "1dMeshIn2dSpace");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();

        TrianglesMeshReader<1,2> mesh_reader2(output_dir + "1dMeshIn2dSpace");

        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 51u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 50u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 51u);
    }

    void TestTriangles1DMeshIn2DSpaceWithDeletedNode()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/semicircle_outline");
        MutableMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.DeleteBoundaryNodeAt(0);

        TrianglesMeshWriter<1,2> mesh_writer("", "1dMeshIn2dSpaceWithDeletedNode");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();

        TrianglesMeshReader<1,2> mesh_reader2(output_dir + "1dMeshIn2dSpaceWithDeletedNode");
        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 50u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 49u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 50u);
    }

    void Test2DClosedMeshIn3DSpace()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/slab_395_elements");

        TetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<2,3> mesh_writer("", "2dClosedMeshIn3dSpace");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,3> mesh_reader2(output_dir+"2dClosedMeshIn3dSpace");

        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 132u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 224u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 0u);
    }

    void Test2DMeshIn3DSpace()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");

        TetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<2,3> mesh_writer("", "2dMeshIn3dSpace");

        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,3> mesh_reader2(output_dir + "2dMeshIn3dSpace");

        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 312u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 522u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), mesh_reader.GetNumNodes());
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), mesh_reader.GetNumElements());

        // Note that this not a straight conversion, since we have culled internal data
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 100u); // culling now occurs in the reader
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 100u);
    }

    void TestCmguiWriter() throw(Exception)
    {
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        CmguiWriter writer("TestCmguiWriter", "cube_2mm_12_elements");

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/cube_2mm_12_elements.exnode mesh/test/data/TestCmguiWriter/cube_2mm_12_elements.exnode").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/cube_2mm_12_elements.exelem mesh/test/data/TestCmguiWriter/cube_2mm_12_elements.exelem").c_str()), 0);
    }

    void TestVtkWriter() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkWriter writer("TestVtkWriter", "cube_2mm_12_elements");

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/cube_2mm_12_elements.vtu mesh/test/data/TestVtkWriter/cube_2mm_12_elements.vtu").c_str()), 0);
#endif //CHASTE_VTK
    }

    void TestVtkWriterWithData() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<3,3> reader("heart/test/data/HeartDecimation_173nodes");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkWriter writer("TestVtkWriter", "heart_decimation", false);

        // Add element quality into the element "cell" data
        std::vector<double> quality;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            quality.push_back(mesh.GetElement(i)->CalculateQuality());
        }
        writer.AddCellData("Quality", quality);

        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(mesh.GetNode(i)->rGetLocation()));
        }
        writer.AddPointData("Distance from origin", distance);

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkWriter/";

        // Only compare the first 531 bytes for now (the offsets and stuff seem to be changing)
        TS_ASSERT_EQUALS(system(("cmp -n 531 " + results_dir + "/heart_decimation.vtu mesh/test/data/TestVtkWriter/heart_decimation.vtu").c_str()), 0);
#endif //CHASTE_VTK
    }

    void TestWriteFilesUsingMeshReaderPermuted()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");

        unsigned num_nodes = mesh_reader.GetNumNodes();

        std::vector<unsigned> permutation;

        //Cover the zero permutation case
        TrianglesMeshWriter<1,1> mesh_writer_straight("", "MeshReaderNotPermuted");
        TS_ASSERT_EQUALS(permutation.size(), 0U);
        mesh_writer_straight.WriteFilesUsingMeshReader(mesh_reader, permutation);
        mesh_reader.Reset();
       
        std::string filename = "MeshReaderPermuted";

        for (unsigned index=0; index<num_nodes; index++)
        {
            permutation.push_back(num_nodes - index - 1);
        }
        TrianglesMeshWriter<1,1> mesh_writer("", filename);
        mesh_writer.WriteFilesUsingMeshReader(mesh_reader, permutation);
 
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<1,1> permuted_mesh_reader(output_dir + filename);

        // This mesh is made of 11 points located spaced 0.1cm each
        mesh_reader.Reset();
        for (unsigned index=0; index<num_nodes; index++)
        {
            std::vector<double> original_node = mesh_reader.GetNextNode();
            TS_ASSERT_DELTA(original_node[0], 0.1*index, 1e-6);

            std::vector<double> permuted_node = permuted_mesh_reader.GetNextNode();
            TS_ASSERT_DELTA(permuted_node[0], 1.0 - 0.1*index, 1e-6);
        }

        for (unsigned elem_index=0; elem_index<mesh_reader.GetNumElements(); elem_index++)
        {
            ElementData original_element = mesh_reader.GetNextElementData();
            ElementData permuted_element = permuted_mesh_reader.GetNextElementData();

            for (unsigned local_node_index=0; local_node_index<original_element.NodeIndices.size(); local_node_index++)
            {
                unsigned original_global_node_index = original_element.NodeIndices[local_node_index];
                TS_ASSERT_EQUALS(permuted_element.NodeIndices[local_node_index],
                                 permutation[original_global_node_index]);
            }
        }
    }
};

#endif //_TESTMEMFEMMESHREADER_HPP_
