/*

Copyright (C) University of Oxford, 2005-2010

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
#include "CmguiMeshWriter.hpp"
#include "CmguiDeformedSolutionsWriter.hpp"
#include "VtkWriter.hpp"
#include "QuadraticMesh.hpp"

class TestMeshWriters : public CxxTest::TestSuite
{
public:

    void TestMemfemToTetgen()
    {
        TrianglesMeshWriter<3,3> mesh_writer("", "MeshFromMemfem");
        MemfemMeshReader<3,3> import_mesh_reader("mesh/test/data/Memfem_slab");

        mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader);
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<3,3>* p_new_mesh_reader;
        p_new_mesh_reader = new TrianglesMeshReader<3,3>(output_dir + "MeshFromMemfem");

        delete p_new_mesh_reader;
    }

    void TestFemlabToTriangles()
    {
        TrianglesMeshWriter<2,2> mesh_writer("", "MeshFromFemlab");

        FemlabMeshReader<2,2> import_mesh_reader("mesh/test/data/",
                                                 "femlab_lshape_nodes.dat",
                                                 "femlab_lshape_elements.dat",
                                                 "femlab_lshape_edges.dat");

        //TS_ASSERT_EQUALS(import_mesh_reader.GetNumFaces(), 54U); //Has internal faces
        TS_ASSERT_EQUALS(import_mesh_reader.GetNumFaces(), 40U); //Has no internal faces
        mesh_writer.WriteFilesUsingMeshReader(import_mesh_reader);
        std::string output_dir = mesh_writer.GetOutputDirectory();

        TrianglesMeshReader<2,2>* p_new_mesh_reader;
        p_new_mesh_reader = new TrianglesMeshReader<2,2>(output_dir + "MeshFromFemlab");
        //TS_ASSERT_EQUALS(p_new_mesh_reader->GetNumFaces(), 54U); //No faces have been culled
        TS_ASSERT_EQUALS(p_new_mesh_reader->GetNumFaces(), 40U); //No faces have been culled

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
    
    void TestTrianglesToMeshalyzer1dIn3d()
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

        TS_ASSERT_EQUALS(import_mesh_reader.GetNumFaces(), 40U);//Has no internal faces
        TetrahedralMesh<2,2> mesh;

        mesh.ConstructFromMeshReader(import_mesh_reader);

        mesh_writer.WriteFilesUsingMesh(mesh);
        std::string output_dir = mesh_writer.GetOutputDirectory();

        TrianglesMeshReader<2,2>* p_new_mesh_reader;
        p_new_mesh_reader = new TrianglesMeshReader<2,2>(output_dir + "MeshFromFemlabViaMesh");
        TS_ASSERT_EQUALS(p_new_mesh_reader->GetNumFaces(), 40U); //Internal faces have been culled
        
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
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 0u);
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
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 2u);
    }

    void TestTriangles1DMeshIn2DSpaceWithDeletedNode() throw (Exception)
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/semicircle_outline");
        MutableMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.DeleteBoundaryNodeAt(0);

        TrianglesMeshWriter<1,2> mesh_writer("", "1dMeshIn2dSpaceWithDeletedNode");
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        TrianglesMeshWriter<1,2> mesh_writer2("", "1dMeshIn2dSpaceWithDeletedNodeConst");
        //mesh_writer2.WriteFilesUsingMesh(static_cast<const MutableMesh<1,2> &>(mesh));
        mesh_writer2.WriteFilesUsingMesh(mesh);

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<1,2> mesh_reader2(output_dir + "1dMeshIn2dSpaceWithDeletedNode");
        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 50u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 49u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 2u);
        
        TrianglesMeshReader<1,2> mesh_reader3(output_dir + "1dMeshIn2dSpaceWithDeletedNodeConst");
        TS_ASSERT_EQUALS(mesh_reader3.GetNumNodes(), 50u);
        TS_ASSERT_EQUALS(mesh_reader3.GetNumElements(), 49u);
        TS_ASSERT_EQUALS(mesh_reader3.GetNumFaces(), 2u);
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
    
    void TestQuadratic1D() throw (Exception)
    {
        QuadraticMesh<1> mesh;
        TrianglesMeshReader<1,1> mesh_reader1("mesh/test/data/1D_0_to_1_10_elements_quadratic", 2, 1, false);
        mesh.ConstructFromMeshReader(mesh_reader1);
        TrianglesMeshWriter<1,1> mesh_writer("", "1d_quadratic");
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<1,1> mesh_reader(output_dir + "1d_quadratic", 2, 2);
        
        //Test that reader is reading correctly
        TS_ASSERT_EQUALS(mesh_reader.GetNextNode().size(), 1u);
        TS_ASSERT_EQUALS(mesh_reader.GetNextElementData().NodeIndices.size(), 3u);
        TS_ASSERT_EQUALS(mesh_reader.GetNextFaceData().NodeIndices.size(), 1u);
        
        //Test that mesh can be reconstructed
        QuadraticMesh<1> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh2.GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(mesh2.GetBoundaryElement(0)->GetNumNodes(), 1U);
        TS_ASSERT_EQUALS(mesh2.GetBoundaryElement(0)->GetNodeGlobalIndex(0), mesh.GetBoundaryElement(0)->GetNodeGlobalIndex(0));
    }

    void TestQuadratic2D() throw (Exception)
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader1("mesh/test/data/square_128_elements_fully_quadratic", 2, 2, false);
        mesh.ConstructFromMeshReader(mesh_reader1);
        TrianglesMeshWriter<2,2> mesh_writer("", "2d_quadratic");
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,2> mesh_reader(output_dir + "2d_quadratic", 2, 2);

        //Test that reader is reading correctly
        TS_ASSERT_EQUALS(mesh_reader.GetNextNode().size(), 2u);
        TS_ASSERT_EQUALS(mesh_reader.GetNextElementData().NodeIndices.size(), 6u);
        TS_ASSERT_EQUALS(mesh_reader.GetNextFaceData().NodeIndices.size(), 3u);
        
        //Test that mesh can be reconstructed
        QuadraticMesh<2> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh2.GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(mesh2.GetBoundaryElement(0)->GetNumNodes(), 3U);
        TS_ASSERT_EQUALS(mesh2.GetBoundaryElement(0)->GetNodeGlobalIndex(2), mesh.GetBoundaryElement(0)->GetNodeGlobalIndex(2));
    }

    void TestQuadratic3D() throw (Exception)
    {
        QuadraticMesh<3> mesh;
        TrianglesMeshReader<3,3> mesh_reader1("mesh/test/data/3D_Single_tetrahedron_element_quadratic", 2, 1, false);
        mesh.ConstructFromMeshReader(mesh_reader1);
        TrianglesMeshWriter<3,3> mesh_writer("", "3d_quadratic");
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(mesh));

        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<3,3> mesh_reader(output_dir + "3d_quadratic", 2, 2);

        //Test that reader is reading correctly
        TS_ASSERT_EQUALS(mesh_reader.GetNextNode().size(), 3u);
        TS_ASSERT_EQUALS(mesh_reader.GetNextElementData().NodeIndices.size(), 10u);
        TS_ASSERT_EQUALS(mesh_reader.GetNextFaceData().NodeIndices.size(), 6u);
        
        //Test that mesh can be reconstructed
        QuadraticMesh<3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh2.GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(mesh2.GetBoundaryElement(0)->GetNumNodes(), 6U);
        TS_ASSERT_EQUALS(mesh2.GetBoundaryElement(0)->GetNodeGlobalIndex(5), mesh.GetBoundaryElement(0)->GetNodeGlobalIndex(5));
    }
    
    void TestCmguiMeshWriter3D() throw(Exception)
    {
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        CmguiMeshWriter<3,3> writer("TestCmguiMeshWriter3D", "cube_2mm_12_elements");

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiMeshWriter3D/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/cube_2mm_12_elements.exnode mesh/test/data/TestCmguiMeshWriter/cube_2mm_12_elements.exnode").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/cube_2mm_12_elements.exelem mesh/test/data/TestCmguiMeshWriter/cube_2mm_12_elements.exelem").c_str()), 0);
        
        //now test the set method for additional fields. We set two fields.
        CmguiMeshWriter<3,3> writer2("TestCmguiMeshWriterAdditionalHeaders3D", "cube_2mm_12_elements");
        
        std::vector<std::string> field_names;
        field_names.push_back("V");
        field_names.push_back("Phi_e");
        std::string results_dir2 = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiMeshWriterAdditionalHeaders3D";
        writer2.SetAdditionalFieldNames(field_names);
        TS_ASSERT_THROWS_NOTHING(writer2.WriteFilesUsingMesh(mesh));
        TS_ASSERT_EQUALS(system(("cmp " + results_dir2 + "/cube_2mm_12_elements.exelem mesh/test/data/TestCmguiMeshWriter/cube_2mm_12_elements_additional_fields.exelem").c_str()), 0);
    }
    
    void TestCmguiMeshWriter2D() throw(Exception)
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        CmguiMeshWriter<2,2> writer("TestCmguiMeshWriter2D", "square_128_elements");

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiMeshWriter2D/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/square_128_elements.exnode mesh/test/data/TestCmguiMeshWriter/square_128_elements.exnode").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/square_128_elements.exelem mesh/test/data/TestCmguiMeshWriter/square_128_elements.exelem").c_str()), 0);
        
        //now test the set method for additional fields. We set two fields.
        CmguiMeshWriter<2,2> writer2("TestCmguiMeshWriterAdditionalHeaders2D", "square_128_elements");
        
        std::vector<std::string> field_names;
        field_names.push_back("V");
        field_names.push_back("Phi_e");
        std::string results_dir2 = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiMeshWriterAdditionalHeaders2D";
        writer2.SetAdditionalFieldNames(field_names);
        TS_ASSERT_THROWS_NOTHING(writer2.WriteFilesUsingMesh(mesh));
        TS_ASSERT_EQUALS(system(("cmp " + results_dir2 + "/square_128_elements.exelem mesh/test/data/TestCmguiMeshWriter/square_128_elements_additional_fields.exelem").c_str()), 0);
    }

    void TestCmguiMeshWriter1D() throw(Exception)
    {
        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        CmguiMeshWriter<1,1> writer("TestCmguiMeshWriter1D", "1D_0_to_1_100_elements");

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiMeshWriter1D/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/1D_0_to_1_100_elements.exnode mesh/test/data/TestCmguiMeshWriter/1D_0_to_1_100_elements.exnode").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/1D_0_to_1_100_elements.exelem mesh/test/data/TestCmguiMeshWriter/1D_0_to_1_100_elements.exelem").c_str()), 0);
        
        //now test the set method for additional fields. We set two fields.
        CmguiMeshWriter<1,1> writer2("TestCmguiMeshWriterAdditionalHeaders1D", "1D_0_to_1_100_elements");
        
        std::vector<std::string> field_names;
        field_names.push_back("V");
        field_names.push_back("Phi_e");
        std::string results_dir2 = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiMeshWriterAdditionalHeaders1D";
        writer2.SetAdditionalFieldNames(field_names);
        TS_ASSERT_THROWS_NOTHING(writer2.WriteFilesUsingMesh(mesh));
        TS_ASSERT_EQUALS(system(("cmp " + results_dir2 + "/1D_0_to_1_100_elements.exelem mesh/test/data/TestCmguiMeshWriter/1D_0_to_1_100_elements_additional_fields.exelem").c_str()), 0);
    }
    
    void TestCmguiDeformedSolutionsWriter() throw(Exception)
    {
        QuadraticMesh<2> mesh(1.0, 2.0, 2, 2);
       
        CmguiDeformedSolutionsWriter<2> writer("TestCmguiDeformedSolutionsWriter", "solution", mesh);
        
        // writes solution_0.exnode and solution_1.node using the mesh
        writer.WriteInitialMesh();
        
        // set up a deformed positions vector
        std::vector<c_vector<double,2> > deformed_positions(mesh.GetNumNodes(), zero_vector<double>(2));
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            deformed_positions[i](0) = 2*mesh.GetNode(i)->rGetLocation()[0];
            deformed_positions[i](1) = 3*mesh.GetNode(i)->rGetLocation()[1];
        }
        // write solution_1.exnode
        writer.WriteDeformationPositions(deformed_positions, 1);

        // ..and again
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            deformed_positions[i](0) = 4*mesh.GetNode(i)->rGetLocation()[0];
            deformed_positions[i](1) = 6*mesh.GetNode(i)->rGetLocation()[1];
        }
        // write solution_1.exnode
        writer.WriteDeformationPositions(deformed_positions, 2);
        // write LoadSolutions.com
        writer.WriteCmguiScript();
        
        deformed_positions.push_back(zero_vector<double>(2));
        
        TS_ASSERT_THROWS_CONTAINS(writer.WriteDeformationPositions(deformed_positions, 3), "The size of rDeformedPositions does not match the number of nodes in the mesh");
        
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiDeformedSolutionsWriter";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/solution_0.exelem mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_0.exelem").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/solution_0.exnode mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_0.exnode").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/solution_1.exnode mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_1.exnode").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/solution_2.exnode mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_2.exnode").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/LoadSolutions.com mesh/test/data/TestCmguiDeformedSolutionsWriter/LoadSolutions.com").c_str()), 0);
    }
    
    void TestCmguiDeformedSolutionsWriterConvertOutput() throw(Exception)
    {
        QuadraticMesh<2> mesh(1.0, 2.0, 2, 2);
        CmguiDeformedSolutionsWriter<2> writer("TestCmguiDeformedSolutionsWriter_ConvertOutput", "solution", mesh);

        // throws as file doesn't exist        
        TS_ASSERT_THROWS_CONTAINS(writer.ConvertOutput("blahblahblah", "blah", 1), "Could not open file:");

        // convert some files
        writer.ConvertOutput("mesh/test/data/TestCmguiDeformedSolutionsWriter", "myoldsolution", 2);

        // myoldsolution_1.nodes and myoldsolution_2.nodes are in chaste format but match the two deformed_positions written in the
        // test above, so the new output should match the good output from the test above.

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestCmguiDeformedSolutionsWriter_ConvertOutput";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/solution_0.exelem mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_0.exelem").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/solution_0.exnode mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_0.exnode").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/solution_1.exnode mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_1.exnode").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/solution_2.exnode mesh/test/data/TestCmguiDeformedSolutionsWriter/solution_2.exnode").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/LoadSolutions.com mesh/test/data/TestCmguiDeformedSolutionsWriter/LoadSolutions.com").c_str()), 0);

        // throws as is incomplete    
        TS_ASSERT_THROWS_CONTAINS(writer.ConvertOutput("mesh/test/data/TestCmguiDeformedSolutionsWriter", "bad_myoldsolution", 1), "Error occurred when reading file");
    } 
    
    void TestVtkWriter() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkWriter<3,3> writer("TestVtkWriter", "cube_2mm_12_elements");

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        //1.6K uncompressed, 1.3K compressed
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/cube_2mm_12_elements.vtu mesh/test/data/TestVtkWriter/cube_2mm_12_elements.vtu").c_str()), 0);
#endif //CHASTE_VTK
    }

    void TestVtkWriter2D() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_200_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkWriter<2,2> writer("TestVtkWriter", "2D_0_to_1mm_200_elements", false);

        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(mesh.GetNode(i)->rGetLocation()));
        }
        writer.AddPointData("Distance from origin", distance);
        
        // Add element quality into the element "cell" data
        std::vector<double> quality;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            quality.push_back(mesh.GetElement(i)->CalculateQuality());
        }
        writer.AddCellData("Quality", quality);


        writer.WriteFilesUsingMesh(mesh);
        //13K uncompressed, 3.7K compressed
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/2D_0_to_1mm_200_elements.vtu mesh/test/data/TestVtkWriter/2D_0_to_1mm_200_elements.vtu").c_str()), 0);
#endif //CHASTE_VTK
    }


    void TestVtkWriterWithData() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<3,3> reader("heart/test/data/HeartDecimation_173nodes");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkWriter<3,3> writer("TestVtkWriter", "heart_decimation", false);

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

        
        //32K uncompressed, 19K compressed
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/heart_decimation.vtu mesh/test/data/TestVtkWriter/heart_decimation.vtu").c_str()), 0);
#endif //CHASTE_VTK
    }

    /**
     * This test is based on TestTrianglesMeshReader.hpp TestReadingElementAttributes.
     */
    void TestWritingElementAttributesInTrianglesFormat() throw (Exception)
    {
        std::string source_mesh = "mesh/test/data/1D_0_to_1_10_elements_with_attributes";
        std::string output_dir = "element_attrs";
        std::string file_from_reader = "from_reader";
        std::string file_from_mesh = "from_mesh";
        
        // Firstly, write directly from a mesh reader
        {
            TrianglesMeshReader<1,1> mesh_reader(source_mesh);
            TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 10u);
            TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);
            
            TrianglesMeshWriter<1,1> mesh_writer(output_dir, file_from_reader, true);
            mesh_writer.WriteFilesUsingMeshReader(mesh_reader);
        }
        
        // Next, write using a mesh object
        {
            TrianglesMeshReader<1,1> mesh_reader(source_mesh);
            TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 10u);
            TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);
            
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            
            TrianglesMeshWriter<1,1> mesh_writer(output_dir, file_from_mesh, false);
            mesh_writer.WriteFilesUsingMesh(mesh);
        }
        
        // Now check the written meshes
        OutputFileHandler handler(output_dir, false);
        TrianglesMeshReader<1,1> reader1(handler.GetOutputDirectoryFullPath() + file_from_reader);
        TS_ASSERT_EQUALS(reader1.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(reader1.GetNumElementAttributes(), 1u);
        TrianglesMeshReader<1,1> reader2(handler.GetOutputDirectoryFullPath() + file_from_mesh);
        TS_ASSERT_EQUALS(reader2.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(reader2.GetNumElementAttributes(), 1u);
        for (unsigned i=0; i<10; i++)
        {
            ElementData next_element_info = reader1.GetNextElementData();
            std::vector<unsigned> nodes = next_element_info.NodeIndices;
            TS_ASSERT_EQUALS(nodes.size(), 2u);
            TS_ASSERT_EQUALS(next_element_info.AttributeValue, i%5 + 1);
            
            next_element_info = reader2.GetNextElementData();
            nodes = next_element_info.NodeIndices;
            TS_ASSERT_EQUALS(nodes.size(), 2u);
            TS_ASSERT_EQUALS(next_element_info.AttributeValue, i%5 + 1);
        }
    }
    void TestWritingBinaryFormat()
    {
        /*Read as ascii*/
        TrianglesMeshReader<3,3> reader("mesh/test/data/simple_cube");

        TrianglesMeshWriter<3,3> writer_from_reader("TestMeshWriter", "simple_cube_binary_from_reader", false);
        writer_from_reader.SetWriteFilesAsBinary();
        writer_from_reader.WriteFilesUsingMeshReader(reader);
        
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);
        TrianglesMeshWriter<3,3> writer_from_mesh("TestMeshWriter", "simple_cube_binary_from_mesh", false);
        writer_from_mesh.SetWriteFilesAsBinary();
        writer_from_mesh.WriteFilesUsingMesh(mesh);
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestMeshWriter/";
        
        //This test is too stringent at the moment, but we'll do something better later
        // -n options ensure that we don't get to the provenance data
        TS_ASSERT_EQUALS(system(("cmp -n 264 " + results_dir + "/simple_cube_binary_from_reader.node mesh/test/data/simple_cube_binary.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp -n 256 " + results_dir + "/simple_cube_binary_from_reader.ele mesh/test/data/simple_cube_binary.ele").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp -n 230 " + results_dir + "/simple_cube_binary_from_reader.face mesh/test/data/simple_cube_binary.face").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp -n 264 " + results_dir + "/simple_cube_binary_from_mesh.node mesh/test/data/simple_cube_binary.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp -n 256 " + results_dir + "/simple_cube_binary_from_mesh.ele mesh/test/data/simple_cube_binary.ele").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp -n 230 " + results_dir + "/simple_cube_binary_from_mesh.face mesh/test/data/simple_cube_binary.face").c_str()), 0);
    }
 
};

#endif //_TESTMEMFEMMESHREADER_HPP_
