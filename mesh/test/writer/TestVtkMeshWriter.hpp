/*

Copyright (C) University of Oxford, 2005-2011

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


#ifndef _TESTVTKMESHWRITER_HPP_
#define _TESTVTKMESHWRITER_HPP_

#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "OutputFileHandler.hpp"
#include "TetrahedralMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include <iostream>

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkVersion.h>
#endif

class TestVtkMeshWriter : public CxxTest::TestSuite
{

public:

    void TestBasicVtkMeshWriter() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "cube_2mm_12_elements");

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        //1.6K uncompressed, 1.3K compressed
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkMeshWriter/";

        std::string target_file;
        if (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION==0)
        {
            target_file = "mesh/test/data/TestVtkMeshWriter/cube_2mm_12_elements.vtu";
        }
        else
        {
            target_file = "mesh/test/data/TestVtkMeshWriter/cube_2mm_12_elements_v52.vtu";
        }
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/cube_2mm_12_elements.vtu " + target_file).c_str()), 0);
#endif //CHASTE_VTK
    }
    
    void TestSequentialMeshCannotWriteParallelFiles() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "cube_2mm_12_elements");
        
        TS_ASSERT_THROWS_THIS( writer.SetParallelFiles(mesh), 
                               "Cannot write parallel files using a sequential mesh");
#endif //CHASTE_VTK        
    }

    void TestParallelVtkMeshWriter() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        {
            TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
            DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
            mesh.ConstructFromMeshReader(reader);
            
            VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "cube_2mm_12_elements");
            
            writer.SetParallelFiles(mesh);

            // Add distance from origin into the node "point" data
            std::vector<double> distance;
            for (DistributedTetrahedralMesh<3,3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                   node_iter != mesh.GetNodeIteratorEnd();
                   ++node_iter)
            {
                distance.push_back(norm_2(node_iter->rGetLocation()));
            }
            writer.AddPointData("Distance from origin", distance);
            
            // Add location (vector) to "point" data
            std::vector< c_vector<double, 3> > location;
            for (DistributedTetrahedralMesh<3,3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                   node_iter != mesh.GetNodeIteratorEnd();
                   ++node_iter)
            {
                location.push_back(node_iter->rGetLocation());
            }
            writer.AddPointData("Location", location);
            
            writer.WriteFilesUsingMesh(mesh);
        }

        if (PetscTools::IsSequential()) ///\todo #1494
        {        
            //1.6K uncompressed, 1.3K compressed
            std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkMeshWriter/";
    
            std::string target_file;
            if (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION==0)
            {
                target_file = "mesh/test/data/TestVtkMeshWriter/cube_2mm_12_elements_with_distance.vtu";
            }
            else
            {
                target_file = "mesh/test/data/TestVtkMeshWriter/cube_2mm_12_elements_v52.vtu";
            }
            TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/cube_2mm_12_elements.vtu " + target_file).c_str()), 0);
        }
        
        {
            TrianglesMeshReader<2,2> reader2("mesh/test/data/2D_0_to_1mm_200_elements");
            DistributedTetrahedralMesh<2,2> mesh2;
            mesh2.ConstructFromMeshReader(reader2);
    
            VtkMeshWriter<2,2> writer2("TestVtkMeshWriter", "2D_0_to_1mm_200_elements_parallel_data", false);
            writer2.SetParallelFiles(mesh2);
            // Add distance from origin into the node "point" data
            std::vector<double> rank;
            //Real rank for the owned nodes
            for (unsigned i=0; i<mesh2.GetNumLocalNodes(); i++)
            {
                rank.push_back(PetscTools::GetMyRank());
            }
            writer2.AddPointData("Process rank", rank);
    
    
            writer2.WriteFilesUsingMesh(mesh2);
        }
#endif //CHASTE_VTK
    }

    void TestVtkMeshWriter2D() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_200_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<2,2> writer("TestVtkMeshWriter", "2D_0_to_1mm_200_elements", false);

        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(mesh.GetNode(i)->rGetLocation()));
        }
        writer.AddPointData("Distance from origin", distance);

        // Add fibre type to "point" data
        std::vector< c_vector<double, 2> > location;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            location.push_back(mesh.GetNode(i)->rGetLocation());
        }
        writer.AddPointData("Location", location);

        // Add element quality into the element "cell" data
        std::vector<double> quality;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            quality.push_back(mesh.GetElement(i)->CalculateQuality());
        }
        writer.AddCellData("Quality", quality);

        // Add fibre type to "cell" data
        std::vector< c_vector<double, 2> > centroid;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            centroid.push_back(mesh.GetElement(i)->CalculateCentroid());
        }
        writer.AddCellData("Centroid", centroid);


        writer.WriteFilesUsingMesh(mesh);
        //13K uncompressed, 3.7K compressed
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkMeshWriter/";

        std::string target_file;
        if (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION==0)
        {
            target_file = "mesh/test/data/TestVtkMeshWriter/2D_0_to_1mm_200_elements.vtu";
        }
        else
        {
            target_file = "mesh/test/data/TestVtkMeshWriter/2D_0_to_1mm_200_elements_v52.vtu";
        }
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/2D_0_to_1mm_200_elements.vtu " + target_file).c_str()), 0);
#endif //CHASTE_VTK
    }


    void TestVtkMeshWriterWithData() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        TrianglesMeshReader<3,3> reader("heart/test/data/HeartDecimation_173nodes");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "heart_decimation", false);

        // Add element quality into the element "cell" data
        std::vector<double> quality;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            quality.push_back(mesh.GetElement(i)->CalculateQuality());
        }
        writer.AddCellData("Quality", quality);

        // Add fibre type to "cell" data
        std::vector< c_vector<double, 3> > centroid;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            centroid.push_back(mesh.GetElement(i)->CalculateCentroid());
        }
        writer.AddCellData("Centroid", centroid);

        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(mesh.GetNode(i)->rGetLocation()));
        }
        writer.AddPointData("Distance from origin", distance);

        // Add fibre type to "point" data
        std::vector< c_vector<double, 3> > location;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            location.push_back(mesh.GetNode(i)->rGetLocation());
        }
        writer.AddPointData("Location", location);

        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));


        //32K uncompressed, 19K compressed
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestVtkMeshWriter/";

        std::string target_file;
        if (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION==0)
        {
            target_file = "mesh/test/data/TestVtkMeshWriter/heart_decimation.vtu";
        }
        else
        {
            target_file = "mesh/test/data/TestVtkMeshWriter/heart_decimation_v52.vtu";
        }
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/heart_decimation.vtu " + target_file).c_str()), 0);
#endif //CHASTE_VTK
    }
    void TestWriteNodesWithoutMesh()
    {
 #ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);
        
        //Note in my version of Paraview, you need data on points before you
        //can view with Glyphs.
        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "just_nodes", false);



        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(mesh.GetNode(i)->rGetLocation()));
        }
        writer.AddPointData("Distance from origin", distance);

        // Add boundary node "point" data
        std::vector<double> boundary;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            boundary.push_back(mesh.GetNode(i)->IsBoundaryNode());
        }
        writer.AddPointData("Boundary", boundary);
        
        // Add fibre type to "point" data
        std::vector< c_vector<double, 3> > location;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            location.push_back(mesh.GetNode(i)->rGetLocation());
        }
        writer.AddPointData("Location", location);

        writer.WriteFilesUsingMesh(mesh);

        
        //When the mesh goes out of scope, then it's a different set of nodes that get destroyed
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
#endif //CHASTE_VTK
    }
};

#endif //_TESTVTKMESHWRITER_HPP_
