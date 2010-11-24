/*

Copyright (C) Fujitsu Laboratories of Europe, 2009

*/

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

#ifndef TESTVTKMESHREADER_
#define TESTVTKMESHREADER_

#include <cxxtest/TestSuite.h>
#include <fstream>

#ifdef CHASTE_VTK
//Requires  "sudo aptitude install libvtk5-dev" or similar
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#endif //CHASTE_VTK

#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "VtkMeshReader.hpp"
#include "GenericMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "PetscSetupAndFinalize.hpp"

#ifdef CHASTE_VTK
typedef VtkMeshReader<3,3> MESH_READER3;
#endif //CHASTE_VTK

class TestVtkMeshReader : public CxxTest::TestSuite
{
//Requires  "sudo aptitude install libvtk5-dev" or similar
public:

    /**
     * Check that input files are opened correctly and non-existent input files throw an Exception.
     */
    void TestFilesOpen(void) throw(Exception)
    {
#ifdef CHASTE_VTK
        TS_ASSERT_THROWS_NOTHING( MESH_READER3 mesh_reader("mesh/test/data/cube_2mm_12_elements.vtu") );
        TS_ASSERT_THROWS_ANYTHING( MESH_READER3 mesh_reader("mesh/test/data/nofile.vtu") );
#endif // CHASTE_VTK
    }

    /**
     * Check outputting as a in VTKUnstructuredGrid format.
     */
    void TestOutputVtkUnstructuredGrid(void) throw(Exception)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements.vtu");

        vtkUnstructuredGrid* vtk_unstructed_grid = mesh_reader.OutputMeshAsVtkUnstructuredGrid();

        TS_ASSERT_EQUALS( vtk_unstructed_grid->GetNumberOfPoints(), 12 );
        TS_ASSERT_EQUALS( vtk_unstructed_grid->GetNumberOfCells(), 12 );
#endif // CHASTE_VTK
    }

    /**
     * Check that the nodes are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing nodes) then an exception is thrown.
     */
    void TestGetNextNode(void) throw(Exception)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements.vtu");

        TS_ASSERT_EQUALS( mesh_reader.GetNumNodes(), 12U);

        std::vector<double> first_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( first_node[0] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( first_node[1] , 0.0, 1e-6 );
        TS_ASSERT_DELTA( first_node[2] , 0.0 , 1e-6 );

        std::vector<double> next_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( next_node[0] , 0.2 , 1e-6 );
        TS_ASSERT_DELTA( next_node[1] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( next_node[2] , 0.0 , 1e-6 );

        for (unsigned node=2; node < mesh_reader.GetNumNodes(); node++)
        {
            next_node = mesh_reader.GetNextNode();
        }

        TS_ASSERT_THROWS_THIS( next_node = mesh_reader.GetNextNode(),
                               "Trying to read data for a node that doesn't exist" );
#endif // CHASTE_VTK
    }

    void TestGetNextElementData(void) throw(Exception)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements.vtu");

        TS_ASSERT_EQUALS( mesh_reader.GetNumElements(), 12U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumElementAttributes(), 0u);

        ElementData first_element_data = mesh_reader.GetNextElementData();
        TS_ASSERT_EQUALS( first_element_data.NodeIndices[0] , 11u );
        TS_ASSERT_EQUALS( first_element_data.NodeIndices[1] ,  3u );
        TS_ASSERT_EQUALS( first_element_data.NodeIndices[2] ,  8u );
        TS_ASSERT_EQUALS( first_element_data.NodeIndices[3] ,  0u );
        TS_ASSERT_EQUALS( first_element_data.AttributeValue, 0u );

        ElementData next_element_data = mesh_reader.GetNextElementData();
        TS_ASSERT_EQUALS( next_element_data.NodeIndices[0] , 10u );
        TS_ASSERT_EQUALS( next_element_data.NodeIndices[1] ,  8u );
        TS_ASSERT_EQUALS( next_element_data.NodeIndices[2] ,  5u );
        TS_ASSERT_EQUALS( next_element_data.NodeIndices[3] , 11u );
        TS_ASSERT_EQUALS( next_element_data.AttributeValue, 0u );

        for(unsigned i=2; i<mesh_reader.GetNumElements(); i++)
        {
            next_element_data = mesh_reader.GetNextElementData();
            TS_ASSERT_EQUALS(next_element_data.AttributeValue, 0u);
        }

        TS_ASSERT_THROWS_THIS( ElementData data = mesh_reader.GetNextElementData(),
                               "Trying to read data for an element that doesn't exist" );

        // Test on a .vtu file where the elements are triangles, rather than tetrahedra
        VtkMeshReader<3,3> invalid_mesh_reader("mesh/test/data/sids.vtu");
        TS_ASSERT_EQUALS( invalid_mesh_reader.GetNumElements(), 736U);
        TS_ASSERT_EQUALS( invalid_mesh_reader.GetNumElementAttributes(), 0u);

        TS_ASSERT_THROWS_THIS( first_element_data = invalid_mesh_reader.GetNextElementData(),
                               "Element is not a vtkTetra" );

#endif // CHASTE_VTK
    }

    void TestGetNextFaceData(void) throw(Exception)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements.vtu");

        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 20U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumEdges(), 20U);

        ElementData first_face_data = mesh_reader.GetNextFaceData();
        TS_ASSERT_EQUALS( first_face_data.NodeIndices[0] , 11u );
        TS_ASSERT_EQUALS( first_face_data.NodeIndices[1] ,  3u );
        TS_ASSERT_EQUALS( first_face_data.NodeIndices[2] ,  0u );

        ElementData next_face_data = mesh_reader.GetNextFaceData();
        TS_ASSERT_EQUALS( next_face_data.NodeIndices[0] ,  3u );
        TS_ASSERT_EQUALS( next_face_data.NodeIndices[1] ,  8u );
        TS_ASSERT_EQUALS( next_face_data.NodeIndices[2] ,  0u );

        mesh_reader.Reset();

        first_face_data = mesh_reader.GetNextEdgeData();
        TS_ASSERT_EQUALS( first_face_data.NodeIndices[0] , 11u );
        TS_ASSERT_EQUALS( first_face_data.NodeIndices[1] ,  3u );
        TS_ASSERT_EQUALS( first_face_data.NodeIndices[2] ,  0u );

        next_face_data = mesh_reader.GetNextEdgeData();
        TS_ASSERT_EQUALS( next_face_data.NodeIndices[0] ,  3u );
        TS_ASSERT_EQUALS( next_face_data.NodeIndices[1] ,  8u );
        TS_ASSERT_EQUALS( next_face_data.NodeIndices[2] ,  0u );

        for(unsigned face=2; face<mesh_reader.GetNumFaces(); face++)
        {
            next_face_data = mesh_reader.GetNextEdgeData();
        }

        TS_ASSERT_THROWS_THIS( next_face_data = mesh_reader.GetNextEdgeData(),
                               "Trying to read data for a boundary element that doesn't exist" );
#endif // CHASTE_VTK
    }

    void TestConstructFromVtkUnstructuredGridObject()
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader_1("mesh/test/data/cube_2mm_12_elements.vtu");
        vtkUnstructuredGrid* vtk_unstructed_grid = mesh_reader_1.OutputMeshAsVtkUnstructuredGrid();

        VtkMeshReader<3,3> mesh_reader(vtk_unstructed_grid);

        TS_ASSERT_EQUALS( mesh_reader.GetNumNodes(), 12U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumElements(), 12U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 20U);

        std::vector<double> first_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( first_node[0] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( first_node[1] , 0.0, 1e-6 );
        TS_ASSERT_DELTA( first_node[2] , 0.0 , 1e-6 );

        std::vector<double> next_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( next_node[0] , 0.2 , 1e-6 );
        TS_ASSERT_DELTA( next_node[1] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( next_node[2] , 0.0 , 1e-6 );
#endif // CHASTE_VTK
    }
    void TestGenericReader()
    {
#ifdef CHASTE_VTK
        GenericMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements.vtu");
       
        TS_ASSERT_EQUALS( mesh_reader.GetNumNodes(), 12U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumElements(), 12U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 20U);

        std::vector<double> first_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( first_node[0] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( first_node[1] , 0.0, 1e-6 );
        TS_ASSERT_DELTA( first_node[2] , 0.0 , 1e-6 );

        std::vector<double> next_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( next_node[0] , 0.2 , 1e-6 );
        TS_ASSERT_DELTA( next_node[1] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( next_node[2] , 0.0 , 1e-6 );
#endif // CHASTE_VTK
    }

    /**
     * Check that we can build a TetrahedralMesh using the mesh reader.
     */
    void TestBuildTetrahedralMeshFromMeshReader(void) throw(Exception)
    {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/TestVtkMeshWriter/heart_decimation.vtu");

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 173u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 610u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 312u);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0963, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.3593, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[2], 0.9925, 1e-4);
        
        TS_ASSERT_DELTA(mesh.GetNode(8)->GetPoint()[0], 1.0969, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(8)->GetPoint()[1], 0.6678, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(8)->GetPoint()[2], 0.7250, 1e-4);

        // Check first element has the right nodes
        TetrahedralMesh<3,3>::ElementIterator it = mesh.GetElementIteratorBegin();
        TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(0), 47u);
        TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(1), 31u);
        TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(2), 131u);
        TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(3), 51u);
        TS_ASSERT_EQUALS((it)->GetNode(1), mesh.GetNode(31));


        /**
         * Check that point and cell data attributes work properly.
         */
        // Check element quality cell attribute is read properly
        for (unsigned i = 0; i < 610; i+=60)
        {
            TS_ASSERT_DELTA( mesh_reader.GetCellData( "Quality" )[i], mesh.GetElement(i)->CalculateQuality(), 1e-4 );

        }

        // Check distance from origin point attribute is read properly
        for (unsigned i = 0; i < 173; i+=17)
        {
            TS_ASSERT_DELTA( mesh_reader.GetPointData( "Distance from origin" )[i], norm_2(mesh.GetNode(i)->rGetLocation()), 1e-4 );
        }

        // Check that we can't ask for cell or point data that doesn't exist
        TS_ASSERT_THROWS_ANYTHING( mesh_reader.GetCellData( "Non-existent data" ) );
        TS_ASSERT_THROWS_ANYTHING( mesh_reader.GetPointData( "Non-existent data" ) );
#endif // CHASTE_VTK
    }

    /**
    * Check that we can build a DistributedTetrahedralMesh using the VTK mesh reader.
    */
   void TestBuildDistributedTetrahedralMeshFromVtkMeshReader(void) throw(Exception)
   {
#ifdef CHASTE_VTK
        VtkMeshReader<3,3> mesh_reader("mesh/test/data/TestVtkMeshWriter/heart_decimation.vtu");

        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 173u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 610u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 312u);

        // Check some node co-ordinates
        try
        {
            Node<3> *node = mesh.GetNode(0);
            TS_ASSERT_DELTA(node->GetPoint()[0], 0.0963, 1e-4);
            TS_ASSERT_DELTA(node->GetPoint()[1], 0.3593, 1e-4);
            TS_ASSERT_DELTA(node->GetPoint()[2], 0.9925, 1e-4);
        }
        catch (Exception& e)
        {
            // Don't own this node
        }

        try
        {
            Node<3> *node = mesh.GetNode(8);
            TS_ASSERT_DELTA(node->GetPoint()[0], 1.0969, 1e-4);
            TS_ASSERT_DELTA(node->GetPoint()[1], 0.6678, 1e-4);
            TS_ASSERT_DELTA(node->GetPoint()[2], 0.7250, 1e-4);
        }
        catch (Exception& e)
        {
            // Don't own this node
        }

        // Check first element has the right nodes
        DistributedTetrahedralMesh<3,3>::ElementIterator it = mesh.GetElementIteratorBegin();
        if ( mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(47) )
        {
            //Owner of node 47 has to own (or part own) element 0
            TS_ASSERT_EQUALS((it)->GetIndex(), 0U);
            TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(0), 47u);
            TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(1), 31u);
            TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(2), 131u);
            TS_ASSERT_EQUALS((it)->GetNodeGlobalIndex(3), 51u);

            Node<3> *mesh_node = mesh.GetNodeOrHaloNode(31);
            Node<3> *iterator_node = (it)->GetNode(1);
            TS_ASSERT_EQUALS(iterator_node, mesh_node);
        }
       
#endif // CHASTE_VTK
   }
};

#endif /*TESTVTKMESHREADER_*/
