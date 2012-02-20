/*

Copyright (C) University of Oxford, 2005-2012

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

#include "CylindricalHoneycombMeshGenerator.hpp"
#include "RandomNumberGenerator.hpp"
#include "MathsCustomFunctions.hpp"

CylindricalHoneycombMeshGenerator::CylindricalHoneycombMeshGenerator(unsigned numNodesAlongWidth, unsigned numNodesAlongLength, unsigned ghosts, double scaleFactor)
{
	mpMesh = NULL;
	mDomainWidth = numNodesAlongWidth*scaleFactor;
	mNumCellWidth = numNodesAlongWidth; //*1 because cells are considered to be size one
	mNumCellLength = numNodesAlongLength;
  // The code below won't work in parallel
    assert(PetscTools::IsSequential());

    // An older version of the constructor might allow the wrong argument through to the scale factor
    assert(scaleFactor > 0.0);

    // Get a unique mesh filename
    std::stringstream pid;
    pid << getpid();
    mMeshFilename = "2D_temporary_honeycomb_mesh_" + pid.str();

    mGhostNodeIndices.empty();

    OutputFileHandler output_file_handler("");
    std::string output_dir = output_file_handler.GetOutputDirectoryFullPath();

    unsigned num_nodes_along_width = mNumCellWidth;
    unsigned num_nodes_along_length = mNumCellLength;
    double horizontal_spacing = mDomainWidth / (double)num_nodes_along_width;
    double vertical_spacing = (sqrt(3)/2)*horizontal_spacing;

    // This line is needed to define ghost nodes later
    mDomainDepth = (double)(num_nodes_along_length) * vertical_spacing;

    // Take account of ghost nodes
    num_nodes_along_length += 2*ghosts;

    unsigned num_nodes            = num_nodes_along_width*num_nodes_along_length;
    unsigned num_elem_along_width = num_nodes_along_width-1;
    unsigned num_elem_along_length = num_nodes_along_length-1;
    unsigned num_elem             = 2*num_elem_along_width*num_elem_along_length;
    unsigned num_edges            = 3*num_elem_along_width*num_elem_along_length + num_elem_along_width + num_elem_along_length;

    double x0 = 0;
    double y0 = -vertical_spacing*ghosts;

    mBottom = -vertical_spacing*ghosts;
    mTop = mBottom + vertical_spacing*(num_nodes_along_length-1);

    // Write node file
    out_stream p_node_file = output_file_handler.OpenOutputFile(mMeshFilename+".node");
    (*p_node_file) << std::scientific;
    (*p_node_file) << std::setprecision(20);
    (*p_node_file) << num_nodes << "\t2\t0\t1" << std::endl;

    unsigned node = 0;
    for (unsigned i=0; i<num_nodes_along_length; i++)
    {
        for (unsigned j=0; j<num_nodes_along_width; j++)
        {
            if (i<ghosts || i>=(ghosts+mNumCellLength))
            {
                mGhostNodeIndices.insert(node);
            }
            unsigned boundary = 0;
            if ((i==0) || (i==num_nodes_along_length-1))
            {
                boundary = 1;
            }

            double x = x0 + horizontal_spacing*((double)j + 0.25*(1.0+ SmallPow(-1,i+1)));
            double y = y0 + vertical_spacing*(double)i;

            // Avoid floating point errors which upset OffLatticeSimulation
            if ( (y<0.0) && (y>-1e-12) )
            {
                // Difficult to cover - just corrects floating point errors that have occurred from time to time!
                #define COVERAGE_IGNORE
                y = 0.0;
                #undef COVERAGE_IGNORE
            }

            (*p_node_file) << node++ << "\t" << x << "\t" << y << "\t" << boundary << std::endl;
        }
    }
    p_node_file->close();

    // Write element file and edge file
    out_stream p_elem_file = output_file_handler.OpenOutputFile(mMeshFilename+".ele");
    (*p_elem_file) << std::scientific;

    out_stream p_edge_file = output_file_handler.OpenOutputFile(mMeshFilename+".edge");
    (*p_node_file) << std::scientific;

    (*p_elem_file) << num_elem << "\t3\t0" << std::endl;
    (*p_edge_file) << num_edges << "\t1" << std::endl;

    unsigned elem = 0;
    unsigned edge = 0;
    for (unsigned i=0; i<num_elem_along_length; i++)
    {
        for (unsigned j=0; j < num_elem_along_width; j++)
        {
            unsigned node0 =     i*num_nodes_along_width + j;
            unsigned node1 =     i*num_nodes_along_width + j+1;
            unsigned node2 = (i+1)*num_nodes_along_width + j;

            if (i%2 != 0)
            {
                node2 = node2 + 1;
            }

            (*p_elem_file) << elem++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;

            unsigned horizontal_edge_is_boundary_edge = 0;
            unsigned vertical_edge_is_boundary_edge = 0;
            if (i==0)
            {
                horizontal_edge_is_boundary_edge = 1;
            }

            (*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << horizontal_edge_is_boundary_edge << std::endl;
            (*p_edge_file) << edge++ << "\t" << node1 << "\t" << node2 << "\t" << 0 << std::endl;
            (*p_edge_file) << edge++ << "\t" << node2 << "\t" << node0 << "\t" << vertical_edge_is_boundary_edge << std::endl;

            node0 = i*num_nodes_along_width + j + 1;

            if (i%2 != 0)
            {
                node0 = node0 - 1;
            }
            node1 = (i+1)*num_nodes_along_width + j+1;
            node2 = (i+1)*num_nodes_along_width + j;

            (*p_elem_file) << elem++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;
        }
    }

    for (unsigned i=0; i<num_elem_along_length; i++)
    {
        unsigned node0, node1;

        if (i%2==0)
        {
             node0 = (i+1)*num_nodes_along_width - 1;
             node1 = (i+2)*num_nodes_along_width - 1;
        }
        else
        {
            node0 = (i+1)*num_nodes_along_width;
            node1 = (i)*num_nodes_along_width;
        }
        (*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << 1 << std::endl;
    }

    for (unsigned j=0; j<num_elem_along_width; j++)
    {
        unsigned node0 = num_nodes_along_width*(num_nodes_along_length-1) + j;
        unsigned node1 = num_nodes_along_width*(num_nodes_along_length-1) + j+1;
        (*p_edge_file) << edge++ << "\t" << node1 << "\t" << node0 << "\t" << 1 << std::endl;
    }

    p_elem_file->close();
    p_edge_file->close();


    // Having written the mesh to file, now construct it using TrianglesMeshReader
    TrianglesMeshReader<2,2> mesh_reader(output_dir + mMeshFilename);
    mpMesh = new Cylindrical2dMesh(mDomainWidth);
    mpMesh->ConstructFromMeshReader(mesh_reader);
    //NodeMap map(mpMesh->GetNumNodes());

    // Make the mesh cylindrical (we use Triangle library mode inside this ReMesh call)
    mpMesh->ReMesh();

    // Delete the temporary files
    std::string command = "rm " + output_dir + mMeshFilename + ".*";
    int return_value = system(command.c_str());
    if (return_value != 0)
    {
        // Can't figure out how to make this throw but seems as if it should be here?
        #define COVERAGE_IGNORE
        EXCEPTION("HoneycombMeshGenerator cannot delete temporary files\n");
        #undef COVERAGE_IGNORE
    }

    // The original files have been deleted, it is better if the mesh object forgets about them
    mpMesh->SetMeshHasChangedSinceLoading();
}

MutableMesh<2,2>* CylindricalHoneycombMeshGenerator::GetMesh()
{
    EXCEPTION("A cylindrical mesh was created but a normal mesh is being requested.");
    return mpMesh; // Not really
}

Cylindrical2dMesh* CylindricalHoneycombMeshGenerator::GetCylindricalMesh()
{
    return (Cylindrical2dMesh*) mpMesh;
}
