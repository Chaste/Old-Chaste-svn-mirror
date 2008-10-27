
/*

Copyright (C) University of Oxford, 2008

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
#ifndef HONEYCOMBMESHGENERATOR_CPP_
#define HONEYCOMBMESHGENERATOR_CPP_

#include "HoneycombMeshGenerator.hpp"


HoneycombMeshGenerator::HoneycombMeshGenerator(unsigned numNodesAlongWidth, unsigned numNodesAlongLength, unsigned ghosts, bool cylindrical, double scaleFactor)
  : mpMesh(NULL),
    mCryptWidth(numNodesAlongWidth*scaleFactor),
    mNumCellWidth(numNodesAlongWidth), //*1 because cells are considered to be size one
    mNumCellLength(numNodesAlongLength),
    mCylindrical(cylindrical)
{
    mGhostNodeIndices.empty();

    std::stringstream pid; // Gives a unique filename
    pid << getpid();
    mMeshFilename = "2D_temporary_periodic_crypt_mesh_" + pid.str();
    Make2dPeriodicCryptMesh(mCryptWidth, ghosts);
    OutputFileHandler output_file_handler("");
    std::string output_dir = output_file_handler.GetOutputDirectoryFullPath();

    TrianglesMeshReader<2,2> mesh_reader(output_dir + mMeshFilename);

    if (!mCylindrical)
    {
        mpMesh = new RefinableMesh<2,2>;
        mpMesh->ConstructFromMeshReader(mesh_reader);
    }
    else
    {
        mpMesh = new Cylindrical2dMesh(mCryptWidth);
        mpMesh->ConstructFromMeshReader(mesh_reader);
        NodeMap map(mpMesh->GetNumNodes());
        mpMesh->ReMesh(map); // This makes the mesh cylindrical (uses Triangle library mode inside this ReMesh call).
    }

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

    CancerParameters::Instance()->SetCryptLength(mCryptDepth);
    CancerParameters::Instance()->SetCryptWidth(mCryptWidth);
}
    
    
HoneycombMeshGenerator::~HoneycombMeshGenerator()
{
    delete mpMesh;
}


RefinableMesh<2,2>* HoneycombMeshGenerator::GetMesh()
{
    if (mCylindrical)
    {
        EXCEPTION("A cylindrical mesh was created but a normal mesh is being requested.");
    }
    return mpMesh;
}


Cylindrical2dMesh* HoneycombMeshGenerator::GetCylindricalMesh()
{
    if (!mCylindrical)
    {
        EXCEPTION("A normal mesh was created but a cylindrical mesh is being requested.");
    }
    return (Cylindrical2dMesh*) mpMesh;
}


std::set<unsigned> HoneycombMeshGenerator::GetGhostNodeIndices()
{
    return mGhostNodeIndices;
}


RefinableMesh<2,2>* HoneycombMeshGenerator::GetCircularMesh(double radius)
{
    assert(!mCylindrical); // Following call only safe if is not a cylindrical mesh

    // Centre the mesh at (0,0)
    c_vector<double,2> centre = zero_vector<double>(2);
    for (unsigned i=0; i<mpMesh->GetNumNodes(); i++)
    {
        centre += mpMesh->GetNode(i)->rGetLocation();
    }
    centre /= (double)mpMesh->GetNumNodes();

    mpMesh->Translate(-centre[0], -centre[1]);

    // Iterate over nodes, deleting any that lie more
    // than the specified radius from (0,0)
    for (unsigned i=0; i<mpMesh->GetNumAllNodes(); i++)
    {
        if ( norm_2(mpMesh->GetNode(i)->rGetLocation()) >= radius)
        {
            mpMesh->DeleteNodePriorToReMesh(i);
        }
        else
        {
            // Jiggle the data
            c_vector<double,2>& r_location = mpMesh->GetNode(i)->rGetModifiableLocation();
            c_vector<double,2> shift;
            RandomNumberGenerator *p_gen = RandomNumberGenerator::Instance();
            double max_jiggle = radius*5e-6;
            shift[0] = max_jiggle*(p_gen->ranf()-0.5);
            shift[1] = max_jiggle*(p_gen->ranf()-0.5);
            r_location += shift;
        }
    }

    // Remesh
    NodeMap map(mpMesh->GetNumNodes());
    mpMesh->ReMesh(map);

    return mpMesh;
}

void HoneycombMeshGenerator::Make2dPeriodicCryptMesh(double width, unsigned ghosts)
{
    OutputFileHandler output_file_handler("");

    if (output_file_handler.IsMaster())
    {
        out_stream p_node_file = output_file_handler.OpenOutputFile(mMeshFilename+".node");
        (*p_node_file) << std::scientific;

        out_stream p_elem_file = output_file_handler.OpenOutputFile(mMeshFilename+".ele");
        (*p_elem_file) << std::scientific;

        unsigned numNodesAlongWidth = mNumCellWidth;
        unsigned numNodesAlongLength = mNumCellLength;
        double horizontal_spacing = width / (double)numNodesAlongWidth;
        double vertical_spacing = (sqrt(3)/2)*horizontal_spacing;

        // This line needed to define ghost nodes later...
        mCryptDepth = (double)(numNodesAlongLength) * vertical_spacing;

        // Add in the ghost nodes...
        if (!mCylindrical)
        {
            numNodesAlongWidth = numNodesAlongWidth + 2*ghosts;
        }
        numNodesAlongLength = numNodesAlongLength + 2*ghosts;

        unsigned num_nodes            = numNodesAlongWidth*numNodesAlongLength;
        unsigned num_elem_along_width = numNodesAlongWidth-1;
        unsigned num_elem_along_length = numNodesAlongLength-1;
        unsigned num_elem             = 2*num_elem_along_width*num_elem_along_length;
        unsigned num_edges            = 3*num_elem_along_width*num_elem_along_length + num_elem_along_width + num_elem_along_length;

        double x0 = -horizontal_spacing*ghosts;
        if (mCylindrical)
        {
            x0 = 0;
        }
        double y0 = -vertical_spacing*ghosts;
        mBottom = -vertical_spacing*ghosts;
        mTop = mBottom + vertical_spacing*(numNodesAlongLength-1);

        (*p_node_file) << num_nodes << "\t2\t0\t1" << std::endl;
        unsigned node = 0;

        for (unsigned i=0; i<numNodesAlongLength; i++)
        {
            for (unsigned j=0; j<numNodesAlongWidth; j++)
            {
                if ( i<ghosts || i>=(ghosts+mNumCellLength))
                {
                    mGhostNodeIndices.insert(node);
                }
                else if ( !mCylindrical && (j < ghosts || j >= (ghosts+mNumCellWidth)))
                {
                    mGhostNodeIndices.insert(node);
                }
                unsigned boundary = 0;
                if ((i==0) || (i==numNodesAlongLength-1))
                {
                    boundary = 1;
                }
                if (!mCylindrical)
                {
                    if ((j==0) || (j==numNodesAlongWidth-1))
                    {
                        boundary = 1;
                    }
                }

                double x = x0 + horizontal_spacing*((double)j + 0.25*(1.0+ pow(-1,i+1)));
                double y = y0 + vertical_spacing*(double)i;

                // Avoid floating point errors which upset CryptSimulation2d
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
                unsigned node0 =     i*numNodesAlongWidth + j;
                unsigned node1 =     i*numNodesAlongWidth + j+1;
                unsigned node2 = (i+1)*numNodesAlongWidth + j;

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
                if (j==0 && i%2==0 && !mCylindrical)
                {
                    vertical_edge_is_boundary_edge = 1;
                }

                (*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 <<  "\t" << horizontal_edge_is_boundary_edge << std::endl;
                (*p_edge_file) << edge++ << "\t" << node1 << "\t" << node2 <<  "\t" << 0 << std::endl;
                (*p_edge_file) << edge++ << "\t" << node2 << "\t" << node0 <<  "\t" << vertical_edge_is_boundary_edge << std::endl;

                node0 = i*numNodesAlongWidth + j + 1;

                if (i%2 != 0)
                {
                    node0 = node0 - 1;
                }
                node1 = (i+1)*numNodesAlongWidth + j+1;
                node2 = (i+1)*numNodesAlongWidth + j;

                (*p_elem_file) << elem++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;
            }
        }

        for (unsigned i=0; i<num_elem_along_length; i++)
        {
            unsigned node0, node1;

            if (i%2==0)
            {
                 node0 = (i+1)*numNodesAlongWidth - 1;
                 node1 = (i+2)*numNodesAlongWidth - 1;
            }
            else
            {
                node0 = (i+1)*numNodesAlongWidth;
                node1 = (i)*numNodesAlongWidth;
            }
            (*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << 1 << std::endl;
        }

        for (unsigned j=0; j<num_elem_along_width; j++)
        {
            unsigned node0 = numNodesAlongWidth*(numNodesAlongLength-1) + j;
            unsigned node1 = numNodesAlongWidth*(numNodesAlongLength-1) + j+1;
            (*p_edge_file) << edge++ << "\t" << node1 << "\t" << node0 << "\t" << 1 << std::endl;
        }

        p_elem_file->close();
        p_edge_file->close();
    }

    // Wait for the new mesh to be available
#ifndef SPECIAL_SERIAL
    PetscTools::Barrier();
#endif //SPECIAL_SERIAL
#define COVERAGE_IGNORE
}
#undef COVERAGE_IGNORE


#endif //HONEYCOMBMESHGENERATOR_CPP_

