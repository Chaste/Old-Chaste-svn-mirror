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


#include "CuboidMeshConstructor.hpp"


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CuboidMeshConstructor<ELEMENT_DIM, SPACE_DIM>::ConstructHyperCube(TetrahedralMesh<1,SPACE_DIM>& rMesh, unsigned width)
{
    rMesh.ConstructLinearMesh(width);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CuboidMeshConstructor<ELEMENT_DIM, SPACE_DIM>::ConstructHyperCube(TetrahedralMesh<2,2>& rMesh, unsigned width)
{
    rMesh.ConstructRectangularMesh(width, width);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CuboidMeshConstructor<ELEMENT_DIM, SPACE_DIM>::ConstructHyperCube(TetrahedralMesh<3,3>& rMesh, unsigned width)
{
    rMesh.ConstructCuboid(width, width, width);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string CuboidMeshConstructor<ELEMENT_DIM, SPACE_DIM>::Construct(unsigned meshNum, double meshWidth)
{
    //The class has only been tested for ELEMENT_DIM == SPACE_DIM or ELEMENT_DIM == 1 && SPACE_DIM == 3
    assert(ELEMENT_DIM == SPACE_DIM || (ELEMENT_DIM == 1 && SPACE_DIM == 3));

    mMeshWidth = meshWidth;
    assert(meshNum < 30); //Sanity
    const std::string mesh_dir = "ConvergenceMesh";
    OutputFileHandler output_file_handler(mesh_dir);

    // Create the mesh
    unsigned mesh_size = (unsigned) pow(2, meshNum+2); // number of elements in each dimension
    double scaling = mMeshWidth/(double) mesh_size;
    TetrahedralMesh<ELEMENT_DIM,SPACE_DIM> mesh;
    ConstructHyperCube(mesh, mesh_size);
    mesh.Scale(scaling, scaling, scaling);
    mNumElements = mesh.GetNumElements();
    mNumNodes = mesh.GetNumNodes();
    std::stringstream file_name_stream;
    file_name_stream << "cube_" << ELEMENT_DIM << "D_2mm_" << mNumElements << "_elements";
    std::string mesh_filename = file_name_stream.str();

    TrianglesMeshWriter<ELEMENT_DIM,SPACE_DIM> mesh_writer(mesh_dir, mesh_filename, false);
    mesh_writer.WriteFilesUsingMesh(mesh);

    std::string mesh_pathname = output_file_handler.GetOutputDirectoryFullPath() + mesh_filename;

    return mesh_pathname;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CuboidMeshConstructor<ELEMENT_DIM, SPACE_DIM>::GetWidth()
{
    return mMeshWidth;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned CuboidMeshConstructor<ELEMENT_DIM, SPACE_DIM>::GetNumElements()
{
    return mNumElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned CuboidMeshConstructor<ELEMENT_DIM, SPACE_DIM>::GetNumNodes()
{
    return mNumNodes;
}



//////////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////////

template class CuboidMeshConstructor<1>;
template class CuboidMeshConstructor<2>;
template class CuboidMeshConstructor<3>;
template class CuboidMeshConstructor<1,3>;
