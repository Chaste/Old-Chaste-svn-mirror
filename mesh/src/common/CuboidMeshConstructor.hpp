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


#ifndef CUBOIDMESHCONSTRUCTOR_HPP_
#define CUBOIDMESHCONSTRUCTOR_HPP_

#include "TrianglesMeshWriter.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "TetrahedralMesh.hpp"

//const double mesh_width = 0.2; // cm

/**
 * Helper class for constructing cuboidal meshes.
 */
template<unsigned DIM>
class CuboidMeshConstructor
{
private:

    /**
     * Construct a one-dimensional linear mesh.
     *
     * @param rMesh  The mesh
     * @param width  Width of the mesh
     */
    void ConstructHyperCube(TetrahedralMesh<1,1>& rMesh, unsigned width);

    /**
     * Construct a two-dimensional rectangular mesh.
     *
     * @param rMesh  The mesh
     * @param width  Width of the mesh
     */
    void ConstructHyperCube(TetrahedralMesh<2,2>& rMesh, unsigned width);

    /**
     * Construct a three-dimensional cuboidal mesh.
     *
     * @param rMesh  The mesh
     * @param width  Width of the mesh
     */
    void ConstructHyperCube(TetrahedralMesh<3,3>& rMesh, unsigned width);

public:

    double mMeshWidth; /**< Width of the mesh. */
    unsigned NumElements; /**< Number of elements in the mesh. \todo Should be mNumElements */
    unsigned NumNodes; /**< Number of nodes in the mesh. \todo Should be mNumNodes  */

    /**
     * Construct the mesh.
     *
     * @param meshNum  Index for the mesh
     * @param meshWidth  Width of the mesh
     */
    std::string Construct(unsigned meshNum, double meshWidth);

    /**
     * Get the width of the mesh.
     */
    double GetWidth();

};

#endif /*CUBOIDMESHCONSTRUCTOR_HPP_*/
