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


#ifndef VORONOITESSELLATION_HPP_
#define VORONOITESSELLATION_HPP_

#include "UblasCustomFunctions.hpp"
#include "TetrahedralMesh.hpp"
#include "VoronoiCell.hpp"


/**
 * Voronoi tessellation class. For use in certain mesh-based tissue simulations.
 */
template<unsigned DIM>
class VoronoiTessellation
{
private:

    friend class TestVoronoiTessellation;
    friend class InventorVoronoiWriter;

    /** A tetrahedral mesh. */
    TetrahedralMesh<DIM,DIM>& mrMesh;

    /** Vertices correspond to elements of the mesh. */
    std::vector< c_vector<double,DIM>* > mVertices;

    /** Faces corespond to edges of the mesh. */
    std::vector< Face<DIM>* > mFaces;

    /** Cells correspond to nodes of the mesh. */
    std::vector< VoronoiCell > mVoronoiCells;

    /** Set of node indices corresponding to non-ghost nodes. */
    std::set<unsigned> mLocationIndices;

    /**
     * Generate the vertices of the tessellation using the
     * circumcentres of the mesh elements.
     */
    void GenerateVerticesFromElementCircumcentres();

    /**
     * Use a tetrahedral mesh to initialise the faces and vertices
     * of the tessellation.
     *
     * @param rMesh the mesh
     */
    void Initialise(TetrahedralMesh<DIM, DIM>& rMesh);

public:

    /**
     * Constructor. Create a tessellation of the given mesh which must be Delaunay
     * (see TetrahedralMesh::CheckIsVoronoi).
     *
     * @param rMesh a tetrahedral mesh
     * @param locationIndices an optional vector of location indices that correspond to non-ghost nodes
     */
    VoronoiTessellation(TetrahedralMesh<DIM,DIM>& rMesh, const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Destructor.
     */
    ~VoronoiTessellation();

    /**
     * Get a VoronoiCell.
     *
     * @param index The index of the cell is the index of the corresponding node in the original mesh.
     * If the corresponding node was on the boundary, this will return a cell with no faces.
     */
    const VoronoiCell& rGetCell(unsigned index) const;

    /**
     * Get the face of the VoronoiCell with a given index.
     *
     * @param index  The index of the cell is the index of the corresponding node in the original mesh.
     */
    const Face<DIM>& rGetFace(unsigned index) const;

    /**
     * Get the number of faces in the tessellation.
     */
    unsigned GetNumFaces() const;

    /**
     * Get the area of the face with a given index.
     *
     * @param index
     */
    double GetAreaOfElement(unsigned index) const;

    /**
     * Get the perimeter of the face with a given index.
     *
     * @param index
     */
    double GetPerimeterOfElement(unsigned index) const;

    /**
     * Get the length of the tessellation edge between two given nodes.
     *
     * @param nodeIndex1
     * @param nodeIndex2
     */
    double GetEdgeLength(unsigned nodeIndex1, unsigned nodeIndex2) const;

    /**
     * Get the number of vertices in the tessellation.
     */
    unsigned GetNumVertices() const;

    /**
     * Get the vertex with a given index.
     *
     * @param index
     */
    c_vector<double,DIM>* GetVertex(unsigned index);

    /**
     * Get the number of VoronoiCells in the tessellation.
     */
    unsigned GetNumCells();

};

#endif /*VORONOITESSELLATION_HPP_*/
