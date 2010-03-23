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
#ifndef VERTEXMESH3D_HPP_
#define VERTEXMESH3D_HPP_

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMeshWriter;

#include <iostream>
#include <map>
#include <algorithm>

#include <climits> // Work around Boost bug

#include "AbstractMesh.hpp"
#include "ArchiveLocationInfo.hpp"
#include "VertexMeshReader.hpp"
#include "VertexMeshWriter.hpp"
#include "VertexElement.hpp"
#include "VertexElementMap.hpp"

/**
 * A 3d vertex-based mesh class, for use in vertex-based tissue simulations.
 */
class VertexMesh3d : public AbstractMesh<3,3>
{
    friend class TestVertexMesh3d;

protected:

    /** Vector of pointers to VertexElements. */
    std::vector<VertexElement<3,3>*> mElements;

    /** Vector of pointers to VertexElements. */
    std::vector<VertexElement<2,3>*> mFaces;

    /**
     * Solve node mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the node
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Solve element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the element
     */
    unsigned SolveElementMapping(unsigned index) const;

    /**
     * Solve boundary element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the boundary element
     */
    unsigned SolveBoundaryElementMapping(unsigned index) const;

public:

    //////////////////////////////////////////////////////////////////////
    //                             Methods                              //
    //////////////////////////////////////////////////////////////////////

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param faces vector of pointer to VertexElements
     * @param vertexElements vector of pointers to VertexElement<3,3>s
     */
    VertexMesh3d(std::vector<Node<3>*> nodes,
                 std::vector<VertexElement<2,3>*> faces,
                 std::vector<VertexElement<3,3>*> vertexElements);

    /**
     * Destructor.
     */
    virtual ~VertexMesh3d();

    /**
     * @return the number of Nodes in the mesh.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * @return the number of VertexFaces in the mesh.
     */
    virtual unsigned GetNumFaces() const;

    /**
     * @return the number of VertexElements in the mesh.
     */
    virtual unsigned GetNumElements() const;

    /**
     * @return the number of VertexElements in the mesh, including those marked as deleted.
     */
    unsigned GetNumAllElements() const;

    /**
     * @param index  the global index of a specified vertex element
     *
     * @return a pointer to the vertex element
     */
    VertexElement<3,3>* GetElement(unsigned index) const;

    /**
     * Compute the volume of an element.
     *
     * This needs to be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the area of the element
     */
    virtual double GetVolumeOfElement(unsigned index);

    /**
     * Compute the surface area of an element.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param index the global index of a specified vertex element
     *
     * @return the surface area of the element
     */
    double GetSurfaceAreaOfElement(unsigned index);

    /**
     * Compute the unit normal vector to a given face in 3D. This is achieved from a triangle
     * of vertices of the face. Note: this may return the outward or inward normal, depending
     * on your point of view.
     * 
     * @param pFace a face in the mesh
     * 
     * @return the unit normal
     */
    c_vector<double,3> GetUnitNormalToFace(VertexElement<2, 3>* pFace);

    /**
     * Get the area of a given face in 3D. This is achieved by projecting the face onto a 2D plane.
     * To avoid degeneracy and optimize robustness, we choose to ignore the dimension of the component
     * of the unit normal to the plane with the greatest absolute value.
     * 
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     * 
     * @param pFace a face in the mesh
     * 
     * @return the area
     */
    double GetAreaOfFace(VertexElement<2, 3>* pFace);

    /**
     * Compute the centroid of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (centroid_x,centroid_y,centroid_z).
     */
    virtual c_vector<double, 3> GetCentroidOfElement(unsigned index);

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<3,3>& rMeshReader);

    /**
     * Delete mNodes, mFaces and mElements.
     */
    virtual void Clear();
};


#endif /*VERTEXMESH3D_HPP_*/
