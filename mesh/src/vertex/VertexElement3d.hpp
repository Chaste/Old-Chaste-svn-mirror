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

#ifndef VERTEXELEMENT3D_HPP_
#define VERTEXELEMENT3D_HPP_

#include "UblasCustomFunctions.hpp"
#include "AbstractElement.hpp"
#include "VertexElement.hpp"
#include <cmath>


/**
 * A 3d Vertex element class
 */
class VertexElement3d : public AbstractElement<3,3>
{
private:

    /**
     * Faces of the Vertex Element, which should be distinct.
     */
    std::vector<VertexElement<2,3>*> mFaces;

    /**
     * How each face is oriented. From the perspective of the centre
     * of the VoronoiCell, the vertices of each face should be ordered
     * anti clockwise. If and only if this is false, the order of vertices
     * in the corresponding face should be reversed.
     *
     * N.B. Most faces belong to two VoronoiCell, but with opposite
     * orientations. This allows us to reuse the face data across the
     * two cells.
     */
    std::vector<bool> mOrientations;

public:
	/**
	 * Default constructor.
	 *
	 * @param index global index of the element
	 * @param nodes vector of pointers to nodes
	 * @param faces vector of pointers to VertexElements
	 * @param orientations vector of how each face is oriented.
	 */
	VertexElement3d(unsigned index,
					std::vector<Node<3>*> nodes,
					std::vector<VertexElement<2, 3>*> faces,
					std::vector<bool> orientations);

	/**
	 * Default constructor for use by serializer.
	 */
	VertexElement3d();

	/**
	 * Destructor.
	 */
	virtual ~VertexElement3d();

    /**
     * Get the number of faces in the VoronoiCell.
     */
    unsigned GetNumFaces() const;

    /**
       * @param index  the global index of a specified face
       *
       * @return a pointer to the face
       */
    VertexElement<2,3>* GetFace(unsigned index) const;

    /**
     * Get whether the face with a given index is oriented clockwise.
     *
     * @param index the index of the face in the VoronoiCell
     */
    bool FaceIsOrientatedClockwise(unsigned index) const;

    /**
	 * Overridden RegisterWithNodes() method.
	 *
	 * Informs all nodes forming this element that they are in this element.
	 */
	void RegisterWithNodes();

	/**
	 * Overridden MarkAsDeleted() method.
	 *
	 * Mark an element as having been removed from the mesh.
	 * Also notify nodes in the element that it has been removed.
	 */
	void MarkAsDeleted();

	/**
	 * Update node at the given index.
	 *
	 * @param rIndex is an local index to which node to change
	 * @param pNode is a pointer to the replacement node
	 */
	void UpdateNode(const unsigned& rIndex, Node<3>* pNode);
};

#endif /*VERTEXELEMENT3D_HPP_*/
