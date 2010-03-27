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


#ifndef FACE_HPP_
#define FACE_HPP_

#include "UblasCustomFunctions.hpp"
#include "Exception.hpp"
#include <vector>


/**
 * A face class for use in the VoronoiTessellation class.
 */
template <unsigned DIM>
class Face
{
private:

    /**
     * The vertices of the face, in anticlockwise order. Each vertex
     * must be distinct.
     *
     * This member variable is public as it is accessed directly by
     * VoronoiTessellation methods.
     */
    std::vector< c_vector<double, DIM>* > mVertices;

    /**
     * Increment the Face vertex iterator.
     *
     * @param rIterator the Face vertex iterator
     * @param rFace the Face
     */
    void Increment(typename std::vector< c_vector<double, DIM>* >::iterator& rIterator,
                   Face<DIM>& rFace) const;

public:

    /**
     * Compare two faces for equality. Two faces are the same if their
     * vertices differ only by cyclic permutation.
     *
     * @param rOtherFace the Face to compare to
     */
    bool operator==(Face<DIM>& rOtherFace);

    /**
     * Compare two faces for inequality.
     *
     * @param rOtherFace the Face to compare to
     */
    bool operator!=(Face<DIM>& rOtherFace);

    /**
     * Return a new face in which the order of the vertices is reversed.
     */
    Face<DIM> operator-();

    /**
     * Return number of vertices of the Face.
     */
    unsigned GetNumVertices() const;

    /**
     * Reorder the vertices of the Face anticlockwise.
     */
    void OrderVerticesAntiClockwise();

    /**
     * Add a vertex to the Face.
     *
     * @param pVertex the location of the new vertex
     */
    void AddVertex(c_vector<double, DIM>* pVertex);

    /**
     * @return the number of vertices in the Face.
     */
    unsigned GetNumVertices();

    /**
     * Get the Vertex with a given index.
     *
     * @param index the index of the Vertex in the Face
     */
    c_vector<double, DIM>& rGetVertex(unsigned index);

    /**
     * Reset the location of the Vertex with a given index.
     *
     * @param index the index of the Vertex in the Face
     * @param pNewLocation the new location of the Vertex
     */
    void SetVertex(unsigned index, c_vector<double, DIM>* pNewLocation);
};

#endif /*FACE_HPP_*/
