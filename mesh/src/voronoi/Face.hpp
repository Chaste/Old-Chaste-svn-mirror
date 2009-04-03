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
     * Increment the Face vertex iterator.
     * 
     * @param rIterator the Face vertex iterator
     * @param rFace the Face
     */
    void Increment(typename std::vector< c_vector<double, DIM>* >::iterator& rIterator,
                   Face<DIM>& rFace) const;

public:

    /**
     * Helper class containing a pointer to a vertex of the face 
     * and the polar angle from the centre of the face to this 
     * vertex.
     * 
     * \todo This is duplicated in the VoronoiTessellation class; move to a separate file?
     */  
    class VertexAndAngle
    {
    public:

        /** Pointer to a vertex. */
        c_vector<double, DIM>* mpVertex;

        /** Polar angle. */
        double mAngle; 

        /**
         * Less-than angle comparison operator.
         * 
         * @param other the VertexAndAngle object to compare to
         */
        bool operator<(const VertexAndAngle& rOther) const
        {
            return mAngle < rOther.mAngle;
        }
    };

    /**
     * The vertices of the face, in anticlockwise order. Each vertex 
     * must be distinct.
     * 
     * This member variable is public as it is accessed directly by 
     * VoronoiTessellation methods.
     */
    std::vector< c_vector<double, DIM>* > mVertices;

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
     * Get the sum of the length of all edges of the Face.
     * 
     * NOTE: Don't use this if you are using a periodic mesh. Use 
     * GetFacePerimeter(face_index) to take into account periodicity.
     */
    double GetPerimeter() const;

    /**
     * Get the area of the Face (works in 2D only).
     * 
     * NOTE: Don't use this if you are using a periodic mesh. Use 
     * GetFaceArea(face_index) to takeinto account periodicity.
     */
    double GetArea() const;

    /**
     * Return number of vertices of the Face.
     */
    unsigned GetNumVertices() const;

    /**
     * Return a vector of vertices owned by the Face.
     */
    std::vector< c_vector<double, DIM>* > GetVertices() const;

    /**
     * Reorder the vertices of the Face anticlockwise.
     */
    void OrderVerticesAntiClockwise();

    /**
     * Return the polar angle of the point (x,y).
     * 
     * \todo This is duplicated in the VoronoiTessellation class; move to a separate file?
     * 
     * @param x x-coordinate
     * @param y y-coordinate
     * @return Polar angle in interval (-PI,PI]
     */
    double ReturnPolarAngle(double x, double y) const;

};

#endif /*FACE_HPP_*/
