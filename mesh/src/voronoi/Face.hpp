/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FACE_HPP_
#define FACE_HPP_

#include "UblasCustomFunctions.hpp"
#include <vector>
#include "Exception.hpp"

template <unsigned DIM>
class Face
{
public:
    /**
     * The vertices of the face, in anticlockwise order. Each vertex must be distinct.
     */
    std::vector< c_vector<double, DIM>* > mVertices;

private:    
    void Increment(typename std::vector< c_vector<double, DIM>* >::iterator& rIterator,
                   Face<DIM>& rFace) const;

public:


    class VertexAndAngle
    {
    public:
        c_vector< double ,DIM >* mpVertex;
        double mAngle;
        bool operator<(const VertexAndAngle& other) const
        {
            return mAngle < other.mAngle;
        }
    };


    /**
     * Compare two faces for equality.
     * Two faces are the same if their vertices differ only by cyclic permutation.
     */
    bool operator==(Face<DIM>& otherFace);
    
    /**
     * Compare two faces for inequality
     */
    bool operator!=(Face<DIM>& otherFace);
    
    /**
     * Return a new face in which the order of the vertices is reversed.
     */
    Face<DIM> operator-();
    
    /**
     * Gets the sum of the length of all edges
     * !!!!! NOTE: Don't use this if you are using a periodic mesh
     * Use GetFacePerimeter(face_index) to takeinto account periodicity.
     */
    double GetPerimeter() const;
    
    /**
     * Gets the area of a voronoi face for 2d space only
     * !!!!! NOTE: Don't use this if you are using a periodic mesh
     * Use GetFaceArea(face_index) to takeinto account periodicity.
     */
    double GetArea() const;
    
    /**
     * Returns number of vertices of a particular face
     */
    unsigned GetNumVertices() const;
    
    /**
     * Returns a vector of vertices
     */
    std::vector< c_vector<double, DIM>* > GetVertices() const;
    
    
    /**
     * Reorders the Vertices of the face anticlockwise
     */
     void OrderVerticesAntiClockwise();
     
     /**
     * @param x x-coordinate
     * @param y y-coordinate
     * @return Polar angle in interval (-PI,PI]
     */
     double ReturnPolarAngle(double x, double y) const;
     
};


#endif /*FACE_HPP_*/
