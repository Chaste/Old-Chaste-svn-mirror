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


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template <unsigned DIM>
void Face<DIM>::Increment(typename std::vector< c_vector<double, DIM>* >::iterator& rIterator,
                          Face<DIM>& rFace) const
{
    rIterator++;
    if (rIterator==rFace.mVertices.end() )
    {
        rIterator = rFace.mVertices.begin();
    }
};

template <unsigned DIM>
bool Face<DIM>::operator==(Face<DIM>& otherFace)
{
    typename std::vector< c_vector<double, DIM>* >::iterator this_iterator=mVertices.begin();
    typename std::vector< c_vector<double, DIM>* >::iterator other_iterator=otherFace.mVertices.begin();
    // find first vertex
    while ( this_iterator!=mVertices.end() &&
            other_iterator!=otherFace.mVertices.end() &&
            norm_2(**this_iterator - **other_iterator) >1e-10 )
    {
        this_iterator++;
    }
    if (this_iterator==mVertices.end() || other_iterator==otherFace.mVertices.end())
    {
        // could not find first vertex
        // faces are distinct unless they are empty
        return this_iterator==mVertices.end() && other_iterator==otherFace.mVertices.end();
    }

    typename std::vector< c_vector<double, DIM>* >::iterator this_start=this_iterator;
    Increment(this_iterator, *this);
    Increment(other_iterator, otherFace);

    // check remanining vertices are equal
    while (this_iterator!=this_start)
    {
        if (norm_2(**this_iterator - **other_iterator) > 1e-10)
        {
            return false;
        }
        else
        {
            Increment(this_iterator, *this);
            Increment(other_iterator, otherFace);
        }
    }
    return other_iterator==otherFace.mVertices.begin();
};

#define COVERAGE_IGNORE //Spuriously not covered
template <unsigned DIM>
bool Face<DIM>::operator!=(Face& otherFace)
{
   return !(*this==otherFace);
};
#undef COVERAGE_IGNORE

template <unsigned DIM>
Face<DIM> Face<DIM>::operator-()
{
   Face<DIM> reversed_face;
   typename std::vector< c_vector<double, DIM>* >::iterator this_iterator=mVertices.end();
   while (this_iterator != mVertices.begin())
   {
       this_iterator--;
       reversed_face.mVertices.push_back(*this_iterator);
   }
   return reversed_face;
};

template <unsigned DIM>
double Face<DIM>::GetPerimeter() const
{
    double perimeter_return = 0;
    for(unsigned i=0; i<mVertices.size(); i++)
    {
        perimeter_return += norm_2(*mVertices[i]-*mVertices[(i+1)%mVertices.size()]);
    }
    return perimeter_return;
};

template <unsigned DIM>
double Face<DIM>::GetArea() const
{
    #define COVERAGE_IGNORE
    assert(DIM==2);
    #undef COVERAGE_IGNORE

    double area_return = 0;
    for(unsigned i=0; i<mVertices.size(); i++)
    {
        //  Area = sum ( x_i * y_i+1 - y_i * x_i+1 )/2.0 over all vertices,
        //      assuming vertices are ordered anti-clockwise
        area_return +=   ( (*mVertices[i])(0) * (*mVertices[(i+1)%mVertices.size()])(1)
                          -(*mVertices[i])(1) * (*mVertices[(i+1)%mVertices.size()])(0) ) / 2.0;
    }
    return area_return;
};

template <unsigned DIM>
unsigned Face<DIM>::GetNumVertices() const
{
    return mVertices.size();
};

template <unsigned DIM>
std::vector< c_vector<double, DIM>* > Face<DIM>::GetVertices() const
{
    return mVertices;
};

template<unsigned DIM>
double Face<DIM>::ReturnPolarAngle(double x, double y) const
{
    if (x==0)
    {
        if (y>0)
        {
            return M_PI/2.0;
        }
        else if (y<0)
        {
            return -M_PI/2.0;
        }
        else
        {
            EXCEPTION("Tried to compute polar angle of (0,0)");
        }
    }

    double angle = atan(y/x);

    if (y >= 0 && x < 0 )
    {
        angle += M_PI;
    }
    else if (y < 0 && x < 0 )
    {
        angle -= M_PI;
    }
    return angle;
};

template <unsigned DIM>
void Face<DIM>::OrderVerticesAntiClockwise()
{
    // Reorder mVertices Anticlockwise
    std::vector< VertexAndAngle > vertices_and_angles;

    c_vector<double,DIM> centre = zero_vector<double>(DIM);

    for(unsigned j=0; j<mVertices.size(); j++)
    {
        //std::cout<< "\n mVertices" << j << " " << (*(mVertices[j]))(0) << std::flush;
        centre += *(mVertices[j]);
    }

    centre /= mVertices.size();
    for(unsigned j=0; j<mVertices.size(); j++)
    {

        VertexAndAngle va;
        c_vector<double, DIM> centre_to_vertex;
        //assert(centre  *(mVertices[j]) );
        centre_to_vertex = *(mVertices[j]) - centre;
        va.mAngle = ReturnPolarAngle(centre_to_vertex(0), centre_to_vertex(1));
        va.mpVertex = mVertices[j];
        vertices_and_angles.push_back(va);
    }

    std::sort(vertices_and_angles.begin(), vertices_and_angles.end());

    // Create face
    mVertices.clear();
    for ( typename std::vector< VertexAndAngle >::iterator vertex_iterator = vertices_and_angles.begin();
          vertex_iterator !=vertices_and_angles.end();
          vertex_iterator++)
    {
        mVertices.push_back(vertex_iterator->mpVertex);
    }
};

#endif /*FACE_HPP_*/
