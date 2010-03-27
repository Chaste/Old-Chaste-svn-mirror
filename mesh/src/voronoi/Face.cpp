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

#include "Face.hpp"
#include <list>

/**
 * Global method allowing alist of pairs (c_vector<double, DIM>*, double) to be compared
 * in terms of their second entry and std::list.sort() to be called.
 */
template<unsigned DIM>
bool VertexAngleComparison(const std::pair<c_vector<double, DIM>*, double> lhs, const std::pair<c_vector<double, DIM>*, double> rhs)
{
    return lhs.second < rhs.second;
}

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void Face<DIM>::Increment(typename std::vector< c_vector<double, DIM>* >::iterator& rIterator,
                          Face<DIM>& rFace) const
{
    rIterator++;
    if (rIterator == rFace.mVertices.end())
    {
        rIterator = rFace.mVertices.begin();
    }
}

template<unsigned DIM>
bool Face<DIM>::operator==(Face<DIM>& rOtherFace)
{
    typename std::vector< c_vector<double, DIM>* >::iterator this_iterator = mVertices.begin();
    typename std::vector< c_vector<double, DIM>* >::iterator other_iterator = rOtherFace.mVertices.begin();

    // Find first vertex
    while ( this_iterator != mVertices.end() &&
            other_iterator != rOtherFace.mVertices.end() &&
            norm_2(**this_iterator - **other_iterator) >1e-10 )
    {
        this_iterator++;
    }
    if (this_iterator==mVertices.end() || other_iterator==rOtherFace.mVertices.end())
    {
        // Could not find first vertex; faces are distinct unless they are empty
        return ( this_iterator==mVertices.end() && other_iterator==rOtherFace.mVertices.end() );
    }

    typename std::vector< c_vector<double, DIM>* >::iterator this_start = this_iterator;
    Increment(this_iterator, *this);
    Increment(other_iterator, rOtherFace);

    // Check remanining vertices are equal
    while (this_iterator != this_start)
    {
        if (norm_2(**this_iterator - **other_iterator) > 1e-10)
        {
            return false;
        }
        else
        {
            Increment(this_iterator, *this);
            Increment(other_iterator, rOtherFace);
        }
    }
    return (other_iterator == rOtherFace.mVertices.begin());
}

#define COVERAGE_IGNORE // Spuriously not covered
template<unsigned DIM>
bool Face<DIM>::operator!=(Face& rOtherFace)
{
   return !(*this == rOtherFace);
}

#undef COVERAGE_IGNORE

template<unsigned DIM>
Face<DIM> Face<DIM>::operator-()
{
   Face<DIM> reversed_face;
   typename std::vector< c_vector<double, DIM>* >::iterator this_iterator = mVertices.end();
   while (this_iterator != mVertices.begin())
   {
       this_iterator--;
       reversed_face.mVertices.push_back(*this_iterator);
   }
   return reversed_face;
}

template<unsigned DIM>
unsigned Face<DIM>::GetNumVertices() const
{
    return mVertices.size();
}

template<unsigned DIM>
void Face<DIM>::OrderVerticesAntiClockwise()
{
    assert(mVertices.size() > 1);

    // Compute the location of the centre of the face
    c_vector<double,DIM> centre = zero_vector<double>(DIM);
    for (unsigned j=0; j<mVertices.size(); j++)
    {
        centre += *(mVertices[j]);
    }
    centre /= mVertices.size();

    // Compute and store the polar angle from the centre to each vertex
    std::list<std::pair<c_vector<double, DIM>*, double> > vertex_angle_list;
    for (unsigned j=0; j<mVertices.size(); j++)
    {
        c_vector<double, DIM> centre_to_vertex = *(mVertices[j]) - centre;
        double angle = atan2(centre_to_vertex(1), centre_to_vertex(0));

        std::pair<c_vector<double, DIM>*, double> pair(mVertices[j], angle);
        vertex_angle_list.push_back(pair);
    }

    // Sort the list in order of increasing angle
    vertex_angle_list.sort(VertexAngleComparison<DIM>);

    // Use polar angles to reorder mVertices anticlockwise
    mVertices.clear();
    for (typename std::list<std::pair<c_vector<double, DIM>*, double> >::iterator list_iter = vertex_angle_list.begin();
         list_iter != vertex_angle_list.end();
         ++list_iter)
    {
        mVertices.push_back(list_iter->first);
    }
}

template<unsigned DIM>
void Face<DIM>::AddVertex(c_vector<double, DIM>* pVertex)
{
    mVertices.push_back(pVertex);
}

template<unsigned DIM>
unsigned Face<DIM>::GetNumVertices()
{
    return mVertices.size();
}

template<unsigned DIM>
c_vector<double, DIM>& Face<DIM>::rGetVertex(unsigned index)
{
    return *(mVertices[index]);
}

template<unsigned DIM>
void Face<DIM>::SetVertex(unsigned index, c_vector<double, DIM>* pNewLocation)
{
    mVertices[index] = pNewLocation;
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class Face<1>;
template class Face<2>;
template class Face<3>;
