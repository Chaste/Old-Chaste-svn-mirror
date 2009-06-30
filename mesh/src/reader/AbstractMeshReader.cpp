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

#include "AbstractMeshReader.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElementAttributes() const
{
    // By default returns 0.  If a concrete class does read attributes
    // it needs to overload this method.
    return 0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumFaceAttributes() const
{
    // By default returns 0.  If a concrete class does read attributes
    // it needs to overload this method.
    return 0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumEdges() const
{
    return GetNumFaces();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextEdge()
{
    return GetNextFaceData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>::GetMeshFileBaseName()
{
    return "";
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class AbstractMeshReader<0,1>;
template class AbstractMeshReader<1,1>;
template class AbstractMeshReader<1,2>;
template class AbstractMeshReader<1,3>;
template class AbstractMeshReader<2,2>;
template class AbstractMeshReader<2,3>;
template class AbstractMeshReader<3,3>;
