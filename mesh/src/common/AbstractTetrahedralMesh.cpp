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

#include "AbstractTetrahedralMesh.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetElementOwnerships(unsigned lo, unsigned hi)
{
    assert(hi >= lo);
    for (unsigned element_index=0; element_index<mElements.size(); element_index++)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element = mElements[element_index];
        p_element->SetOwnership(false);
        for (unsigned local_node_index=0; local_node_index< p_element->GetNumNodes(); local_node_index++)
        {
            unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
            if (lo<=global_node_index && global_node_index<hi)
            {
                p_element->SetOwnership(true);
                break;
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::AbstractTetrahedralMesh()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~AbstractTetrahedralMesh()
{
    // Iterate over elements and free the memory
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    // Iterate over boundary elements and free the memory
    for (unsigned i=0; i<mBoundaryElements.size(); i++)
    {
        delete mBoundaryElements[i];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements() const
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllBoundaryElements() const
{
    return mBoundaryElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements() const
{
    return mBoundaryElements.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<ELEMENT_DIM, SPACE_DIM>* AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    unsigned local_index = SolveElementMapping(index);
    return mElements[local_index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElement(unsigned index) const
{
    unsigned local_index = SolveBoundaryElementMapping(index);
    return mBoundaryElements[local_index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElementIteratorBegin() const
{
    return mBoundaryElements.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElementIteratorEnd() const
{
    return mBoundaryElements.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetInverseJacobianForElement(
        unsigned elementIndex,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian,
        double& rJacobianDeterminant,
        c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian) const
{
    mElements[SolveElementMapping(elementIndex)]->CalculateInverseJacobian(rJacobian, rJacobianDeterminant, rInverseJacobian);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWeightedDirectionForBoundaryElement(
        unsigned elementIndex,
        c_vector<double, SPACE_DIM>& rWeightedDirection,
        double& rJacobianDeterminant) const
{
    mBoundaryElements[SolveBoundaryElementMapping(elementIndex)]->CalculateWeightedDirection(rWeightedDirection, rJacobianDeterminant );
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class AbstractTetrahedralMesh<1,1>;
template class AbstractTetrahedralMesh<1,2>;
template class AbstractTetrahedralMesh<1,3>;
template class AbstractTetrahedralMesh<2,2>;
template class AbstractTetrahedralMesh<2,3>;
template class AbstractTetrahedralMesh<3,3>;
