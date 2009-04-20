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

#include "NodeBox.hpp"

template<unsigned DIM>
NodeBox<DIM>::NodeBox(c_vector<double, 2*DIM> minAndMaxValues)
{
    mMinAndMaxValues = minAndMaxValues;
}

template<unsigned DIM>
c_vector<double, 2*DIM>& NodeBox<DIM>::rGetMinAndMaxValues()
{
    return mMinAndMaxValues;
}

template<unsigned DIM>
void NodeBox<DIM>::AddNode(Node<DIM>* p_node)
{
    mNodesContained.insert(p_node);
}
   
   
template<unsigned DIM>
void NodeBox<DIM>::RemoveNode(Node<DIM>* p_node)
{
    mNodesContained.erase(p_node);
}

template<unsigned DIM>
std::set< Node<DIM>* >& NodeBox<DIM>::rGetNodesContained()
{
    return mNodesContained;
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class NodeBox<1>;
template class NodeBox<2>;
template class NodeBox<3>;
