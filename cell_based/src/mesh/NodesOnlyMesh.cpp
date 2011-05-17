/*

Copyright (C) University of Oxford, 2005-2011

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

#include "NodesOnlyMesh.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ConstructNodesWithoutMesh(const std::vector<Node<SPACE_DIM>*> & rNodes)
{
    this->Clear();
    for (unsigned i=0; i<rNodes.size(); i++)
    {
        assert(!rNodes[i]->IsDeleted());
        bool boundary = rNodes[i]->IsBoundaryNode();
        c_vector<double, SPACE_DIM> location=rNodes[i]->rGetLocation();
        
        Node<SPACE_DIM>* p_node_copy = new Node<SPACE_DIM>(i, location, boundary);
        this->mNodes.push_back( p_node_copy );
        mCellRadii.push_back(1.0);
    }
}

template<unsigned SPACE_DIM>
double NodesOnlyMesh<SPACE_DIM>::GetCellRadius(unsigned index)
{
    return mCellRadii[index];    
}
    
template<unsigned SPACE_DIM>     
void NodesOnlyMesh<SPACE_DIM>::SetCellRadius(unsigned index, double radius)
{
    mCellRadii[index] = radius;
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class NodesOnlyMesh<1>;
template class NodesOnlyMesh<2>;
template class NodesOnlyMesh<3>;
