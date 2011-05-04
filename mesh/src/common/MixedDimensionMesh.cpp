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

#include "MixedDimensionMesh.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::MixedDimensionMesh(DistributedTetrahedralMeshPartitionType::type partitioningMethod)
    : DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::DistributedTetrahedralMesh(partitioningMethod)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::~MixedDimensionMesh()
{
    for (unsigned i=0; i<mCableElements.size(); i++)
    {
        delete mCableElements[i];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader)
{
    DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ConstructFromMeshReader(rMeshReader);
    
    // Add cable elements
    mNumCableElements = rMeshReader.GetNumCableElements();
    //this->mCableElements.reserve(mNumCableElements);
    
    for (unsigned element_index=0; element_index < mNumCableElements; element_index++)
    {
        ElementData element_data = rMeshReader.GetNextCableElementData();

        //Determine if we own any nodes on this cable element
        bool node_owned = false;
        for (unsigned j=0; j<2; j++) // cables are always 1d
        {
            try
            {
                this->SolveNodeMapping(element_data.NodeIndices[j]);
                node_owned = true;
                break;
            }
            catch (Exception &e)
            {
                //We deal with non-owned nodes in the next part
            }
        }
        
        //If we don't locally own either node, then we don't construct the cable      
        if (node_owned)
        {           
            std::vector<Node<SPACE_DIM>*> nodes;
            nodes.reserve(2u);
    
            for (unsigned j=0; j<2; j++) // cables are always 1d
            {
                //Note (#1760) that if we own one node on a cable element then we are likely to own the other.
                //If not, we are likely to have a halo.
                //If not, (free-running Purkinje with monodomain mesh?), then this will throw.
                nodes.push_back(this->GetNodeOrHaloNode(element_data.NodeIndices[j]) );
            }
    
            Element<1u,SPACE_DIM>* p_element = new Element<1u,SPACE_DIM>(element_index, nodes);
            this->mCableElements.push_back(p_element);
    
            if (rMeshReader.GetNumCableElementAttributes() > 0)
            {
                assert(rMeshReader.GetNumCableElementAttributes() == 1);
                unsigned attribute_value = element_data.AttributeValue;
                p_element->SetRegion(attribute_value);
            }
        }
    }

    rMeshReader.Reset();
}
    
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::GetNumCableElements() const
{
   return mNumCableElements;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalCableElements() const
{
   return mCableElements.size();
}
    
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<1u, SPACE_DIM>* MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::GetCableElement(unsigned index) const
{
    assert(index < mCableElements.size());
    return mCableElements[index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::CalculateDesignatedOwnershipOfCableElement( unsigned elementIndex )
{
    //This may throw in the distributed parallel case
    unsigned tie_break_index = this->GetCableElement(elementIndex)->GetNodeGlobalIndex(0);

    //if it is in my range
    if (this->GetDistributedVectorFactory()->IsGlobalIndexLocal(tie_break_index))
    {
        return true;
    }
    else
    {
        return false;
    }
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::CableElementIterator MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::GetCableElementIteratorBegin() const
{
    return mCableElements.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::CableElementIterator MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::GetCableElementIteratorEnd() const
{
    return mCableElements.end();
}



/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class MixedDimensionMesh<1,1>;
template class MixedDimensionMesh<1,2>;
template class MixedDimensionMesh<1,3>;
template class MixedDimensionMesh<2,2>;
template class MixedDimensionMesh<2,3>;
template class MixedDimensionMesh<3,3>;
