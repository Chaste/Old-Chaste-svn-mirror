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
void MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader)
{
    rMeshReader.Reset();
    this->mMeshFileBaseName = rMeshReader.GetMeshFileBaseName();

    // Record number of corner nodes
    const unsigned num_nodes = rMeshReader.GetNumNodes();

    // Reserve memory for nodes, so we don't have problems with pointers stored in elements becoming invalid.
    this->mNodes.reserve(num_nodes);

    // Add nodes
    std::vector<double> coords;
    for (unsigned i=0; i < num_nodes; i++)
    {
        coords = rMeshReader.GetNextNode();
        Node<SPACE_DIM>* p_node =  new Node<SPACE_DIM>(i, coords, false);

        for (unsigned i = 0; i < rMeshReader.GetNodeAttributes().size(); i++)
        {
            double attribute = rMeshReader.GetNodeAttributes()[i];
            p_node->AddNodeAttribute(attribute);
        }

        this->mNodes.push_back(p_node);
    }

    // Add elements
    this->mElements.reserve(rMeshReader.GetNumElements());
    for (unsigned element_index=0; element_index < rMeshReader.GetNumElements(); element_index++)
    {
        ElementData element_data = rMeshReader.GetNextElementData();
        std::vector<Node<SPACE_DIM>*> nodes;

        // NOTE: currently just reading element vertices from mesh reader - even if it
        // does contain information about internal nodes (ie for quadratics) this is
        // ignored here and used elsewhere: ie don't do this:
        //   unsigned nodes_size = node_indices.size();

        for (unsigned j=0; j<ELEMENT_DIM+1; j++) // num vertices=ELEMENT_DIM+1, may not be equal to nodes_size.
        {
            assert(element_data.NodeIndices[j] <  this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        Element<ELEMENT_DIM,SPACE_DIM>* p_element = new Element<ELEMENT_DIM,SPACE_DIM>(element_index, nodes);
        this->mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = element_data.AttributeValue;
            p_element->SetRegion(attribute_value);
        }
    }

    // Add boundary elements and nodes
    for (unsigned face_index=0; face_index<(unsigned)rMeshReader.GetNumFaces(); face_index++)
    {
        ElementData face_data = rMeshReader.GetNextFaceData();
        std::vector<unsigned> node_indices = face_data.NodeIndices;

        // NOTE: as above just read boundary element *vertices* from mesh reader - even if
        // it is a quadratic mesh with internal elements, the extra nodes are ignored here
        // and used elsewhere: ie, we don't do this:
        //   unsigned nodes_size = node_indices.size();

        std::vector<Node<SPACE_DIM>*> nodes;
        for (unsigned node_index=0; node_index<ELEMENT_DIM; node_index++) // node_index from 0 to DIM-1, not 0 to node.size()-1
        {
            assert(node_indices[node_index] < this->mNodes.size());
            // Add Node pointer to list for creating an element
            nodes.push_back(this->mNodes[node_indices[node_index]]);
        }

        // This is a boundary face
        // Ensure all its nodes are marked as boundary nodes

        assert(nodes.size()==ELEMENT_DIM); // just taken vertices of boundary node from
        for (unsigned j=0; j<nodes.size(); j++)
        {
            if (!nodes[j]->IsBoundaryNode())
            {
                nodes[j]->SetAsBoundaryNode();
                this->mBoundaryNodes.push_back(nodes[j]);
            }
            // Register the index that this bounday element will have with the node
            nodes[j]->AddBoundaryElement(face_index);
        }

        // The added elements will be deleted in our destructor
        BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* p_boundary_element = new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(face_index, nodes);
        this->mBoundaryElements.push_back(p_boundary_element);

        if (rMeshReader.GetNumFaceAttributes() > 0)
        {
            assert(rMeshReader.GetNumFaceAttributes() == 1);
            unsigned attribute_value = face_data.AttributeValue;
            p_boundary_element->SetRegion(attribute_value);
        }
    }

    // Add cable elements
    this->mCableElements.reserve(rMeshReader.GetNumCableElements());
    for (unsigned element_index=0; element_index < rMeshReader.GetNumCableElements(); element_index++)
    {
        ElementData element_data = rMeshReader.GetNextCableElementData();
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.reserve(2u);

        for (unsigned j=0; j<2; j++) // cables are always 1d
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
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

    rMeshReader.Reset();
}
    
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::GetNumCableElements() const
{
   return mCableElements.size();
}
    
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<1u, SPACE_DIM>* MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::GetCableElement(unsigned index) const
{
    assert(index < mCableElements.size());
    return mCableElements[index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size() );
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size() );
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    assert(index < this->mBoundaryElements.size() );
    return index;
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class MixedDimensionMesh<2,2>;
template class MixedDimensionMesh<3,3>;
