/*

Copyright (C) University of Oxford, 2008

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

#ifndef PARALLELTETRAHEDRALMESH_HPP_
#define PARALLELTETRAHEDRALMESH_HPP_

#include "AbstractMesh.hpp"
#include "AbstractMeshReader.hpp"
#include "Element.hpp"
#include "BoundaryElement.hpp"
#include "Node.hpp"
#include "DistributedVector.hpp"
#include "PetscTools.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ParallelTetrahedralMesh : public AbstractMesh< ELEMENT_DIM, SPACE_DIM>
{
    
private: 

    unsigned mTotalNumElements;
    unsigned mTotalNumBoundaryElements; 
    unsigned mTotalNumNodes;
    
    std::vector<Node<SPACE_DIM>* > mGhostNodes;
    
    std::map<unsigned, unsigned> mNodesMapping;
    std::map<unsigned, unsigned> mGhostNodesMapping;    
    std::map<unsigned, unsigned> mElementsMapping;        
    
public:

//    ParallelTetrahedralMesh();

    virtual ~ParallelTetrahedralMesh();

    void ComputeMeshPartioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
                               std::set<unsigned>& rNodesOwned,
                               std::set<unsigned>& rGhostNodesOwned,
                               std::set<unsigned>& rElementsOwned) const;    

    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader,
                                 bool cullInternalFaces=false);

    unsigned GetNumLocalNodes() const;

    unsigned GetNumLocalElements() const;

    unsigned GetNumNodes() const;

    unsigned GetNumElements() const;
    
    unsigned GetNumBoundaryElements() const;
    
    BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* GetBoundaryElement(unsigned index) const;
    
    void SetElementOwnerships(unsigned lo, unsigned hi);

private:

    void RegisterNode(unsigned index);
    void RegisterGhostNode(unsigned index);
    void RegisterElement(unsigned index);

    unsigned SolveNodeMapping(unsigned index);
    unsigned SolveGhostNodeMapping(unsigned index);
    unsigned SolveElementMapping(unsigned index);            
};

//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ParallelTetrahedralMesh()
//{
//}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~ParallelTetrahedralMesh()
{
    for (unsigned i=0; i<this->mGhostNodes.size(); i++)
    {
        delete this->mGhostNodes[i];
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ComputeMeshPartioning(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
    std::set<unsigned>& rNodesOwned,
    std::set<unsigned>& rGhostNodesOwned,
    std::set<unsigned>& rElementsOwned) const
{
    ///\todo: add a timing event for the partitioning
    
    // Not calling ParMETIS yet, dumb partitioning (for the moment)
    DistributedVector::SetProblemSize(mTotalNumNodes);
    for(DistributedVector::Iterator node_number = DistributedVector::Begin(); node_number != DistributedVector::End(); ++node_number)
    {
         rNodesOwned.insert(node_number.Global);
    }
    
    for(unsigned element_number = 0; element_number < mTotalNumElements; element_number++)
    {
        ElementData element_data = rMeshReader.GetNextElementData();

        for(unsigned i=0; i<ELEMENT_DIM+1; i++)
        {
            if (rNodesOwned.find(element_data.NodeIndices[i]) != rNodesOwned.end())
            {
                rElementsOwned.insert(element_number);
                
                std::set<unsigned> temp; /// \todo: there should be a way of avoiding the use of temp

                std::set_difference(element_data.NodeIndices.begin(), element_data.NodeIndices.end(),
                                    rNodesOwned.begin(), rNodesOwned.end(),
                                    std::inserter(temp, temp.begin()) );                              
                               
                rGhostNodesOwned.insert(temp.begin(), temp.end());
            }
        }
    }
    
    rMeshReader.Reset();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
    bool cullInternalFaces)
{
    
    if(ELEMENT_DIM==1)
    {
        cullInternalFaces = true;
    }    
    
    mTotalNumElements = rMeshReader.GetNumElements();
    mTotalNumNodes = rMeshReader.GetNumNodes();
    mTotalNumBoundaryElements = rMeshReader.GetNumFaces();   
        
    std::set<unsigned> nodes_owned;
    std::set<unsigned> ghost_nodes_owned;
    std::set<unsigned> elements_owned;
    
    ComputeMeshPartioning(rMeshReader, nodes_owned, ghost_nodes_owned, elements_owned);
    
    // Reserve memory
    this->mElements.reserve(elements_owned.size());
    this->mNodes.reserve(nodes_owned.size());
        
    // Load the nodes owned by the processor
    std::vector<double> coords;
    for (unsigned node_index=0; node_index < mTotalNumNodes; node_index++)
    {
        /// \todo: assert the node is not considered both owned and ghost-owned. Remove continue statement few lines below then.
        coords = rMeshReader.GetNextNode();
        
        // The node is owned by the processor
        if (nodes_owned.find(node_index) != nodes_owned.end())
        {
            RegisterNode(node_index);
            this->mNodes.push_back(new Node<SPACE_DIM>(node_index, coords, false));
            continue;
        }

        // The node is a ghost node in this processor
        if (ghost_nodes_owned.find(node_index) != ghost_nodes_owned.end())
        {
            RegisterGhostNode(node_index);
            mGhostNodes.push_back(new Node<SPACE_DIM>(node_index, coords, false));
        }

    }

    // Load the elements owned by the processor
    for (unsigned element_index=0; element_index < mTotalNumElements; element_index++)
    {
        ElementData element_data = rMeshReader.GetNextElementData();

        // The element is owned by the processor
        if (elements_owned.find(element_index) != elements_owned.end())
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            unsigned node_local_index;
            for (unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                if (nodes_owned.find(element_data.NodeIndices[j]) != nodes_owned.end())
                {
                    node_local_index = SolveNodeMapping(element_data.NodeIndices[j]);
                    nodes.push_back(this->mNodes[node_local_index]);
                }
                else
                {
                    node_local_index = SolveGhostNodeMapping(element_data.NodeIndices[j]);
                    nodes.push_back(this->mGhostNodes[node_local_index]);
                }                    
            }
            
            RegisterElement(element_index);
            
            Element<ELEMENT_DIM,SPACE_DIM>* p_element = new Element<ELEMENT_DIM,SPACE_DIM>(element_index, nodes);
            this->mElements.push_back(p_element);
            
            if (rMeshReader.GetNumElementAttributes() > 0)
            {
                assert(rMeshReader.GetNumElementAttributes() == 1);
                unsigned attribute_value = element_data.AttributeValue;
                p_element->SetRegion(attribute_value);
            }
        }
    }
    
    // Boundary nodes and elements    
    unsigned actual_face_index = 0; 
    for (unsigned face_index=0; face_index<mTotalNumBoundaryElements; face_index++)
    {
        std::vector<unsigned> node_indices = rMeshReader.GetNextFace();

        bool own = false;

        for (unsigned node_index=0; node_index<node_indices.size(); node_index++)
        {
            // if I own this node
            if(mNodesMapping.find(node_indices[node_index]) != mNodesMapping.end())
            {
                own = true;
            }
        }
        
        if (!own )
        {
            continue;
        }
            
        
        bool is_boundary_face = true;

        // Determine if this is a boundary face
        std::set<unsigned> containing_element_indices; // Elements that contain this face
        std::vector<Node<SPACE_DIM>*> nodes;

        for (unsigned node_index=0; node_index<node_indices.size(); node_index++)
        {
            // if I own this node
            if(mNodesMapping.find(node_indices[node_index]) != mNodesMapping.end())
            {
                // Add Node pointer to list for creating an element
                unsigned node_local_index = SolveNodeMapping(node_indices[node_index]);
                nodes.push_back(this->mNodes[node_local_index]);
            }

            // if I ghost-own this node
            if(mGhostNodesMapping.find(node_indices[node_index]) != mGhostNodesMapping.end())
            {
                // Add Node pointer to list for creating an element
                unsigned node_local_index = SolveGhostNodeMapping(node_indices[node_index]);                
                nodes.push_back(this->mGhostNodes[node_local_index]); 
            }
            
            if(cullInternalFaces)
            {
                // Work out what elements contain this face, by taking the intersection
                // of the sets of elements containing each node in the face.
                if (node_index == 0)
                {
                    containing_element_indices = nodes[node_index]->rGetContainingElementIndices();
                }
                else
                {
                    std::set<unsigned> temp;
                    std::set_intersection(nodes[node_index]->rGetContainingElementIndices().begin(),
                                          nodes[node_index]->rGetContainingElementIndices().end(),
                                          containing_element_indices.begin(), containing_element_indices.end(),
                                          std::inserter(temp, temp.begin()));
                    containing_element_indices = temp;
                }
            }
        }
       

        if(cullInternalFaces)
        {
            // only if not 1D as this assertion does not apply to quadratic 1D meshes
            if(ELEMENT_DIM!=1)
            {
                //If the following assertion is thrown, it means that the .edge/.face file does not
                //match the .ele file -- they were generated at separate times.  Simply remove the internal
                //edges/faces by hand.
                assert(containing_element_indices.size() != 0);
            }

            // if num_containing_elements is greater than 1, it is not an boundary face
            if(containing_element_indices.size() > 1)
            {
                is_boundary_face = false;
            }
            
            // in 1D QUADRATICS, all nodes are faces, so internal nodes which don't have any
            // containing elements must also be unmarked as a boundary face
            if( (ELEMENT_DIM==1) && (containing_element_indices.size()==0))
            {
                is_boundary_face = false;
            }
        }

        if (is_boundary_face)
        {
            // This is a boundary face
            // Ensure all its nodes are marked as boundary nodes
            for (unsigned j=0; j<nodes.size(); j++)
            {
                if (!nodes[j]->IsBoundaryNode())
                {
                    nodes[j]->SetAsBoundaryNode();
                    this->mBoundaryNodes.push_back(nodes[j]);
                }
                //Register the index that this bounday element will have
                //with the node
                nodes[j]->AddBoundaryElement(actual_face_index);
            }

            // The added elements will be deleted in our destructor
            this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(actual_face_index, nodes));
            actual_face_index++;
        }
    }        
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalNodes() const
{
    return this->mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalElements() const
{
    return this->mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mTotalNumNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mTotalNumElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements() const
{
    return mTotalNumBoundaryElements;
}

//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//Element<ELEMENT_DIM, SPACE_DIM>* ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
//{
//    /// \todo: assert ownership
//    return (this->mElements[index]);
//}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElement(unsigned index) const
{
    /// \todo: assert ownership
    return (this->mBoundaryElements[index]);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetElementOwnerships(unsigned lo, unsigned hi)
{
    // all the local nodes are owned by the processor (obviously...)    
    assert(hi>=lo);
    for (unsigned element_index=0; element_index<this->mElements.size(); element_index++)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element=this->mElements[element_index];
        p_element->SetOwnership(true);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RegisterNode(unsigned index)
{
    mNodesMapping[index] = this->mNodes.size();
}    

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RegisterGhostNode(unsigned index)
{
    mGhostNodesMapping[index] = mGhostNodes.size();    
}    

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RegisterElement(unsigned index)
{
    mElementsMapping[index] = this->mElements.size();    
}    

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index)
{
    std::map<unsigned, unsigned>::const_iterator node_position = mNodesMapping.find(index); 
        
    if(node_position == mNodesMapping.end())
    {
        std::stringstream message;
        message << "Requested node " << index << " does not belong to processor " << PetscTools::GetMyRank();
        EXCEPTION(message.str().c_str());        
    }
    return node_position->second;    
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveGhostNodeMapping(unsigned index)
{
    assert(mGhostNodesMapping.find(index) != mGhostNodesMapping.end());
    return mGhostNodesMapping[index];    
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index)
{
    std::map<unsigned, unsigned>::const_iterator element_position = mElementsMapping.find(index); 
    
    if(element_position == mElementsMapping.end())
    {
        std::stringstream message;
        message << "Requested element " << index << " does not belong to processor " << PetscTools::GetMyRank();
        EXCEPTION(message.str().c_str());        
    }

    return element_position->second;    
}            


#endif /*PARALLELTETRAHEDRALMESH_HPP_*/
