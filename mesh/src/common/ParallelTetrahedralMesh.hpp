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

#include "TetrahedralMesh.hpp"
#include "Element.hpp"
#include "DistributedVector.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ParallelTetrahedralMesh : public TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>
{
    
private: 

    unsigned mTotalNumElements; 
    unsigned mTotalNumNodes;
    
    std::map< unsigned, Element<ELEMENT_DIM, SPACE_DIM>* > mElements;
    std::map< unsigned, Node<SPACE_DIM>* > mNodes;
    std::map< unsigned, Node<SPACE_DIM>* > mGhostNodes;
    std::map< unsigned, BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* > mBoundaryElements;

    std::map< unsigned,Node<SPACE_DIM> *> mBoundaryNodes;
    

public:

//    ParallelTetrahedralMesh();

    virtual ~ParallelTetrahedralMesh();

    void ComputeMeshPartioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
                               std::set<unsigned>& nodesOwned,
                               std::set<unsigned>& ghostNodesOwned,
                               std::set<unsigned>& elementsOwned) const;    

    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader,
                                 bool cullInternalFaces=false);

    unsigned GetNumLocalNodes();

    unsigned GetNumLocalElements();

    unsigned GetNumNodes();

    unsigned GetNumElements();

    
};

//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ParallelTetrahedralMesh()
//{
//}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~ParallelTetrahedralMesh()
{
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ComputeMeshPartioning(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
    std::set<unsigned>& nodesOwned,
    std::set<unsigned>& ghostNodesOwned,
    std::set<unsigned>& elementsOwned) const
{
    ///\todo: add a timing event for the partitioning
    
    // Not calling ParMETIS yet, dumb partitioning (for the moment)
    DistributedVector::SetProblemSize(mTotalNumNodes);
    for(DistributedVector::Iterator node_number = DistributedVector::Begin(); node_number != DistributedVector::End(); ++node_number)
    {
         nodesOwned.insert(node_number.Global);
    }
    
    std::vector<unsigned> node_indices;
    for(unsigned element_number = 0; element_number < mTotalNumElements; element_number++)
    {
        node_indices = rMeshReader.GetNextElement();
        for(unsigned node_index = 0; node_index < node_indices.size(); node_index++)
        {
            if (nodesOwned.find(node_indices[node_index]) != nodesOwned.end())
            {
                elementsOwned.insert(element_number);
                
                std::set<unsigned> temp;
                std::set_difference(node_indices.begin(), node_indices.end(),
                                      nodesOwned.begin(), nodesOwned.end(),
                                      std::inserter(temp, temp.begin()) );                              
                               
                ghostNodesOwned.insert(temp.begin(), temp.end()); //need to remove the owned nodes from here
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
    mTotalNumElements = rMeshReader.GetNumElements();
    mTotalNumNodes = rMeshReader.GetNumNodes();    
    
    std::set<unsigned> nodes_owned;
    std::set<unsigned> ghost_nodes_owned;
    std::set<unsigned> elements_owned;
    
    ComputeMeshPartioning(rMeshReader, nodes_owned, ghost_nodes_owned, elements_owned);
    
    // Load the nodes owned by the processor
    std::vector<double> coords;
    for (unsigned node_index=0; node_index < mTotalNumNodes; node_index++)
    {
        coords = rMeshReader.GetNextNode();
        
        // The node is owned by the processor
        if (nodes_owned.find(node_index) != nodes_owned.end())
        {
            mNodes[node_index] = new Node<SPACE_DIM>(node_index, coords, false);
            //continue;
        }

        // The node is a ghost node in this processor
        if (ghost_nodes_owned.find(node_index) != ghost_nodes_owned.end())
        {
            mGhostNodes[node_index] = new Node<SPACE_DIM>(node_index, coords, false);
        }

    }

    // Load the elements owned by the processor
    std::vector<unsigned> node_indices;
    for (unsigned element_index=0; element_index < mTotalNumElements; element_index++)
    {
        node_indices = rMeshReader.GetNextElement();

        // The element is owned by the processor
        if (elements_owned.find(element_index) != elements_owned.end())
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            for (unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                if (nodes_owned.find(node_indices[j]) != nodes_owned.end())
                {
                    nodes.push_back(mNodes[node_indices[j]]);
                }
                else
                {
                    nodes.push_back(mGhostNodes[node_indices[j]]);
                }                    
            }
            mElements[element_index] = new Element<ELEMENT_DIM,SPACE_DIM>(element_index, nodes);
        }
    }
    
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalNodes()
{
    return mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalElements()
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes()
{
    return mTotalNumNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements()
{
    return mTotalNumElements;
}



#endif /*PARALLELTETRAHEDRALMESH_HPP_*/
