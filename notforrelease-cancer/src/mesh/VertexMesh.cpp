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

#include "VertexMesh.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(std::vector<Node<SPACE_DIM>*> nodes, 
                                               std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> vertexElements,
                                               double thresholdDistance)
    : mThresholdDistance(thresholdDistance),
      mAddedNodes(true)
{
    Clear();
    for (unsigned node_index=0; node_index<nodes.size(); node_index++)
    {
        Node<SPACE_DIM>* temp_node = nodes[node_index];
        mNodes.push_back(temp_node);
    }
    
    for (unsigned elem_index=0; elem_index<vertexElements.size(); elem_index++)
    {
        VertexElement<ELEMENT_DIM,SPACE_DIM>* temp_vertex_element = vertexElements[elem_index];
        mElements.push_back(temp_vertex_element);
    }
    
    SetupVertexElementsOwnedByNodes();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(double thresholdDistance)
    : mThresholdDistance(thresholdDistance),
      mAddedNodes(false)
{
    assert(thresholdDistance > 0.0);
    Clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(unsigned numAcross, unsigned numUp, double thresholdDistance)
    : mThresholdDistance(thresholdDistance),
      mAddedNodes(true)
{
    if (SPACE_DIM==2)
    {    
        assert(numAcross > 1);
        unsigned node_index = 0;
        
        // Create the nodes
        for (unsigned j=0; j<=2*numUp+1; j++)
        {
            if (j%2 == 0)
            {
                for (unsigned i=1; i<=3*numAcross+1; i+=2)
                {
                    if (j!=0 || i!= 3*numAcross+1)
                    {
                        if (i%3 != 2)
                        {
                            Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index, false, i/(2.0*sqrt(3)),j/2.0);
                            mNodes.push_back(p_node);
                            node_index++;
                        }
                    }
                }
            }
            else 
            {
                for (unsigned i=0; i<=3*numAcross+1; i+=2)
                {
                    if ((j!=2*numUp+1 || i != 0) && (j!=2*numUp+1 || i!= 3*numAcross+1))
                    {
                        if (i%3 != 2)
                        {
                            Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index, false, i/(2.0*sqrt(3)), j/2.0);
                            mNodes.push_back(p_node);
                            node_index++;
                        }
                    }
                }
            }
        }  
        
        // Create the elements. The array node_indices contains the 
        // global node indices from bottom left, going anticlockwise.
        
        unsigned node_indices[6];
        unsigned element_index;
        
        for (unsigned j=0; j<numUp; j++)
        {
            for (unsigned i=0; i<numAcross; i++)
            {
                element_index = j*numAcross + i;
                
                if (numAcross%2==0) // numAcross is even
                {
                    if (j == 0)     // bottom row
                    {
                        if (i%2 == 0) // even
                        {
                            node_indices[0] = i;
                        }
                        else // odd
                        {
                            node_indices[0] = numAcross+i;
                        }                                           
                    }                       
                    else    // not on the bottom row 
                    {
                         if (i%2 == 0) // even
                        {
                            node_indices[0] = (2*numAcross+1)+2*(j-1)*(numAcross+1)+i;
                        }
                        else // odd
                        {
                            node_indices[0] = (2*numAcross+1)+(2*j-1)*(numAcross+1)+i;
                        }                        
                    }
                        
                }
                else // numAcross is odd
                {
                    if (i%2 == 0) // even
                    {
                        node_indices[0] = 2*j*(numAcross+1)+i;
                    }
                    else // odd
                    {
                        node_indices[0] = (2*j+1)*(numAcross+1)+i;
                    }
                }
                node_indices[1] = node_indices[0] + 1;
                node_indices[2] = node_indices[0] + numAcross + 2;
                node_indices[3] = node_indices[0] + 2*numAcross + 3;
                node_indices[4] = node_indices[0] + 2*numAcross + 2;
                node_indices[5] = node_indices[0] + numAcross + 1;
                 
                if ((j==numUp-1)&&(i%2 == 1))
                {
                    // On top row and its an odd column nodes 
                    node_indices[3] -= 1;
                    node_indices[4] -= 1;
                }
                  
                if ((j==0)&&(i%2 == 0)&&(numAcross%2==0))
                {
                    // On bottom row and its an even column and there is
                    // an even number of columns in total, (i.e. the very bottom) 
                    node_indices[2] -= 1;
                    node_indices[3] -= 1;
                    node_indices[4] -= 1;
                    node_indices[5] -= 1;
                }

                std::vector<Node<SPACE_DIM>*> element_nodes;
                
                for (int i=0; i<6; i++)
                {
                   element_nodes.push_back(mNodes[node_indices[i]]);
                }
                VertexElement<ELEMENT_DIM,SPACE_DIM>* p_element = new VertexElement<ELEMENT_DIM,SPACE_DIM>(element_index, element_nodes);
                mElements.push_back(p_element);
            }
        }  
    
    //    HoneycombMeshGenerator generator(numAcross+1,numUp+1,0,false);
    //    MutableMesh<2,2>* p_mesh = generator.GetMesh();
    //    VoronoiTessellation<2> tessellation(*p_mesh);
    //    
    //    for (unsigned i = 0;i<tessellation.GetNumVertices();i++)
    //    {
    //        c_vector<double,2>* position = tessellation.GetVertex(i);
    //        Node<2>* p_node = new Node<2>(0, false, (*position)(0), (*position)(1));
    //        mNodes.push_back(p_node);
    //    }    
    
        /// \todo: loop over the p_mesh's nodes, and if it is a non-boundary node create a VertexElement using
        //         the corresponding cell. Then get rid of the nodes in mNodes that do not belong in any cell.
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::~VertexMesh()
{
    Clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetThresholdDistance() const
{
    return mThresholdDistance;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetThresholdDistance(double thresholdDistance)
{
    mThresholdDistance = thresholdDistance;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetupVertexElementsOwnedByNodes()
{
    for (unsigned index=0; index<mElements.size(); index++)
    {
        VertexElement<ELEMENT_DIM,SPACE_DIM>* p_temp_vertex_element = mElements[index];
        for (unsigned node_index=0; node_index<p_temp_vertex_element->GetNumNodes(); node_index++)
        {
            Node<SPACE_DIM>* p_temp_node = p_temp_vertex_element->GetNode(node_index);
            p_temp_node->AddElement(p_temp_vertex_element->GetIndex());
        }
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    mDeletedNodeIndices.clear();
    mAddedNodes = false;
    
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        delete mNodes[i];
    }

    mNodes.clear();
    mElements.clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mNodes.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElements.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index) const
{
    assert(index < mNodes.size());
    return mNodes[index];
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM,SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM> *pNewNode)
{
    if (mDeletedNodeIndices.empty())
    {
        pNewNode->SetIndex(mNodes.size());
        mNodes.push_back(pNewNode);
    }
    else
    {
        unsigned index = mDeletedNodeIndices.back();
        pNewNode->SetIndex(index);
        mDeletedNodeIndices.pop_back();
        delete mNodes[index];
        mNodes[index] = pNewNode;
    }
    mAddedNodes = true;
    return pNewNode->GetIndex();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::AddElement(VertexElement<ELEMENT_DIM,SPACE_DIM> *pNewElement)
{
    if (mDeletedElementIndices.empty())
    {
        pNewElement->SetIndex(mElements.size());
        mElements.push_back(pNewElement);
    }
    else
    {
        unsigned index = mDeletedElementIndices.back();
        pNewElement->SetIndex(index);
        mDeletedElementIndices.pop_back();
        delete mElements[index];
        mElements[index] = pNewElement;
    }
    
    mAddedElements = true;
    pNewElement->RegisterWithNodes();
    
    return pNewElement->GetIndex();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point)
{
    mNodes[nodeIndex]->SetPoint(point);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideEdge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();
    
    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));
    
    // Check that the nodes have a common edge
    assert(shared_elements.size()>0);
    
    // Create a new node (position is not important as it will be changed)
    Node<SPACE_DIM>* p_new_node = new Node<SPACE_DIM>(GetNumNodes(), false, 0.0, 0.0);

    // Update the node location
    c_vector<double, SPACE_DIM> new_node_position = 0.5*(pNodeA->rGetLocation() + pNodeB->rGetLocation());
    ChastePoint<SPACE_DIM> point(new_node_position);    
    p_new_node->SetPoint(new_node_position);
    
    // Add node to mesh
    mNodes.push_back(p_new_node);

    // Iterate over common elements    
    for (std::set<unsigned>::iterator iter=shared_elements.begin();
         iter!=shared_elements.end();
         ++iter)
    {
        // Find which node has the lower local index in this element
        /// \todo tidy this code up (see #885)
        unsigned local_indexA = GetElement(*iter)->GetNodeLocalIndex(pNodeA->GetIndex());
        unsigned local_indexB = GetElement(*iter)->GetNodeLocalIndex(pNodeB->GetIndex());
        
        unsigned index = local_indexB;
        if ( (local_indexA == 0) || (local_indexB == 0) || (local_indexB > local_indexA) )
        {
            index = local_indexA;                 
        }
        if ( (local_indexA == 0) && (local_indexB == GetElement(*iter)->GetNumNodes()-1))
        if (local_indexA == 0)
        {
            index = local_indexB;
        }

        // Add new node to this element
        GetElement(*iter)->AddNode(index, p_new_node);
        
        // Add this element to new node
        p_new_node->AddElement(*iter);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh(NodeMap& elementMap)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert( SPACE_DIM==2 || SPACE_DIM==3 );
    assert( ELEMENT_DIM == SPACE_DIM );
    #undef COVERAGE_IGNORE

    // Make sure the map is big enough
    elementMap.Resize(GetNumElements());
    
    if (SPACE_DIM==2)
    {    
        unsigned new_index = 0;
        for (unsigned i=0; i<GetNumElements(); i++)
        {
            if (mElements[i]->IsDeleted())
            {
                elementMap.SetDeleted(i);
            }
            else
            {
                elementMap.SetNewIndex(i, new_index);
                new_index++;
            }
        }

        /*
         * We do not need to call Clear() and remove all current data, since
         * cell birth, rearrangement and death result only in local remeshing
         * of a vertex-based mesh.
         * 
         * Instead, we should now remove any deleted nodes and elements.
         * 
         * We should then construct any new nodes, including boundary nodes; 
         * then new elements; then new edges.
         * 
         * Finally (or should this be at the start?), we should perform any 
         * cell rearrangements.
         */

        // Start of element rearrangement code...
                
        // Restart check after each T1Swap as it changes elements
        bool recheck_mesh = true;
        while (recheck_mesh == true)
        {
            recheck_mesh = false;

            // Loop over elements
            for (unsigned elem_index=0; elem_index<mElements.size(); elem_index++)
            {
                if (!recheck_mesh)
                {
                    unsigned num_nodes = mElements[elem_index]->GetNumNodes();
                    assert(num_nodes>0); // if not element should be deleted
                    
                    // Loop over element vertices
                    for (unsigned local_index=0; local_index<num_nodes; local_index++)
                    {
                        // Find locations of current node and anticlockwise node
                        Node<SPACE_DIM>* p_current_node = mElements[elem_index]->GetNode(local_index);
                        unsigned local_index_plus_one = (local_index+1)%num_nodes; //TODO Should use iterators to tidy this up
                        Node<SPACE_DIM>* p_anticlockwise_node = mElements[elem_index]->GetNode(local_index_plus_one);
                        
                        // Find distance between nodes
                        double distance_between_nodes = norm_2(p_current_node->rGetLocation() - p_anticlockwise_node->rGetLocation());
        
                        if (distance_between_nodes < mThresholdDistance)
                        {
                            // Identify the type of node swap/merge needed then call method to perform swap/merge
                            IdentifySwapType(p_current_node, p_anticlockwise_node);
                            
                            recheck_mesh = true;
                            break;
                        } 
                    }
                }
                else
                {
                    break;
                }
            } 
        }        
        // ... end of element rearrangement code
        
        // areas and perimeters of elements are sorted in T1Swap Method.
    }
    else // 3D
    {
        #define COVERAGE_IGNORE
        EXCEPTION("Remeshing has not been implemented in 3D (see #827 and #860)\n");
        #undef COVERAGE_IGNORE
        /// \todo put code for remeshing in 3D here (see also the paper doi:10.1016/j.jtbi.2003.10.001)
    }    
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh()
{
    NodeMap map(GetNumElements());
    ReMesh(map);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
    // Create a set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices; 
    
    // Find the indices of the elements owned by this node
    std::set<unsigned> containing_elem_indices = GetNode(nodeIndex)->rGetContainingElementIndices();

    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter=containing_elem_indices.begin();
         elem_iter != containing_elem_indices.end();
         ++elem_iter)
    {
        // Find the local index of this node in this element
        unsigned local_index = GetElement(*elem_iter)->GetNodeLocalIndex(nodeIndex);
        
        // Find the global indices of the preceding and successive nodes in this element
        unsigned num_nodes = GetElement(*elem_iter)->GetNumNodes();
        unsigned previous_local_index = (local_index - 1)%num_nodes;
        unsigned next_local_index = (local_index + 1)%num_nodes;
        
        // Add the global indices of these two nodes to the set of neighbouring node indices
        neighbouring_node_indices.insert(GetElement(*elem_iter)->GetNodeGlobalIndex(previous_local_index));
        neighbouring_node_indices.insert(GetElement(*elem_iter)->GetNodeGlobalIndex(next_local_index));
    }
    
    return neighbouring_node_indices;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodeNotAlsoInElement(unsigned nodeIndex, unsigned elemIndex)
{
    /// \todo We should probably assert here that the node is in fact contained in the element (see #827)

    // Create a set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices_not_in_this_element;
    
    // Get the indices of this node's neighbours
    std::set<unsigned> node_neighbours = GetNeighbouringNodeIndices(nodeIndex);
    
    // Get the indices of the nodes contained in this element
    std::set<unsigned> node_indices_in_this_element;
    for (unsigned local_index=0; local_index<GetElement(elemIndex)->GetNumNodes(); local_index++)
    {
        unsigned global_index = GetElement(elemIndex)->GetNodeGlobalIndex(local_index);
        node_indices_in_this_element.insert(global_index);
    }

    // Check if each neighbour is also in this element; if not, add it to the set
    for (std::set<unsigned>::iterator iter=node_neighbours.begin();
         iter!=node_neighbours.end();
         ++iter)
    {
        if (node_indices_in_this_element.find(*iter) == node_indices_in_this_element.end())
        {
            neighbouring_node_indices_not_in_this_element.insert(*iter);
        }
    }
    
    return neighbouring_node_indices_not_in_this_element;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::IdentifySwapType(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert( SPACE_DIM==2 ); // only works in 2D at present
    assert( ELEMENT_DIM == SPACE_DIM );
    #undef COVERAGE_IGNORE

    c_vector<double, SPACE_DIM> nodeA_location = pNodeA->rGetLocation();
    c_vector<double, SPACE_DIM> nodeB_location = pNodeB->rGetLocation(); 
    
    // Find the sets of elements containing nodes A and B
    std::set<unsigned> nodeA_elem_indices = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> nodeB_elem_indices = pNodeB->rGetContainingElementIndices();    
    
    // Form the set union
    std::set<unsigned> all_indices, temp_set;
    std::set_union(nodeA_elem_indices.begin(), nodeA_elem_indices.end(), 
                   nodeB_elem_indices.begin(), nodeB_elem_indices.end(), 
                   std::inserter(temp_set, temp_set.begin()));
    all_indices.swap(temp_set); // temp_set will be deleted
    
    if (all_indices.size()==1) // nodes are only in one elment hence on boundary so merge nodes
    {
        /*
         * Looks like 
         * 
         *    A   B
         * ---o---o---
         * 
         * on the boundray of the tissue
         */
        PerformNodeMerge(pNodeA,pNodeB, all_indices);
    }
    else if (all_indices.size()==2) // nodes are in two elments hence on and interior boundary so merge nodes
    {
        if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)
        {
            /*
             * Looks like 
             * 
             *    A   B
             * ---o---o---
             * 
             * on an internal edge  
             */
             PerformNodeMerge(pNodeA,pNodeB, all_indices); 
        }
        else
        {
            /*
             * Looks like 
             *
             * Outside
             *         /
             *   --o--o (2)
             *     (1) \
             * 
             * Here we employ a PartialT1Swap 
             */
             PerformT1Swap(pNodeA, pNodeB, all_indices);
        }
    }
    else if (all_indices.size()==3) // nodes are contained in three elments 
    {
       /*
        * Looks like  
        * 
        *     A  B             A  B
        *   \                       /
        *    \  (1)           (1)  /
        * (3) o--o---   or  ---o--o (3)    Element number in brackets
        *    /  (2)           (2)  \
        *   /                       \
        *
        * Perform a PartialT1Swap 
        */
        PerformT1Swap(pNodeA, pNodeB, all_indices);
    }
    else if (all_indices.size()==4) // Correct set up for T1Swap 
    {
        /*
         * Looks like this
         * 
         *   \(1)/
         *    \ / Node A
         * (2) |   (4)     elements in Brackets
         *    / \ Node B
         *   /(3)\
         * 
         * Perform a T1Swap
         * 
         */  
        PerformT1Swap(pNodeA, pNodeB, all_indices);
    }
    else
    {
        std::cout << "\n nodes are in more than 4 elements so we can't remesh\n";
        assert(0);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformNodeMerge(Node<SPACE_DIM>* pNodeA, 
                                                          Node<SPACE_DIM>* pNodeB,
                                                          std::set<unsigned> elementsContainingNodes)
{   
    c_vector<double, SPACE_DIM> node_midpoint = 0.5*pNodeA->rGetLocation() + 0.5*pNodeB->rGetLocation();
    
    if (pNodeA->GetIndex()<pNodeB->GetIndex())
    {
        // Remove node B
        c_vector<double, SPACE_DIM>& r_nodeA_location = pNodeA->rGetModifiableLocation();
        r_nodeA_location = node_midpoint;
    
        for (std::set<unsigned>::const_iterator it = elementsContainingNodes.begin();
             it != elementsContainingNodes.end();
             ++it)
        {
            unsigned nodeB_local_index =  mElements[*it]->GetNodeLocalIndex(pNodeB->GetIndex());
            assert(nodeB_local_index < UINT_MAX); // this element should contain node B
        
            mElements[*it]->DeleteNode(nodeB_local_index); 
        }
    }
    else
    {
        // Remove node A
        c_vector<double, SPACE_DIM>& r_nodeB_location = pNodeB->rGetModifiableLocation();
        r_nodeB_location = node_midpoint;
    
        for (std::set<unsigned>::const_iterator it = elementsContainingNodes.begin();
             it != elementsContainingNodes.end();
             ++it)
        {
            unsigned nodeA_local_index =  mElements[*it]->GetNodeLocalIndex(pNodeA->GetIndex());
            assert(nodeA_local_index < UINT_MAX); // this element should contain node A
        
            mElements[*it]->DeleteNode(nodeA_local_index); 
        }
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformT1Swap(Node<SPACE_DIM>* pNodeA,
                                                       Node<SPACE_DIM>* pNodeB,
                                                       std::set<unsigned> elementsContainingNodes)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert( SPACE_DIM==2 ); // only works in 2D at present
    assert( ELEMENT_DIM == SPACE_DIM );
    #undef COVERAGE_IGNORE
    
    /*
     * Restructure elements - remember to update nodes and elements.
     * 
     * We need to implement the following changes:
     * 
     * The element whose index was in nodeA_elem_indices but not nodeB_elem_indices,
     * and the element whose index was in nodeB_elem_indices but not nodeA_elem_indices,
     * should now both contain nodes A and B. 
     * 
     * The element whose index was in nodeA_elem_indices and nodeB_elem_indices, and which 
     * node C lies inside, should now only contain node A. 
     * 
     * The element whose index was in nodeA_elem_indices and nodeB_elem_indices, and which 
     * node D lies inside, should now only contain node B.
     * 
     * Iterate over all elements involved and identify which element they are 
     * in the diagram then update the nodes as necessary.
     * 
     *   \(1)/
     *    \ / Node A
     * (2) |   (4)     elements in brackets
     *    / \ Node B
     *   /(3)\
     * 
     */  

    /*
     * Compute the locations of two new nodes C, D, placed on either side of the 
     * edge E_old formed by nodes current_node and anticlockwise_node, such 
     * that the edge E_new formed by the new nodes is the perpendicular bisector 
     * of E_old, with |E_new| 'just larger' than mThresholdDistance.
     */
          
    // Find the sets of elements containing nodes A and B
    std::set<unsigned> nodeA_elem_indices = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> nodeB_elem_indices = pNodeB->rGetContainingElementIndices();  

    double distance_between_nodes_CD = 2*mThresholdDistance; /// \todo Decide what this should really be (see #860)

    c_vector<double, SPACE_DIM> nodeA_location = pNodeA->rGetLocation();
    c_vector<double, SPACE_DIM> nodeB_location = pNodeB->rGetLocation();

    c_vector<double, SPACE_DIM> a_to_b = nodeB_location - nodeA_location;    
    c_vector<double, SPACE_DIM> perpendicular_vector;
    perpendicular_vector(0) = -a_to_b(1);
    perpendicular_vector(1) = a_to_b(0);

    c_vector<double, SPACE_DIM> c_to_d = distance_between_nodes_CD / norm_2(nodeB_location - nodeA_location) * perpendicular_vector;    
    c_vector<double, SPACE_DIM> nodeC_location = nodeA_location + 0.5*a_to_b - 0.5*c_to_d;
    c_vector<double, SPACE_DIM> nodeD_location = nodeC_location + c_to_d;

    /*
     * Move node A to C and node B to D
     */
     
    c_vector<double, SPACE_DIM>& r_nodeA_location = pNodeA->rGetModifiableLocation();
    r_nodeA_location = nodeC_location;
    
    c_vector<double, SPACE_DIM>& r_nodeB_location = pNodeB->rGetModifiableLocation();
    r_nodeB_location = nodeD_location;
    
    for (std::set<unsigned>::const_iterator it = elementsContainingNodes.begin();
         it != elementsContainingNodes.end();
         ++it)
    {
        if (nodeA_elem_indices.find(*it) == nodeA_elem_indices.end()) // not in nodeA_elem_indices so element 3
        {
            /*
             * In this case the element index was not in 
             * nodeA_elem_indices, so this element 
             * does not contain node A. Therefore we must add node A
             * (which has been moved to node C) to this element.
             *
             * Locate local index of node B in element then add node A after 
             * in anticlockwise direction. 
             */  
            
            unsigned nodeB_local_index =  mElements[*it]->GetNodeLocalIndex(pNodeB->GetIndex());
            assert(nodeB_local_index < UINT_MAX); // this element should contain node B
    
            mElements[*it]->AddNode(nodeB_local_index,pNodeA);
        }
        else if (nodeB_elem_indices.find(*it) == nodeB_elem_indices.end()) // not in nodeB_elem_indices so element 1
        {
            /*
             * In this case the element index was not in 
             * nodeB_elem_indices, so this element
             * does not contain node B. Therefore we must add node B
             * (which has been moved to node D) to this element.
             *
             * Locate local index of node A in element then add node B after 
             * in anticlockwise direction. 
             */  
            unsigned nodeA_local_index =  mElements[*it]->GetNodeLocalIndex(pNodeA->GetIndex());
            assert(nodeA_local_index < UINT_MAX); // this element should contain node A
            mElements[*it]->AddNode(nodeA_local_index,pNodeB); 
        }    
        else
        {
            /*
             * In this case the element index was in both nodeB_elem_indices and nodeB_elem_indices
             * so is element 2 or 4 
             */
            
            /*
             * Locate local index of nodeA and nodeB and use the oredering to 
             * identify the element if nodeB_index > nodeA_index then element 4
             * and if nodeA_index > nodeB_index then element 2 
             */  
            unsigned nodeA_local_index =  mElements[*it]->GetNodeLocalIndex(pNodeA->GetIndex());
            assert(nodeA_local_index < UINT_MAX); // this element should contain node A
            
            unsigned nodeB_local_index =  mElements[*it]->GetNodeLocalIndex(pNodeB->GetIndex());
            assert(nodeB_local_index < UINT_MAX); // this element should contain node B
     
            unsigned nodeB_local_index_plus_one = (nodeB_local_index + 1)%(mElements[*it]->GetNumNodes());
            
            if (nodeA_local_index == nodeB_local_index_plus_one)
            {
                /*
                 * In this case the local index of nodeA is the local index of 
                 * nodeB plus one so we are in element 2 so we remove nodeB
                 */
                 mElements[*it]->DeleteNode(nodeB_local_index); 
            }
            else
            {
                assert(nodeB_local_index == (nodeA_local_index + 1)%(mElements[*it]->GetNumNodes())); // as A and B are next to each other
                /*
                 * In this case the local index of nodeA is the local index of 
                 * nodeB minus one so we are in element 4 so we remove nodeA
                 */             
                 mElements[*it]->DeleteNode(nodeA_local_index); 
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert( SPACE_DIM==2 ); // only works in 2D at present
    assert( ELEMENT_DIM == SPACE_DIM );
    #undef COVERAGE_IGNORE
    
    // Find short axis
    c_vector<double, SPACE_DIM> centroid = pElement->CalculateCentroid();
    c_vector<double, SPACE_DIM> short_axis = pElement->CalculateShortAxis();
    c_vector<double, SPACE_DIM> long_axis; // this is perpendicular to the short axis
    long_axis(0) = -short_axis(1);
    long_axis(1) = short_axis(0);
    
    /// \todo Remove this temporary bool
    
    unsigned num_nodes = pElement->GetNumNodes();

    // Store if the node is on the side of the short axis which the long axis points to 
    bool is_on_left[num_nodes];
    
    for (unsigned i=0; i<num_nodes; i++)
    {
        c_vector<double, SPACE_DIM> node_location_from_centroid = pElement->GetNodeLocation(i)- centroid; 
        
        if (inner_prod(node_location_from_centroid, long_axis) >= 0)
        {
            is_on_left[i] = true;
        }
        else // inner_prod(node_location_from_centroid,long_axis)<0
        {
            is_on_left[i] = false;
        }
    }
    
    std::vector<unsigned> intersecting_nodes;
        
    for (unsigned i=0; i<num_nodes-1; i++)
    {
        if (is_on_left[i]!=is_on_left[i+1])
        {
            intersecting_nodes.push_back(i);
        }
    }
    if (is_on_left[0]!=is_on_left[num_nodes-1])
    {
        intersecting_nodes.push_back(num_nodes-1);
    }
    
    /// \todo remove assert and make an if statement returning an error if needed
    assert(intersecting_nodes.size()==2); // only divide if 2 intersections 
    
    std::vector<unsigned> new_node_global_indices;  
    
    for (unsigned i=0; i<intersecting_nodes.size(); i++)
    {
        // Find intersections between edges and short_axis
        c_vector<double, SPACE_DIM> position_a = pElement->GetNodeLocation(intersecting_nodes[i]);
        c_vector<double, SPACE_DIM> position_b;
        if (intersecting_nodes[i] < num_nodes-1)
        {
            position_b = pElement->GetNodeLocation(intersecting_nodes[i]+1);
        }
        else
        {
            position_b = pElement->GetNodeLocation(0);
        }
        
        c_vector<double, SPACE_DIM> a_to_b = position_b - position_a;
        
        /*
         * Let the first one on edge be a and the second one be b, 
         * then we are interested in the intersection of
         *  
         * position_a + alpha * a_to_b 
         * and
         * centroid + beta * short_axis 
         * 
         */

        double determinant = a_to_b[0]*short_axis[1] - a_to_b[1]*short_axis[0];
         
        double alpha = (centroid[0]*a_to_b[1]-position_a[0]*a_to_b[1]
                        -centroid[1]*a_to_b[0]+position_a[1]*a_to_b[0])/determinant;

        c_vector<double, SPACE_DIM> intersection = centroid + alpha*short_axis;
        
        // Create New node with location intersection and add it to all correct elements
        unsigned node_global_index = this->AddNode(new Node<SPACE_DIM>(0, false, intersection[0], intersection[1]));
        new_node_global_indices.push_back(node_global_index);
        
        /*
         * 
         * 
         *   Now need to add node to correct elements
         * 
         * 
         */
    }
    
//    std::cout << "\n global indices" << new_node_global_indices.size() << "\t" << new_node_global_indices[0] << "\t" << new_node_global_indices[1] << std::flush; 
//    std::cout << "\n element nodes \t" << pElement->GetNumNodes() << "\t" << std::flush; 
//    // NOTE pElement->GetNodeLocalIndex will return UINT_MAX if node is not in this element.
//    std::cout << "\n local indices" << new_node_global_indices.size() << "\t" << pElement->GetNodeLocalIndex(new_node_global_indices[0]) << "\t" << pElement->GetNodeLocalIndex(new_node_global_indices[1]) << std::flush; 
//        
//    // Now call DivideElement(.,.,.) to divide the elemenent using the new nodes 
//    DivideElement(pElement, pElement->GetNodeLocalIndex(new_node_global_indices[0]),pElement->GetNodeLocalIndex(new_node_global_indices[1]));
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned NodeAIndex, unsigned NodeBIndex)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert( SPACE_DIM==2 ); // only works in 2D at present
    assert( ELEMENT_DIM == SPACE_DIM );
    #undef COVERAGE_IGNORE
   
    assert(NodeBIndex!=NodeAIndex);
    // sort NodeA and NodeB such that NodeBIndex>NodeAindex
    unsigned Node1Index, Node2Index;
    if (NodeAIndex<NodeBIndex)
    {
        Node1Index=NodeAIndex;
        Node2Index=NodeBIndex;
    }
    else
    {
        Node1Index=NodeBIndex;
        Node2Index=NodeAIndex;
    }
    
//    std::cout << "\n node1 " << Node1Index << std::flush;
//    std::cout << "\n node2 " << Node2Index << std::flush;
    
    // copy element 
//    std::cout<<"\nCopy Element \n"<< std::flush; 
    std::vector<Node<SPACE_DIM>*> nodes_elem;
    unsigned num_nodes = pElement->GetNumNodes(); //Store this as it changes when you delete nodes from element
    
    for(unsigned i=0; i<num_nodes; i++)
    {
        nodes_elem.push_back(pElement->GetNode(i));
    }
    unsigned new_element_index = AddElement(new VertexElement<ELEMENT_DIM,SPACE_DIM>(0, nodes_elem));
    
    // Remove nodes  # < node1 and # > node2 from pElement  
    // Remove nodes node1 < # < node2 from new_element
    
//    std::cout << "\n "<<  mElements[new_element_index]->GetNumNodes() << "\n\n" << std::flush;
//    for(unsigned i=0; i<mElements[new_element_index]->GetNumNodes(); i++)
//    {   
//        std::cout << "\n containing elements of node " <<  mElements[new_element_index]->GetNodeGlobalIndex(i) << "\t ("<<  mElements[new_element_index]->GetNode(i)->GetNumContainingElements() << ")" <<std::flush;
//    }
    
    for(unsigned i=num_nodes; i>0; i--)
    {
        
        if (i-1<Node1Index || i-1>Node2Index)
        {
            pElement->DeleteNode(i-1);
        }
        else if (i-1>Node1Index && i-1<Node2Index)
        {
            mElements[new_element_index]->DeleteNode(i-1);
        }
    }
    return new_element_index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(VertexMeshReader2d& rMeshReader)
{
    // Store numbers of nodes and elements
    unsigned num_nodes = rMeshReader.GetNumNodes();
    unsigned num_elements = rMeshReader.GetNumElements();
    
    // Reserve memory for nodes
    mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    // Add nodes
    std::vector<double> coords;
    for (unsigned i=0; i<num_nodes; i++)
    {
        coords = rMeshReader.GetNextNode();
        mNodes.push_back(new Node<SPACE_DIM>(i, coords, false));
    }

    rMeshReader.Reset();

    // Reserve memory for nodes
    mElements.reserve(rMeshReader.GetNumElements());

    // Add elements
    for (unsigned elem_index=0; elem_index<num_elements; elem_index++)
    {
        VertexElementData element_data = rMeshReader.GetNextElementData();
        std::vector<Node<SPACE_DIM>*> nodes;

        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j=0; j<num_nodes_in_element; j++)
        {
            assert(element_data.NodeIndices[j] < mNodes.size());
            nodes.push_back(mNodes[element_data.NodeIndices[j]]);
        }

        VertexElement<ELEMENT_DIM,SPACE_DIM>* p_element = new VertexElement<ELEMENT_DIM,SPACE_DIM>(elem_index, nodes);
        mElements.push_back(p_element);
    
        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = element_data.AttributeValue;
            p_element->SetRegion(attribute_value);
        }        
    }
}


/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class VertexMesh<1,1>;
template class VertexMesh<1,2>;
template class VertexMesh<2,2>;
template class VertexMesh<2,3>;
template class VertexMesh<3,3>;
