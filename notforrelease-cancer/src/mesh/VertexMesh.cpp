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


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
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
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point)
{
    mNodes[nodeIndex]->SetPoint(point);
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
        /// \todo put code for remeshing in 2D here (see #827)
    
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
        
        // Loop over elements
        for (unsigned elem_index=0; elem_index<mElements.size(); elem_index++)
        {
            unsigned num_nodes = mElements[elem_index]->GetNumNodes();

            // Loop over element vertices
            for (unsigned local_index=0; local_index<num_nodes; local_index++)
            {
                // Find locations of current node and anticlockwise node
                Node<SPACE_DIM>* p_current_node = mElements[elem_index]->GetNode(local_index);
                Node<SPACE_DIM>* p_anticlockwise_node = mElements[elem_index]->GetNode(local_index);
                
                // Find distance between nodes
                double distance_between_nodes = norm_2(p_current_node->rGetLocation() - p_anticlockwise_node->rGetLocation());

                if (distance_between_nodes < mThresholdDistance)
                {
                    PerformT1Swap(p_current_node, p_anticlockwise_node);
                } 
            }
        }        
        // ... end of element rearrangement code
    }
    else // 3D
    {
        #define COVERAGE_IGNORE
        EXCEPTION("Remeshing has not tes been implemented in 3D (see #827 and #860)\n");
        #undef COVERAGE_IGNORE
        /// \todo put code for remeshing in 3D here (see also the paper doi:10.1016/j.jtbi.2003.10.001)
    }    
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh()
{
    NodeMap map(GetNumElements());
    ReMesh(map);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformT1Swap(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert( SPACE_DIM==2 ); // only works in 2D at present
    assert( ELEMENT_DIM == SPACE_DIM );
    #undef COVERAGE_IGNORE

    /*
     * Compute the locations of two new nodes C, D, placed on either side of the 
     * edge E_old formed by nodes current_node and anticlockwise_node, such 
     * that the edge E_new formed by the new nodes is the perpendicular bisector 
     * of E_old, with |E_new| 'just larger' than mThresholdDistance.
     */
    c_vector<double, SPACE_DIM> nodeA_location = pNodeA->rGetLocation();
    c_vector<double, SPACE_DIM> nodeB_location = pNodeB->rGetLocation(); 
    
    double distance_between_nodes_AB = norm_2(nodeB_location - nodeA_location);
    double distance_between_nodes_CD = 2*mThresholdDistance; /// \todo Decide what this should really be (see #860)
    
    double hypotenuse = 0.5*sqrt(pow(distance_between_nodes_CD, 2.0) + pow(distance_between_nodes_AB, 2.0));
    
    double angle_ABC = atan(distance_between_nodes_CD/distance_between_nodes_AB);
    double angle_between_AB_and_horizontal = atan2(nodeB_location[1]-nodeA_location[1], 
                                                   nodeB_location[0]-nodeA_location[0]);
    
    double x_offset_from_nodeB = hypotenuse*sin(angle_ABC + angle_between_AB_and_horizontal);
    double y_offset_from_nodeB = hypotenuse*cos(angle_ABC + angle_between_AB_and_horizontal);
    
    c_vector<double, SPACE_DIM> nodeC_location;
    nodeC_location[0] = nodeB_location[0] - x_offset_from_nodeB;
    nodeC_location[1] = nodeB_location[1] - y_offset_from_nodeB;
    
    c_vector<double, SPACE_DIM> nodeD_location = nodeA_location + nodeB_location - nodeC_location;
     
    /*
     * Move node A to C and node B to D
     */
    c_vector<double, SPACE_DIM>& r_nodeA_location = pNodeA->rGetModifiableLocation();
    r_nodeA_location = nodeC_location;
    
    c_vector<double, SPACE_DIM>& r_nodeB_location = pNodeA->rGetModifiableLocation();
    r_nodeB_location = nodeD_location;

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
     */

    /*
     * Find elements containing nodes A and B
     * 
     * \todo try using set_union instead
     */
    std::set<unsigned> nodeA_elem_indices = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> nodeB_elem_indices = pNodeB->rGetContainingElementIndices();    
    std::set<unsigned> neighbouring_element_indices;    
    std::pair<std::set<unsigned>::iterator, bool> pr;
    
    for (std::set<unsigned>::const_iterator it = nodeA_elem_indices.begin();
         it != nodeA_elem_indices.end();
         ++it)
    {
        neighbouring_element_indices.insert(*it);
    }
    for (std::set<unsigned>::const_iterator it = nodeB_elem_indices.begin();
         it != nodeB_elem_indices.end();
         ++it)
    {
        pr = neighbouring_element_indices.insert(*it);
        if (pr.second)
        {
            /*
             * In this case the element index was not already in 
             * the set neighbouring_element_indices, so this element 
             * does not contain node A. Therefore we must add node A
             * (which has been moved to node C) to this element.
             */ 
        }
        else
        {
            /*
             * In this case the element index was already in the 
             * set neighbouring_element_indices, so this element 
             * also contains node A. Therefore we must either node A
             * (which has been moved to node C) or node B (which has
             * been moved to node D) from this element.
             */
        }
    }
    /*
     * The final element to deal with is the one which contains node A
     * but not node B. We must add node B (which has been moved to node D)
     * to this element.
     */
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
    for (unsigned i=0; i < num_nodes; i++)
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
            assert(element_data.NodeIndices[j] <  mNodes.size());
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
