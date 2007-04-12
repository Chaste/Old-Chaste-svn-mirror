#ifndef _CONFORMINGTETRAHEDRALMESH_CPP_
#define _CONFORMINGTETRAHEDRALMESH_CPP_

#include "ConformingTetrahedralMesh.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConformingTetrahedralMesh()
{}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConformingTetrahedralMesh(unsigned numElements)
{
    mElements.reserve(numElements);
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
    unsigned orderOfBasisFunctions)
{

    // We only use linear or quadratic basis functions
    assert(orderOfBasisFunctions == 1 || orderOfBasisFunctions == 2);
    
    // Record number of corner nodes
    mNumCornerNodes = rMeshReader.GetNumNodes();
    
    // Reserve memory for nodes, so we don't have problems with pointers stored in
    // elements becoming invalid.
    // For quadratic basis functions we need to go through and count them first
    
    if (orderOfBasisFunctions == 1)
    {
        mNodes.reserve(mNumCornerNodes);
    }
    else if (orderOfBasisFunctions == 2)
    {
        typename std::map<std::pair<unsigned,unsigned>,unsigned>::const_iterator iterator;
        std::map<std::pair<unsigned,unsigned>,unsigned> internal_nodes_map;
        
        std::vector<Node<SPACE_DIM> > temp_nodes;
        temp_nodes.reserve(mNumCornerNodes);
        
        // Add corner nodes to temporary nodes vector
        std::vector<double> coords;
        for (unsigned i=0; i < mNumCornerNodes; i++)
        {
            coords = rMeshReader.GetNextNode();
            temp_nodes.push_back(Node<SPACE_DIM>(i, coords, false));
        }
        
        unsigned new_node_index = mNumCornerNodes;
        
        // Create the map of internal node indices
        for (unsigned i=0; i < rMeshReader.GetNumElements(); i++)
        {
            std::vector<unsigned> node_indices = rMeshReader.GetNextElement();
            std::vector<const Node<SPACE_DIM>*> nodes;
            
            //Debugging code
            unsigned nodes_size = node_indices.size();
            
            for (unsigned j=0; j < nodes_size; j++)
            {
                assert(node_indices[j] < temp_nodes.size());
                nodes.push_back(&temp_nodes[node_indices[j]]);
            }
            
            for (unsigned j=0; j < ELEMENT_DIM + 1; j++)
            {
                for (unsigned k=j+1; k < ELEMENT_DIM + 1; k++)
                {
                    assert(j < nodes_size && k < nodes_size);
                    unsigned node_i = nodes[j]->GetIndex();
                    unsigned node_j = nodes[k]->GetIndex();
                    if (node_j < node_i)
                    {
                        unsigned temp = node_i;
                        node_i = node_j;
                        node_j = temp;
                    }
                    iterator = internal_nodes_map.find(std::pair<unsigned,unsigned>(node_i, node_j));
                    if (iterator == internal_nodes_map.end())
                    {
                        // add node to map
                        internal_nodes_map[(std::pair<unsigned,unsigned>(node_i, node_j))] = new_node_index;
                        new_node_index++;
                    }
                }
            }
        }
        mNodes.reserve(new_node_index);
    }
    
    rMeshReader.Reset();
    
    typename std::map<std::pair<unsigned,unsigned>,unsigned>::const_iterator iterator;
    std::map<std::pair<unsigned,unsigned>,unsigned> internal_nodes_map;
    
    // Add corner nodes
    std::vector<double> coords;
    for (unsigned i=0; i < mNumCornerNodes; i++)
    {
        coords = rMeshReader.GetNextNode();
        mNodes.push_back(new Node<SPACE_DIM>(i, coords, false));
    }
    
    unsigned new_node_index = mNumCornerNodes;
    
    if (orderOfBasisFunctions == 2)
    {
        for (unsigned i=0; i < rMeshReader.GetNumElements(); i++)
        {
            std::vector<unsigned> node_indices = rMeshReader.GetNextElement();
            std::vector<const Node<SPACE_DIM>*> nodes;
            
            //Debugging code
            unsigned nodes_size = node_indices.size();
            
            for (unsigned j=0; j < nodes_size; j++)
            {
                assert(node_indices[j] < mNodes.size());
                nodes.push_back(mNodes[node_indices[j]]);
            }
            
            for (unsigned j=0; j < ELEMENT_DIM + 1; j++)
            {
                for (unsigned k=j+1; k < ELEMENT_DIM + 1; k++)
                {
                    assert(j < nodes_size && k < nodes_size);
                    unsigned node_i = nodes[j]->GetIndex();
                    unsigned node_j = nodes[k]->GetIndex();
                    if (node_j < node_i)
                    {
                        unsigned temp = node_i;
                        node_i = node_j;
                        node_j = temp;
                    }
                    iterator = internal_nodes_map.find(std::pair<unsigned,unsigned>(node_i, node_j));
                    if (iterator == internal_nodes_map.end())
                    {
                        // add node to map
                        internal_nodes_map[(std::pair<unsigned,unsigned>(node_i, node_j))] = new_node_index;
                        // add node to mesh
                        const Node<SPACE_DIM>* p_node1 = GetNode(node_i);
                        const Node<SPACE_DIM>* p_node2 = GetNode(node_j);
                        c_vector<double,SPACE_DIM> mid_point
                        =  0.5 * (p_node1->rGetLocation() + p_node2->rGetLocation());
                        
                        Node<SPACE_DIM> *p_new_node=new Node<SPACE_DIM>(new_node_index,
                                                                        mid_point,
                                                                        p_node1->IsBoundaryNode() && p_node2->IsBoundaryNode());
                        mNodes.push_back(p_new_node);
                        new_node_index++;
                    }
                }
            }
        }
    }
    rMeshReader.Reset();
    // Add elements
    new_node_index = mNumCornerNodes;
    mElements.reserve(rMeshReader.GetNumElements());
    
    for (unsigned element_index=0; element_index < (unsigned) rMeshReader.GetNumElements(); element_index++)
    {
        std::vector<unsigned> node_indices = rMeshReader.GetNextElement();
        std::vector<Node<SPACE_DIM>*> nodes;
        unsigned nodes_size = node_indices.size();
        for (unsigned j=0; j<nodes_size; j++)
        {
            assert(node_indices[j] <  mNodes.size());
            nodes.push_back(mNodes[node_indices[j]]);
        }
        
        if (orderOfBasisFunctions == 2)
        {
            for (unsigned j=0; j < ELEMENT_DIM + 1; j++)
            {
                for (unsigned k=j+1; k < ELEMENT_DIM + 1; k++)
                {
                    assert(j < nodes_size && k < nodes_size);
                    unsigned node_i = nodes[j]->GetIndex();
                    unsigned node_j = nodes[k]->GetIndex();
                    if (node_j < node_i)
                    {
                        unsigned temp = node_i;
                        node_i = node_j;
                        node_j = temp;
                    }
                    iterator = internal_nodes_map.find(std::pair<unsigned,unsigned>(node_i, node_j));
                    assert(iterator != internal_nodes_map.end());
                    // add node to element
                    nodes.push_back(this->GetNode(iterator->second));
                    new_node_index++;
                }
            }
        }
        mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(element_index, nodes,orderOfBasisFunctions));
    }
    
    // Add boundary elements & nodes
    unsigned actual_face_index=0;
    for (unsigned face_index=0; face_index<(unsigned)rMeshReader.GetNumFaces(); face_index++)
    {
        std::vector<unsigned> node_indices = rMeshReader.GetNextFace();
        
        // Determine if this is a boundary face
        std::set<unsigned> containing_element_indices; // Elements that contain this face
        std::vector<Node<SPACE_DIM>*> nodes;
        for (unsigned node_index=0; node_index<node_indices.size(); node_index++)
        {
            assert(node_indices[node_index] <  mNodes.size());
            // Add Node pointer to list for creating an element
            nodes.push_back(mNodes[node_indices[node_index]]);
            
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
        
        //If the following assertion is thrown, it means that the .edge/.face file does not
        //match the .ele file -- they were generated at separate times.  Simply remove the internal
        //edges/faces by hand.
        assert(containing_element_indices.size() != 0);
        
        if (containing_element_indices.size() == 1)
        {
            // This is a boundary face
            // Ensure all its nodes are marked as boundary nodes
            for (unsigned j=0; j<nodes.size(); j++)
            {
                if (!nodes[j]->IsBoundaryNode())
                {
                    nodes[j]->SetAsBoundaryNode();
                    mBoundaryNodes.push_back(nodes[j]);
                }
                //Register the index that this bounday element will have
                //with the node
                nodes[j]->AddBoundaryElement(actual_face_index);
            }
            
            
            if (orderOfBasisFunctions == 2)
            {
                unsigned nodes_size = node_indices.size();
                UNUSED_OPT(nodes_size);
                for (unsigned j=0; j < ELEMENT_DIM; j++)
                {
                    for (unsigned k=j+1; k < ELEMENT_DIM; k++)
                    {
                        assert(j < nodes_size && k < nodes_size);
                        unsigned node_i = nodes[j]->GetIndex();
                        unsigned node_j = nodes[k]->GetIndex();
                        if (node_j < node_i)
                        {
                            unsigned temp = node_i;
                            node_i = node_j;
                            node_j = temp;
                        }
                        iterator = internal_nodes_map.find(std::pair<unsigned,unsigned>(node_i, node_j));
                        assert(iterator != internal_nodes_map.end());
                        // add node to element
                        nodes.push_back(this->GetNode(iterator->second));
                    }
                }
            }
            
            // The added elements will be deleted in our destructor
            mBoundaryElements.push_back(
                new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(actual_face_index,nodes,orderOfBasisFunctions));
            actual_face_index++;
            
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~ConformingTetrahedralMesh()
{
    // Iterate over nodes and free the memory
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        delete mNodes[i];
    }
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
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM> *pNewNode)
{

    if (mDeletedNodeIndices.empty())
    {
        pNewNode->SetIndex(mNodes.size());
        mNodes.push_back(pNewNode);
    }
    else
    {
        unsigned index=mDeletedNodeIndices.back();
        pNewNode->SetIndex(index);
        mDeletedNodeIndices.pop_back();
        delete mNodes[index];
        mNodes[index] = pNewNode;
    }
    
    return pNewNode->GetIndex();
}


/**
 * Get a node reference from the mesh.
 *
 * Note that this may become invalid if nodes are subsequently added to the mesh.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM> *ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index)
{
    assert(index < mNodes.size());
    return (mNodes[index]);
}

/// Returns the number of nodes that are actually in use
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes()
{
    return mNodes.size() - mDeletedNodeIndices.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllNodes()
{
    return mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements()
{
    return mElements.size() - mDeletedElementIndices.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements()
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryNodes()
{
    return mBoundaryNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements()
{
    return mBoundaryElements.size() - mDeletedBoundaryElementIndices.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllBoundaryElements()
{
    return mBoundaryElements.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumCornerNodes()
{
    return mNumCornerNodes;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RescaleMeshFromBoundaryNode(Point<1> updatedPoint, unsigned boundaryNodeIndex)
{
    assert(GetNode(boundaryNodeIndex)->IsBoundaryNode());
    double scaleFactor = updatedPoint[0] / GetNode(boundaryNodeIndex)->GetPoint()[0];
    double temp;
    for (unsigned i=0; i < boundaryNodeIndex+1; i++)
    {
        temp = scaleFactor * mNodes[i]->GetPoint()[0];
        Point<1> newPoint(temp);
        mNodes[i]->SetPoint(newPoint);
    }
    RefreshMesh();
}

/** SetNode moves the node with a particular index to a new point in space and
  * verifies that the signed areas of the supporting Elements are positive
  * @param index is the index of the node to be moved
  * @param point is the new target location of the node
  * @param concreteMove is set to false if we want to skip the signed area tests
  *
  */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetNode(unsigned index,
        Point<SPACE_DIM> point,
        bool concreteMove)
{
    mNodes[index]->SetPoint(point);
    if (concreteMove)
    {
        for (unsigned i=0; i<mNodes[index]->GetNumContainingElements(); i++)
        {
            try
            {
                GetElement(mNodes[index]->GetNextContainingElementIndex())->RefreshJacobianDeterminant();
            }
            catch (Exception e)
            {
                if (ELEMENT_DIM == SPACE_DIM)
                {
                    EXCEPTION("Moving node caused an element to have a non-positive Jacobian determinant");
                }
                else
                {
                    EXCEPTION("Moving node caused an subspace element to change direction");
                }
            }
        }
        for (unsigned i=0; i<mNodes[index]->GetNumBoundaryElements(); i++)
        {
            try
            {
                GetBoundaryElement(mNodes[index]->GetNextBoundaryElementIndex())->RefreshJacobianDeterminant();
            }
            catch (Exception e)
            {
                EXCEPTION("Moving node caused a boundary element to have a non-positive Jacobian determinant");
            }
        }
    }
}

/**
 * DeleteNode deletes a node from the mesh by finding an appropriate neighbour node
 * to merge it with.
 *
 * @param index is the index of the node to be deleted
 *
**/
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::DeleteNode(unsigned index)
{
    unsigned target_index;
    bool found_target=false;
    if (mNodes[index]->IsDeleted())
    {
        EXCEPTION("Trying to delete a deleted node");
    }
    while (!found_target)
    {
        Element <ELEMENT_DIM,SPACE_DIM> *p_element=
            mElements[mNodes[index]->GetNextContainingElementIndex()];
        for (unsigned i=0; i<=ELEMENT_DIM && !found_target; i++)
        {
            target_index=p_element->GetNodeGlobalIndex(i);
            try
            {
                MoveMergeNode(index, target_index, false);
                found_target=true;
            }
            catch (Exception e)
            {
                //Just go round the loops and try again
            }
        }
    }
    
    MoveMergeNode(index, target_index);
}
/**
 * MoveMergeNode moves one node to another (i.e. merges the nodes), refreshing/deleting elements as
 * appropriate.
 *
 * @param index is the index of the node to be moved
 * @param targetIndex is the index of the node to move to
 * @param concreteMove can be set to false if you just want to check whether this will work.
 *     Set it to true if you're doing the merger for real, in order to do all the bookkeeping.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MoveMergeNode(unsigned index,
        unsigned targetIndex,
        bool concreteMove)
{

    if (mNodes[index]->IsDeleted())
    {
        EXCEPTION("Trying to move a deleted node");
    }
    
    if (index == targetIndex)
    {
        EXCEPTION("Trying to merge a node with itself");
    }
    if (mNodes[index]->IsBoundaryNode())
    {
        if (!mNodes[targetIndex]->IsBoundaryNode())
        {
            EXCEPTION("A boundary node can only be moved on to another boundary node");
        }
    }
    std::set<unsigned> unshared_element_indices;
    std::set_difference(mNodes[index]->rGetContainingElementIndices().begin(),
                        mNodes[index]->rGetContainingElementIndices().end(),
                        mNodes[targetIndex]->rGetContainingElementIndices().begin(),
                        mNodes[targetIndex]->rGetContainingElementIndices().end(),
                        std::inserter(unshared_element_indices, unshared_element_indices.begin()));
                        
                        
    if (unshared_element_indices.size() == mNodes[index]->rGetContainingElementIndices().size())
    {
        EXCEPTION("These nodes cannot be merged since they are not neighbours");
    }
    
    std::set<unsigned> unshared_boundary_element_indices;
    std::set_difference(mNodes[index]->rGetContainingBoundaryElementIndices().begin(),
                        mNodes[index]->rGetContainingBoundaryElementIndices().end(),
                        mNodes[targetIndex]->rGetContainingBoundaryElementIndices().begin(),
                        mNodes[targetIndex]->rGetContainingBoundaryElementIndices().end(),
                        std::inserter(unshared_boundary_element_indices, unshared_boundary_element_indices.begin()));
                        
                        
    if (mNodes[index]->IsBoundaryNode())
    {
        if (unshared_boundary_element_indices.size()
            == mNodes[index]->rGetContainingBoundaryElementIndices().size())
        {
            EXCEPTION("These nodes cannot be merged since they are not neighbours on the boundary");
        }
    }
    
    mNodes[index]->rGetModifiableLocation() = mNodes[targetIndex]->rGetLocation();
    
    for (std::set<unsigned>::const_iterator element_iter=unshared_element_indices.begin();
             element_iter != unshared_element_indices.end();
             element_iter++)
        {
            try
            {
            
                GetElement(*element_iter)->RefreshJacobianDeterminant(concreteMove);
                if (concreteMove)
                {
                    GetElement(*element_iter)->ReplaceNode(mNodes[index], mNodes[targetIndex]);
                }
                
            }
            catch (Exception e)
            {
                EXCEPTION("Moving node caused an element to have a non-positive Jacobian determinant");
            }
        }
    for (std::set<unsigned>::const_iterator boundary_element_iter=
                 unshared_boundary_element_indices.begin();
             boundary_element_iter != unshared_boundary_element_indices.end();
             boundary_element_iter++)
        {
        
            GetBoundaryElement(*boundary_element_iter)->RefreshJacobianDeterminant(concreteMove);
            if (concreteMove)
            {
                GetBoundaryElement(*boundary_element_iter)->ReplaceNode(mNodes[index], mNodes[targetIndex]);
            }
        }
        
    std::set<unsigned> shared_element_indices;
    std::set_intersection(mNodes[index]->rGetContainingElementIndices().begin(),
                          mNodes[index]->rGetContainingElementIndices().end(),
                          mNodes[targetIndex]->rGetContainingElementIndices().begin(),
                          mNodes[targetIndex]->rGetContainingElementIndices().end(),
                          std::inserter(shared_element_indices, shared_element_indices.begin()));
    for (std::set<unsigned>::const_iterator element_iter=shared_element_indices.begin();
             element_iter != shared_element_indices.end();
             element_iter++)
        {
            if (concreteMove)
            {
                GetElement(*element_iter)->MarkAsDeleted();
                mDeletedElementIndices.push_back(*element_iter);
            }
            else
            {
                GetElement(*element_iter)->ZeroJacobianDeterminant();
            }
        }
        
        
    std::set<unsigned> shared_boundary_element_indices;
    std::set_intersection(mNodes[index]->rGetContainingBoundaryElementIndices().begin(),
                          mNodes[index]->rGetContainingBoundaryElementIndices().end(),
                          mNodes[targetIndex]->rGetContainingBoundaryElementIndices().begin(),
                          mNodes[targetIndex]->rGetContainingBoundaryElementIndices().end(),
                          std::inserter(shared_boundary_element_indices, shared_boundary_element_indices.begin()));
    for (std::set<unsigned>::const_iterator boundary_element_iter=shared_boundary_element_indices.begin();
             boundary_element_iter != shared_boundary_element_indices.end();
             boundary_element_iter++)
        {
            if (concreteMove)
            {
                GetBoundaryElement(*boundary_element_iter)->MarkAsDeleted();
                mDeletedBoundaryElementIndices.push_back(*boundary_element_iter);
            }
            else
            {
                GetBoundaryElement(*boundary_element_iter)->ZeroJacobianDeterminant();
                GetBoundaryElement(*boundary_element_iter)->ZeroWeightedDirection();
            }
        }
        
    if (concreteMove)
    {
        mNodes[index]->MarkAsDeleted();
        mDeletedNodeIndices.push_back(index);
    }
}


/**
 * This method allows the mesh properties to be re-calculated after one
 * or more node have been moved.
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefreshMesh()
{
    for (unsigned i=0; i<mElements.size();i++)
    {
        if (!mElements[i]->IsDeleted())
        {
            mElements[i]->RefreshJacobianDeterminant();
        }
    }
    
    //Refresh each boundary element
    for (unsigned i=0; i<mBoundaryElements.size();i++)
    {
        if (!mBoundaryElements[i]->IsDeleted())
        {
            try
            {
                mBoundaryElements[i]->RefreshJacobianDeterminant();
            }
            catch (Exception e)
            {
                //Since we may have rotated the mesh, it's okay for normals to swing round
            }
        }
    }
    
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefineElement(
    Element<ELEMENT_DIM,SPACE_DIM>* pElement,
    Point<SPACE_DIM> point)
{

    //Check that the point is in the element
    if (pElement->IncludesPoint(point, true) == false)
    {
        EXCEPTION("RefineElement could not be started (point is not in element)");
    }
    
    // Add a new node from the point that is passed to RefineElement
    unsigned new_node_index = AddNode(new Node<SPACE_DIM>(0, point.rGetLocation()));
    // Note: the first argument is the index of the node, which is going to be
    //       overriden by AddNode, so it can safely be ignored
    
    //This loop constructs the extra elements which are going to fill the space
    for (unsigned i = 0; i < ELEMENT_DIM; i++)
    {
    
        // First, make a copy of the current element making sure we update its index
        unsigned new_elt_index;
        if (mDeletedElementIndices.empty())
        {
            new_elt_index = mElements.size();
        }
        else
        {
            new_elt_index = mDeletedElementIndices.back();
            mDeletedElementIndices.pop_back();
        }
        
        Element<ELEMENT_DIM,SPACE_DIM>* p_new_element=
            new Element<ELEMENT_DIM,SPACE_DIM>(*pElement, new_elt_index);
            
        // Second, update the node in the element with the new one
        p_new_element->UpdateNode(ELEMENT_DIM-1-i, mNodes[new_node_index]);
        
        try
        {
            p_new_element->RefreshJacobianDeterminant();
        }
        catch (Exception e)
        {
            assert(0); //This code should now be redundant and fail the coverage test
            //Clean up this element, because it's not going to be added to the mesh
            delete p_new_element;
            EXCEPTION("RefineElement could not be completed (point was not in element)");
        }
        
        // Third, add the new element to the set
        if ((unsigned) new_elt_index == mElements.size())
        {
            mElements.push_back(p_new_element);
        }
        else
        {
            delete mElements[new_elt_index];
            mElements[new_elt_index] = p_new_element;
        }
        
    }
    
    // Lastly, update the last node in the element to be refined
    pElement->UpdateNode(ELEMENT_DIM, mNodes[new_node_index]);
    try
    {
        pElement->RefreshJacobianDeterminant();
    }
    catch (Exception e)
    {
        //There's no memory management to do here, because this original element is
        //already referenced by the mesh
        EXCEPTION("RefineElement could not be completed (point was not in element)");
    }
    
    
    return new_node_index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateMeshVolume()
{
    double mesh_volume = 0.0;
    
    ElementIterator it = GetElementIteratorBegin();
    
    while (it != GetElementIteratorEnd())
    {
        mesh_volume += (*it)->GetJacobianDeterminant();
        it++;
    }
    
    if (ELEMENT_DIM == 2)
    {
        mesh_volume /= 2.0;
    }
    else if (ELEMENT_DIM == 3)
    {
        mesh_volume /= 6.0;
    }
    return mesh_volume;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateMeshSurface()
{
    //ELEMENT_DIM-1 is the dimension of the boundary element
    unsigned bound_element_dim=ELEMENT_DIM-1;
    assert(bound_element_dim < 3 &&  bound_element_dim>=0);
    if ( bound_element_dim == 0)
    {
        return 0.0;
    }
    
    double mesh_surface = 0.0;
    BoundaryElementIterator it = GetBoundaryElementIteratorBegin();
    
    while (it != GetBoundaryElementIteratorEnd())
    {
        mesh_surface += (*it)->GetJacobianDeterminant();
        it++;
    }
    
    if ( bound_element_dim == 2)
    {
        mesh_surface /= 2.0;
    }
    
    return mesh_surface;
}



/**
 * Scale the mesh.
 * @param xFactor is the scale in the x-direction,
 * @param yFactor is the scale in the y-direction,
 * @param zFactor is the scale in the z-direction
 **/
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Scale(
    const double xScale,
    const double yScale,
    const double zScale)
{
    unsigned num_nodes=GetNumAllNodes();
    
    for (unsigned i=0; i<num_nodes; i++)
    {
        c_vector<double, SPACE_DIM>& r_location = mNodes[i]->rGetModifiableLocation();
        if (SPACE_DIM>=3)
        {
            r_location[2] *= zScale;
        }
        if (SPACE_DIM>=2)
        {
            r_location[1] *= yScale;
        }
        r_location[0] *= xScale;
    }
    
    RefreshMesh();
}

/**
 * Translate the mesh.
 * @param xMovement is the x-displacement,
 * @param yMovement is the y-displacement,
 * @param zMovement is the z-displacement,
 **/
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Translate(
    const double xMovement,
    const double yMovement,
    const double zMovement)
{
    c_vector<double , SPACE_DIM> displacement;
    
    switch (SPACE_DIM)
    {
        case 3:
            displacement[2]=zMovement;
        case 2:
            displacement[1]=yMovement;
        case 1:
            displacement[0]=xMovement;
    }
    
    Translate(displacement);
}


/**
 * Translate mesh using the BOOST ublas library - this is the method that actually does the work
 * @param transVec is a translation vector of the correct size
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Translate(c_vector<double, SPACE_DIM> transVec)
{
    unsigned num_nodes=GetNumAllNodes();
    
    for (unsigned i=0; i<num_nodes; i++)
    {
        c_vector<double, SPACE_DIM>& r_location = mNodes[i]->rGetModifiableLocation();
        r_location += transVec;
    }
    
    RefreshMesh();
}




/**
 * Do a general mesh rotation with an +ve determinant orthonormal rotation_matrix
 * This is the method that actually does the work
 * @param rotation_matrix is a Ublas rotation matrix of the correct form
 **/
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Rotate(
    c_matrix<double , SPACE_DIM, SPACE_DIM> rotation_matrix)
{
    unsigned num_nodes=GetNumAllNodes();
    for (unsigned i=0; i<num_nodes; i++)
    {
        c_vector<double, SPACE_DIM>& r_location = mNodes[i]->rGetModifiableLocation();
        r_location = prod(rotation_matrix, r_location);
    }
    
    RefreshMesh();
}

/**
* Do an angle axis rotation
* @param axis is the axis of rotation (does not need to be normalised)
* @param angle is the angle in radians
 **/
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Rotate(c_vector<double,3> axis, double angle)
{
    assert(SPACE_DIM == 3);
    double norm = norm_2(axis);
    c_vector<double,3> unit_axis=axis/norm;
    
    c_matrix<double, SPACE_DIM,SPACE_DIM> rotation_matrix;
    
    double c = cos(angle);
    double s = sin(angle);
    
    rotation_matrix(0,0) = unit_axis(0)*unit_axis(0)+c*(1-unit_axis(0)*unit_axis(0));
    rotation_matrix(0,1) = unit_axis(0)*unit_axis(1)*(1-c) - unit_axis(2)*s;
    rotation_matrix(1,0) = unit_axis(0)*unit_axis(1)*(1-c) + unit_axis(2)*s;
    rotation_matrix(1,1) = unit_axis(1)*unit_axis(1)+c*(1-unit_axis(1)*unit_axis(1));
    rotation_matrix(0,2) = unit_axis(0)*unit_axis(2)*(1-c)+unit_axis(1)*s;
    rotation_matrix(1,2) = unit_axis(1)*unit_axis(2)*(1-c)-unit_axis(0)*s;
    rotation_matrix(2,0) = unit_axis(0)*unit_axis(2)*(1-c)-unit_axis(1)*s;
    rotation_matrix(2,1) = unit_axis(1)*unit_axis(2)*(1-c)+unit_axis(0)*s;
    rotation_matrix(2,2) = unit_axis(2)*unit_axis(2)+c*(1-unit_axis(2)*unit_axis(2));
    
    Rotate(rotation_matrix);
}


/**
 * Rotate the mesh about the x-axis
 * @param theta is the angle in radians
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RotateX(const double theta)
{
    if (SPACE_DIM != 3)
    {
        EXCEPTION("This rotation is only valid in 3D");
    }
    c_matrix<double , SPACE_DIM, SPACE_DIM> x_rotation_matrix=identity_matrix<double>(SPACE_DIM);
    
    x_rotation_matrix(1,1) = cos(theta);
    x_rotation_matrix(1,2) = sin(theta);
    x_rotation_matrix(2,1) = -sin(theta);
    x_rotation_matrix(2,2) = cos(theta);
    Rotate(x_rotation_matrix);
}


/**
 * Rotate the mesh about the y-axis
 * @param theta is the angle in radians
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RotateY(const double theta)
{
    if (SPACE_DIM != 3)
    {
        EXCEPTION("This rotation is only valid in 3D");
    }
    c_matrix<double , SPACE_DIM, SPACE_DIM> y_rotation_matrix=identity_matrix<double>(SPACE_DIM);
    
    y_rotation_matrix(0,0) = cos(theta);
    y_rotation_matrix(0,2) = -sin(theta);
    y_rotation_matrix(2,0) = sin(theta);
    y_rotation_matrix(2,2) = cos(theta);
    
    
    Rotate(y_rotation_matrix);
}

/**
 * Rotate the mesh about the z-axis
 * @param theta is the angle in radians
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RotateZ(const double theta)
{
    if (SPACE_DIM < 2)
    {
        EXCEPTION("This rotation is not valid in less than 2D");
    }
    c_matrix<double , SPACE_DIM, SPACE_DIM> z_rotation_matrix=identity_matrix<double>(SPACE_DIM);
    
    
    z_rotation_matrix(0,0) = cos(theta);
    z_rotation_matrix(0,1) = sin(theta);
    z_rotation_matrix(1,0) = -sin(theta);
    z_rotation_matrix(1,1) = cos(theta);
    
    Rotate(z_rotation_matrix);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::DeleteBoundaryNodeAt(unsigned index)
{
    if (!mNodes[index]->IsBoundaryNode() )
    {
        EXCEPTION(" You may only delete a boundary node ");
    }
    
    mNodes[index]->MarkAsDeleted();
    mDeletedNodeIndices.push_back(index);
    // Update the boundary node vector
    typename std::vector<Node<SPACE_DIM>*>::iterator b_node_iter
    = std::find(mBoundaryNodes.begin(), mBoundaryNodes.end(), mNodes[index]);
    mBoundaryNodes.erase(b_node_iter);
    
    // Remove boundary elements containing this node
    std::set<unsigned> boundary_element_indices = mNodes[index]->rGetContainingBoundaryElementIndices();
    std::set<unsigned>::const_iterator boundary_element_indices_iterator = boundary_element_indices.begin();
    while (boundary_element_indices_iterator != boundary_element_indices.end())
    {
        BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* p_boundary_element = GetBoundaryElement(*boundary_element_indices_iterator);
        p_boundary_element->MarkAsDeleted();
        mDeletedBoundaryElementIndices.push_back(*boundary_element_indices_iterator);
        boundary_element_indices_iterator++;
    }
    
    // Remove elements containing this node
    std::set<unsigned> element_indices = mNodes[index]->rGetContainingElementIndices();
    std::set<unsigned>::const_iterator element_indices_iterator = element_indices.begin();
    while (element_indices_iterator != element_indices.end())
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(*element_indices_iterator);
        for (unsigned i=0 ; i< p_element->GetNumNodes();i++)
        {
            Node<SPACE_DIM>* p_node = p_element->GetNode(i);
            if (!p_node->IsDeleted())
            {
                p_node->SetAsBoundaryNode();
                // Update the boundary node vector
                mBoundaryNodes.push_back(p_node);
            }
        }
        p_element->MarkAsDeleted();
        mDeletedElementIndices.push_back(p_element->GetIndex());
        element_indices_iterator++;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReIndex()
{

    unsigned to_go=0;
    while (to_go<mElements.size())
    {
        while (to_go<mElements.size() && mElements[to_go]->IsDeleted() == false)
        {
            to_go++;
        }
        if (to_go < mElements.size())
        {
            unsigned to_stay=to_go;
            while (to_stay<mElements.size() && mElements[to_stay]->IsDeleted())
            {
                to_stay++;
            }
            if (to_stay<mElements.size())
            {
                //Swap the elements
                Element<ELEMENT_DIM,SPACE_DIM> *p_element=mElements[to_go];
                mElements[to_go]=mElements[to_stay];
                mElements[to_stay]=p_element;
                to_go++;
            }
            else
            {
                break;//Nothing left to swap with
            }
        }
    }
    while (mDeletedElementIndices.size() != 0)
    {
        mDeletedElementIndices.pop_back();
        assert(mElements.back()->IsDeleted());
        mElements.pop_back();
    }
    assert(mElements.back()->IsDeleted() == false );
    
    
    to_go=0;
    while (to_go<mNodes.size())
    {
        while (to_go<mNodes.size() && mNodes[to_go]->IsDeleted() == false)
        {
            to_go++;
        }
        if (to_go < mNodes.size())
        {
            unsigned to_stay=to_go;
            while (to_stay<mNodes.size() && mNodes[to_stay]->IsDeleted())
            {
                to_stay++;
            }
            if (to_stay<mNodes.size())
            {
                //Swap the elements
                Node<SPACE_DIM> *p_node=mNodes[to_go];
                mNodes[to_go]=mNodes[to_stay];
                mNodes[to_stay]=p_node;
                to_go++;
            }
            else
            {
                break;//Nothing left to swap with
            }
        }
    }
    while (mDeletedNodeIndices.size() != 0)
    {
        mDeletedNodeIndices.pop_back();
        assert(mNodes.back()->IsDeleted());
        mNodes.pop_back();
    }
    assert(mNodes.back()->IsDeleted() == false );
    
    to_go=0;
    while (to_go<mBoundaryElements.size())
    {
        while (to_go<mBoundaryElements.size() && mBoundaryElements[to_go]->IsDeleted() == false)
        {
            to_go++;
        }
        if (to_go < mBoundaryElements.size())
        {
            unsigned to_stay=to_go;
            while (to_stay<mBoundaryElements.size() && mBoundaryElements[to_stay]->IsDeleted())
            {
                to_stay++;
            }
            if (to_stay<mBoundaryElements.size())
            {
                //Swap the elements
                BoundaryElement <ELEMENT_DIM-1,SPACE_DIM> *p_el=mBoundaryElements[to_go];
                mBoundaryElements[to_go]=mBoundaryElements[to_stay];
                mBoundaryElements[to_stay]=p_el;
                to_go++;
            }
            else
            {
                break;//Nothing left to swap with
            }
        }
    }
    while (mDeletedBoundaryElementIndices.size() != 0)
    {
        mDeletedBoundaryElementIndices.pop_back();
        assert(mBoundaryElements.back()->IsDeleted());
        mBoundaryElements.pop_back();
    }
    assert(mBoundaryElements.back()->IsDeleted() == false );
    
    
    for (unsigned i=0; i<mNodes.size();i++)
    {
    
        mNodes[i]->SetIndex(i);
    }
    for (unsigned i=0; i<mElements.size();i++)
    {
    
        mElements[i]->ResetIndex(i);
        
    }
    
    for (unsigned i=0; i<mBoundaryElements.size();i++)
    {
    
        mBoundaryElements[i]->ResetIndex(i);
        
    }
    
}
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh(NodeMap &map)
{
    //Make sure the map is big enough
    map.Reserve(GetNumAllNodes());
    
    //Make sure that we are in the correct dimension
    assert( SPACE_DIM==2 || SPACE_DIM==3 );
    assert( ELEMENT_DIM == SPACE_DIM );
    OutputFileHandler handler("");
    out_stream node_file=handler.OpenOutputFile("temp.node");
    
    (*node_file)<<GetNumNodes()<<"\t" << SPACE_DIM << "\t0\t0\n";
    
    unsigned new_index = 0;
    
    for (unsigned i=0; i<GetNumAllNodes(); i++)
    {
        if (mNodes[i]->IsDeleted())
        {
            map.SetDeleted(i);
        }
        else
        {
            map.SetNewIndex(i,new_index);
            new_index++;
            const c_vector<double, SPACE_DIM> node_loc = mNodes[i]->rGetLocation();
            (*node_file)<<i<<"\t"<<node_loc[0]<<"\t"<<node_loc[1];
            if (SPACE_DIM ==3)
            {
                (*node_file)<<"\t"<<node_loc[2];
            }
            (*node_file)<<"\n";
        }
    }
    
    node_file->close();
    
    //system("cat /tmp/chaste/testoutput/temp.node");
    
    
    std::string full_name = handler.GetTestOutputDirectory("")+"temp.";
    std::string binary_name;
    if (SPACE_DIM==2)
    {
        binary_name="triangle";
    }
    else
    {
        binary_name="tetgen";
    }
    std::string command   = "./bin/"+ binary_name +" -e "
                            + full_name + "node"
                            + " > /dev/null";
    system(command.c_str());
    
    //Read the new mesh back from file
    TrianglesMeshReader<SPACE_DIM,SPACE_DIM> mesh_reader(full_name+"1");
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> temporary_mesh;
    temporary_mesh.ConstructFromMeshReader(mesh_reader);
    
    //Delete current data
    
    // Iterate over nodes and free the memory
    
    
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        delete mNodes[i];
    }
    mNodes.clear();
    // Iterate over elements and free the memory
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    mElements.clear();
    // Iterate over boundary elements and free the memory
    for (unsigned i=0; i<mBoundaryElements.size(); i++)
    {
        delete mBoundaryElements[i];
    }
    mBoundaryElements.clear();
    
    mDeletedElementIndices.clear();
    mDeletedBoundaryElementIndices.clear();
    mDeletedNodeIndices.clear();
    mBoundaryNodes.clear();
    
    for (unsigned i=0; i<temporary_mesh.GetNumNodes(); i++)
    {
        c_vector<double, SPACE_DIM> node_loc = temporary_mesh.GetNode(i)->rGetLocation();
        bool is_boundary=temporary_mesh.GetNode(i)->IsBoundaryNode();
        Node<SPACE_DIM>* p_node=new Node<SPACE_DIM>(i,node_loc,is_boundary);
        mNodes.push_back(p_node);
        if (is_boundary)
        {
            mBoundaryNodes.push_back(p_node);
        }
    }
    
    for (unsigned i=0; i<temporary_mesh.GetNumElements(); i++)
    {
        std::vector<Node<SPACE_DIM>* > nodes;
        for (unsigned j=0; j<SPACE_DIM+1; j++)
        {
            unsigned index=(temporary_mesh.GetElement(i))->GetNodeGlobalIndex(j);
            nodes.push_back(mNodes[index]);
        }
        Element<ELEMENT_DIM,SPACE_DIM> *p_element=
            new Element<ELEMENT_DIM,SPACE_DIM>(i, nodes);
        mElements.push_back(p_element);
    }
    
    for (unsigned i=0; i<temporary_mesh.GetNumBoundaryElements(); i++)
    {
        std::vector<Node<SPACE_DIM>* > nodes;
        for (unsigned j=0; j<SPACE_DIM; j++)
        {
            unsigned index=(temporary_mesh.GetBoundaryElement(i))->GetNodeGlobalIndex(j);
            nodes.push_back(mNodes[index]);
        }
        BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> *p_b_element=
            new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(i, nodes);
        mBoundaryElements.push_back(p_b_element);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes()
{
    RandomNumberGenerator *p_rng=RandomNumberGenerator::Instance();
    
    //Working from the back, each node is swapped with a random node that precedes it in the array
    for (unsigned index=mNodes.size()-1; index>0; index--)
    {
        unsigned  other=p_rng->randMod(index+1); //includes the possibility of rolling "index"
        //Swap index and other
        Node<SPACE_DIM> *temp=mNodes[index];
        mNodes[index]=mNodes[other];
        mNodes[other]=temp;
    }
    
    //Update indices
    for (unsigned index=0; index<mNodes.size(); index++)
    {
        mNodes[index]->SetIndex(index);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes(std::vector<unsigned> perm)
{
    //Let's not do this if there are any deleted nodes
    assert( GetNumAllNodes() == GetNumNodes());
    
    //Copy the node pointers
    std::vector <Node <SPACE_DIM> *> copy_m_nodes;
    copy_m_nodes.assign(mNodes.begin(), mNodes.end());
    
    
    for (unsigned i=0;i<mNodes.size();i++)
    {
        mNodes[ perm[i] ] = copy_m_nodes[i];
    }
    
    //Update indices
    for (unsigned index=0; index<mNodes.size(); index++)
    {
        mNodes[index]->SetIndex(index);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodesWithMetisBinaries()
{
    assert( ELEMENT_DIM==2 || ELEMENT_DIM==3 );
    assert( GetNumAllElements() == GetNumElements());
    assert( GetNumAllNodes() == GetNumNodes());
    
    //Open a file for the elements
    OutputFileHandler handler("");
    out_stream metis_file=handler.OpenOutputFile("metis.mesh");
    
    (*metis_file)<<GetNumElements()<<"\t";
    if (ELEMENT_DIM==2)
    {
        (*metis_file)<<1<<"\n"; //1 is Metis speak for triangles
    }
    else
    {
        (*metis_file)<<2<<"\n"; //2 is Metis speak for tetrahedra
    }
    
    for (unsigned i=0; i<(unsigned)GetNumElements(); i++)
    {
        for (unsigned j=0; j<ELEMENT_DIM+1; j++)
        {
            //Note the +1 since Metis wants meshes indexed from 1
            (*metis_file)<<mElements[i]->GetNode(j)->GetIndex() + 1<<"\t";
        }
        (*metis_file)<<"\n";
    }
    metis_file->close();
    
    
    std::string convert_command   = "./bin/mesh2nodal "+handler.GetTestOutputDirectory("")
                                    + "metis.mesh"
                                    + " > /dev/null";
    system(convert_command.c_str());
    
    std::string permute_command   = "./bin/onmetis "+handler.GetTestOutputDirectory("")
                                    + "metis.mesh.ngraph"
                                    + " > /dev/null";
    system(permute_command.c_str());
    
    //Read the permutation back into a std::vector
    std::string perm_file_name   = handler.GetTestOutputDirectory("")
                                   + "metis.mesh.ngraph.iperm";
    std::ifstream perm_file(perm_file_name.c_str());
    std::vector<unsigned> perm;
    for (unsigned i=0; i<(unsigned)GetNumNodes(); i++)
    {
        unsigned new_index;
        perm_file>>new_index;
        perm.push_back(new_index);
    }
    perm_file.close();
    
    PermuteNodes(perm);
    
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CheckVoronoi(Element<ELEMENT_DIM, SPACE_DIM> *pElement, double maxPenetration)
{
    assert (ELEMENT_DIM == SPACE_DIM);
    unsigned num_nodes = pElement->GetNumNodes();
    std::set<unsigned> neighbouring_elements_indices;
    std::set< Element<ELEMENT_DIM,SPACE_DIM> *> neighbouring_elements;
    std::set<unsigned> neighbouring_nodes_indices;
    
    //Form a set of neighbouring elements via the nodes
    for (unsigned i = 0 ; i < num_nodes; i++)
    {
        Node<SPACE_DIM>* node = pElement->GetNode(i);
        neighbouring_elements_indices = node->rGetContainingElementIndices();
        ///\todo Should use a set union operation here
        for (std::set<unsigned>::const_iterator it = neighbouring_elements_indices.begin();
                 it != neighbouring_elements_indices.end(); ++it)
            {
                neighbouring_elements.insert(GetElement(*it));
            }
    }
    neighbouring_elements.erase(pElement);
    
    //For each neighbouring element find the supporting nodes
    typedef typename std::set<Element<ELEMENT_DIM,SPACE_DIM> *>::const_iterator ElementIterator;
    
    for (ElementIterator it = neighbouring_elements.begin();
         it != neighbouring_elements.end(); ++it)
    {
        for (unsigned i = 0 ; i < num_nodes; i++)
        {
            neighbouring_nodes_indices.insert((*it)->GetNodeGlobalIndex(i));
        }
    }
    //Remove the nodes that support this element
    for (unsigned i = 0 ; i < num_nodes; i++)
    {
        neighbouring_nodes_indices.erase(pElement->GetNodeGlobalIndex(i));
    }
    
    //Get the circumsphere information
    c_vector <double, ELEMENT_DIM+1> this_circum_centre;
    this_circum_centre = pElement->CalculateCircumsphere();
    
    //Copy the actualy circumcentre into a smaller vector
    c_vector <double, ELEMENT_DIM> circum_centre;
    for (unsigned i=0;i<ELEMENT_DIM;i++)
    {
        circum_centre[i]=this_circum_centre[i];
    }
    
    for (std::set<unsigned>::const_iterator it = neighbouring_nodes_indices.begin();
             it != neighbouring_nodes_indices.end(); ++it)
        {
            c_vector <double, ELEMENT_DIM> node_location = GetNode(*it)->rGetLocation();
            
            // Calculate vector from circumcenter to node
            node_location -= circum_centre;
            // This is to calculate the squared distance betweeen them
            double squared_distance = inner_prod(node_location, node_location);
            
            // If the squared idstance is less than the elements circum-radius(squared),
            // then the voronoi property is violated.
            
            if (squared_distance < this_circum_centre[ELEMENT_DIM])
            {
                // We know the node is inside the circumsphere, but we don't know how far
                double radius = sqrt(this_circum_centre[ELEMENT_DIM]);
                double distance = radius - sqrt(squared_distance);
                
                // If the node penetration is greater than supplied maximum penetration factor
                if (distance/radius > maxPenetration)
                {
                    return false;
                }
            }
        }
    return true;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CheckVoronoi(double maxPenetration)
{
    // Looping through all the elements in the mesh
    for (unsigned i=0; i < mElements.size();i++)
    {
        // Check if the element is not deleted
        if (!mElements[i]->IsDeleted())
        {
            // Checking the Voronoi of the Element
            if (CheckVoronoi(mElements[i], maxPenetration) == false)
            {
                return false;
            }
        }
    }
    return true;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructRectangularMesh(unsigned width, unsigned height, bool stagger)
{
    assert(SPACE_DIM == 2);
    assert(ELEMENT_DIM == 2);
    
    //Construct the nodes
    unsigned node_index=0;
    for (int j=(int)height;j>=0;j--) //j must be signed for this loop to terminate
    {
        for (unsigned i=0;i<width+1;i++)
        {
            bool is_boundary=false;
            if (i==0 || j==0 || i==width || j==(int)height)
            {
                is_boundary=true;
            }
            Node<2>* p_node = new Node<2>(node_index++, is_boundary, i, j);
            mNodes.push_back(p_node);
            if(is_boundary)
            {
                mBoundaryNodes.push_back(p_node);
            }
        }
    }
    
    //Construct the boundary elements
    unsigned belem_index=0;
    //Top
    for (unsigned i=0;i<width;i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(mNodes[i]);
        nodes.push_back(mNodes[i+1]);
        mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Right
    for (unsigned i=1;i<height+1;i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(mNodes[(width+1)*i-1]);
        nodes.push_back(mNodes[(width+1)*(i+1)-1]);
        mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Bottom
    for (unsigned i=0;i<width;i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(mNodes[height*(width+1)+i+1]);
        nodes.push_back(mNodes[height*(width+1)+i]);
        mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Left
    for (unsigned i=0;i<height;i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(mNodes[(width+1)*(i+1)]);
        nodes.push_back(mNodes[(width+1)*(i)]);
        mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    
    //Construct the elements
    unsigned elem_index=0;
    for (unsigned j=0;j<height;j++)
    {
        for (unsigned i=0; i<width; i++)
        {
            unsigned parity=(i+j)%2;
            std::vector<Node<SPACE_DIM>*> upper_nodes;
            upper_nodes.push_back(mNodes[j*(width+1)+i]);
            upper_nodes.push_back(mNodes[j*(width+1)+i+1]);
            if (stagger==false  || parity == 0)
            {
                upper_nodes.push_back(mNodes[(j+1)*(width+1)+i+1]);
            }
            else
            {
                upper_nodes.push_back(mNodes[(j+1)*(width+1)+i]);
            }
            mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++,upper_nodes));
            std::vector<Node<SPACE_DIM>*> lower_nodes;
            lower_nodes.push_back(mNodes[(j+1)*(width+1)+i+1]);
            lower_nodes.push_back(mNodes[(j+1)*(width+1)+i]);
            if (stagger==false  ||parity == 0)
            {
                lower_nodes.push_back(mNodes[j*(width+1)+i]);
            }
            else
            {
                lower_nodes.push_back(mNodes[j*(width+1)+i+1]);
            }
            mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++,lower_nodes));
        }
    }
    
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndex(Point<SPACE_DIM> testPoint, bool strict)
{
    ///\todo This ought to return a set of all elements that contain the point (if the point is a node in the mesh then it's contained in multiple elements)
    ///\todo Polling every element is unnecessary.  We ought to start from a likely place and hill climb
    for (unsigned i=0; i < mElements.size();i++)
    {
        ///\todo What if the element is deleted?
        if (mElements[i]->IncludesPoint(testPoint, strict))
        {
            return i;
        }
    }
    
    //If it's in none of the elements, then throw
    EXCEPTION("Point is not in mesh");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNearestElementIndex(Point<SPACE_DIM> testPoint)
{
    ///\todo This ought to return a set of all elements that contain the point (if the point is a node in the mesh then it's contained in multiple elements)
    ///\todo Polling every element is unnecessary.  We ought to start from a likely place and hill climb
    
    double max_min_weight=-INFINITY;
    unsigned closest_index=0;
    for (unsigned i=0; i < mElements.size();i++)
    {
        ///\todo What if the element is deleted?
        c_vector<double, ELEMENT_DIM+1> weight=mElements[i]->CalculateInterpolationWeights(testPoint);
        double min_weight=1.0;
        for (unsigned j=0; j<=ELEMENT_DIM; j++)
        {
            if (weight[j]<min_weight)
            {
                min_weight=weight[j];
            }
        }
        if (min_weight > max_min_weight)
        {
            max_min_weight = min_weight;
            closest_index=i;
        }
        
    }
    return closest_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndices(Point<SPACE_DIM> testPoint)
{
    std::vector<unsigned> element_indices;
    for (unsigned i=0; i < mElements.size();i++)
    {
        ///\todo What if the element is deleted?
        if (mElements[i]->IncludesPoint(testPoint))
        {
            element_indices.push_back(i);
        }
    }
    return element_indices;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetElementOwnerships(unsigned lo, unsigned hi)
{
    assert(hi>=lo);
    for (unsigned element_index=0; element_index<mElements.size(); element_index++)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element=mElements[element_index];
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

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructCuboid(unsigned width,
        unsigned height,
        unsigned depth,
        bool stagger)
{
    assert(SPACE_DIM == 3);
    assert(ELEMENT_DIM == 3);
    //Construct the nodes
    
    unsigned node_index=0;
    for (unsigned k=0;k<depth+1;k++)
    {
        for (unsigned j=0;j<height+1;j++)
        {
            for (unsigned i=0;i<width+1;i++)
            {
                bool is_boundary=false;
                if (i==0 || j==0 || k==0 || i==width || j==height || k==depth)
                {
                    is_boundary=true;
                }
                
                Node<3>* p_node = new Node<3>(node_index++, is_boundary, i, j, k);

                mNodes.push_back(p_node);
                if(is_boundary)
                {
                    mBoundaryNodes.push_back(p_node);
                }
            }
        }
    }
    
    // Construct the elements
    
    unsigned elem_index=0;
    unsigned belem_index=0;
    unsigned element_nodes[4][6][4] = {{{0, 1, 5, 7}, {0, 1, 3, 7},
                                    {0, 2, 3, 7}, {0, 2, 6, 7},
                                    {0, 4, 6, 7}, {0, 4, 5, 7}},
                                    {{1, 0, 2, 6}, {1, 0, 4, 6},
                                    {1, 5, 4, 6}, {1, 5, 7, 6},
                                    {1, 3, 2, 6}, {1, 3, 7, 6}},
                                    {{2, 0, 1, 5}, {2, 0, 4, 5},
                                    {2, 3, 1, 5}, {2, 3, 7, 5},
                                    {2, 6, 4, 5}, {2, 6, 7, 5}},
                                    {{3, 1, 0, 4}, {3, 1, 5, 4},
                                    {3, 2, 0, 4}, {3, 2, 6, 4},
                                    {3, 7, 5, 4}, {3, 7, 6, 4}}};     
                                    
    std::vector<Node<SPACE_DIM>*> tetrahedra_nodes;
    
    for (unsigned k=0;k<depth;k++)
    {
        for (unsigned j=0;j<height;j++)
        {
            for (unsigned i=0;i<width;i++)
            {
                // Compute the nodes' index
                unsigned global_node_indices[8];
                unsigned local_node_index = 0;
                
                for (unsigned z = 0; z < 2; z++)
                {
                    for (unsigned y = 0; y < 2; y++)
                    {
                        for (unsigned x = 0; x < 2; x++)
                        {
                            global_node_indices[local_node_index] = i+x+(width+1)*(j+y+(height+1)*(k+z));
                            
                            local_node_index++;
                        }
                    }
                }
                
                for (unsigned m = 0; m < 6; m++)
                {
                    // Tetrahedra #m
                    
                    tetrahedra_nodes.clear();
                                        
                    for (unsigned n = 0; n < 4; n++)
                    {
                        if (stagger)
                        {
                            if (i%2==0)
                            {
                                if (j%2==0)
                                {
                                    if (k%2==0)
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[0][m][n]]]);    
                                    }
                                    else
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[3][m][n]]]);
                                    }
                                }
                                else
                                {
                                    if (k%2==0)
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[2][m][n]]]);
                                    }
                                    else
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[1][m][n]]]);
                                    }
                                }
                            }
                            else
                            {
                               if (j%2==0)
                               {
                                    if (k%2==0)
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[3][m][n]]]);
                                    }
                                    else
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[1][m][n]]]);
                                    }
                               }
                               else
                               {
                                    if (k%2==0)
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[2][m][n]]]);
                                    }
                                    else
                                    {
                                        tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[0][m][n]]]);
                                    } 
                               }
                            }
                        
                        }
                                    
                        else
                        {
                            tetrahedra_nodes.push_back(mNodes[global_node_indices[element_nodes[0][m][n]]]);
                        }
                    }
                    
                    mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++, tetrahedra_nodes));
                }
                
                //Are we at a boundary?
                std::vector<Node<SPACE_DIM>*> triangle_nodes;
                if (i == 0) //low face at x==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[2]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[6]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[6]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[4]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (i == width-1) //high face at x=width
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[1]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[5]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[7]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[1]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[7]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[3]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (j == 0) //low face at y==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[5]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[1]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[4]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[5]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (j == height-1) //high face at y=height
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[2]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[3]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[7]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[2]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[7]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[6]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (k == 0) //low face at z==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[3]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[2]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[1]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[3]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (k == depth-1) //high face at z=depth
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[4]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[7]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[5]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(mNodes[global_node_indices[4]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[6]]);
                    triangle_nodes.push_back(mNodes[global_node_indices[7]]);
                    mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                
            }//i
        }//j
    }//k
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    // three loops, just like the destructor. note we don't delete boundary nodes.
    for (unsigned i=0; i<mBoundaryElements.size(); i++)
    {
        delete mBoundaryElements[i];
    }
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
    mBoundaryElements.clear();
    mDeletedElementIndices.clear();
    mDeletedBoundaryElementIndices.clear();
    mDeletedNodeIndices.clear();
    mBoundaryNodes.clear();
    
    mNumCornerNodes = 0;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateBoundaryOfFlaggedRegion()
{
    // a set of nodes which lie on the face, size 3 in 2D, size 4 in 3D
    typedef std::set<unsigned> FaceNodes; 

    // face maps to true the first time it is encountered, and false subsequent 
    // times. Thus, faces mapping to true at the end are boundary faces
    std::map<FaceNodes,bool> face_on_boundary;

    // loop over all elements
    ElementIterator iter = GetElementIteratorBegin();
    while (iter != GetElementIteratorEnd())
    {
        if((*iter)->IsFlagged())
        {
            // to get faces, initially start with all nodes..
            std::set<unsigned> all_nodes;
            for(unsigned i=0; i<ELEMENT_DIM+1; i++)
            {
                all_nodes.insert( (*iter)->GetNodeGlobalIndex(i) );
            }
            
            // remove one node in turn to obtain each face
            for(unsigned i=0; i<ELEMENT_DIM+1; i++)
            {
                FaceNodes face_nodes = all_nodes;
                face_nodes.erase( (*iter)->GetNodeGlobalIndex(i) );

                // search the map of faces to see if it contains this face                 
                std::map<FaceNodes,bool>::iterator it=face_on_boundary.find(face_nodes);

                if(it == face_on_boundary.end())
                {
                    // face not found, add and assume on boundary
                    face_on_boundary[face_nodes]=true;
                }
                else
                {
                    // face found in map, so not on boundary
                    it->second = false;
                }
            }
            
        }
        iter++;
    }
    
    // boundary nodes to be returned
    std::set<unsigned> boundary_of_flagged_region;
    
    // get all faces in the map
    std::map<FaceNodes,bool>::iterator it=face_on_boundary.begin();
    while(it!=face_on_boundary.end())
    {
        // if the face maps to true it is on the boundary
        if(it->second==true)
        {
            // get all nodes in the face and put in set to be returned
            boundary_of_flagged_region.insert(it->first.begin(),it->first.end());
        }
        it++;
    }
        
    return boundary_of_flagged_region;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocationA, const c_vector<double, SPACE_DIM>& rLocationB)
{
    c_vector<double, SPACE_DIM> vector = rLocationB - rLocationA;
    
    return vector;
}



#endif // _CONFORMINGTETRAHEDRALMESH_CPP_

