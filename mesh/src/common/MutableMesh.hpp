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

#ifndef MUTABLEMESH_HPP_
#define MUTABLEMESH_HPP_

//#include <boost/serialization/access.hpp>
//
////Jonathan Shewchuk's triangle
//#define REAL double
//#define VOID void
//#include "triangle.h"
//#undef REAL
//
//#include <boost/serialization/export.hpp>

#include "TetrahedralMesh.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MutableMesh : public TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>
{
protected:

    /// Indices of elements/nodes that have been deleted - these indices can be reused when adding
    /// new elements/nodes
    std::vector<unsigned> mDeletedElementIndices;
    std::vector<unsigned> mDeletedBoundaryElementIndices;
    std::vector<unsigned> mDeletedNodeIndices;
    bool mAddedNodes;


private:

    /**
     * Check whether any neighbouring node is inside the circumsphere of this element.
     * @param pointer to an element
     * @param maxPenetration is the maximum distance a node is allowed to be inside the
     * circumsphere of the element, as a proportion of the circumsphere radius.
     */
    bool CheckVoronoi(Element<ELEMENT_DIM, SPACE_DIM>  *pElement, double maxPenetration);

public:
    MutableMesh();
    MutableMesh(unsigned numElements);
    MutableMesh(std::vector<Node<SPACE_DIM> *> nodes);

    virtual ~MutableMesh();

    void Clear();

    unsigned GetNumNodes() const;
    unsigned GetNumElements() const;
    unsigned GetNumBoundaryElements() const;
    ///should unsigned GetNumBoundaryNodes() be overloaded too??
    
    void RescaleMeshFromBoundaryNode(ChastePoint<1> updatedPoint, unsigned boundaryNodeIndex);
    
    /**
     * Add a node to the mesh.
     *
     * NB. After calling this one or more times, you must then call ReMesh
     *
     */
    virtual unsigned AddNode(Node<SPACE_DIM> *pNewNode);

    virtual void SetNode(unsigned index, ChastePoint<SPACE_DIM> point, bool concreteMove=true);
    void MoveMergeNode(unsigned index, unsigned targetIndex, bool concreteMove=true);

    void DeleteNode(unsigned index);

    void DeleteNodePriorToReMesh(unsigned index);

    unsigned RefineElement(Element<ELEMENT_DIM,SPACE_DIM>* pElement, ChastePoint<SPACE_DIM> Point);

    /**
     * Remove a boundary node, and update all the appropriate data structures.
     *
     * The deleted node is not removed from the list, merely marked as deleted,
     * and can be reused when a new node is added to the mesh.
     *
     * Any elements or boundary elements containing this node will be removed.
     * The boundary nodes information will be updated with new boundary node(s).
     * NB: New boundary elements WILL NOT be added.
     *
     * @param index  The index of the node to remove.
     */
    void DeleteBoundaryNodeAt(unsigned index);

    /**
     * Re-index a mesh so that it has no deleted elements or nodes
     */
    void ReIndex(NodeMap& map);

    /**
     * Re-mesh a mesh using triangle (via library calls) or tetgen
     * @param map is a NodeMap which associates the indices of nodes in the old mesh
     * with indices of nodes in the new mesh.  This should be created with the correct size (NumAllNodes)
     */
    virtual void ReMesh(NodeMap& map);

    /**
     * Alternative version of remesh which takes no parameters does not require a NodeMap. Note: inherited
     * classes should overload ReMesh(NodeMap&)
     */
    void ReMesh();
    
    /**
     * Checks the entire mesh element by element and checks whether any neighbouring node
     * is inside the circumsphere of this element.
     * @param maxPenetration is the maximum distance a node is allowed to be inside the
     * circumsphere of an element that it is not a member of, as a proportion of the
     * circumsphere radius.
     */
    bool CheckVoronoi(double maxPenetration=0.0);    

};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableMesh<ELEMENT_DIM, SPACE_DIM>::MutableMesh()
{
    mAddedNodes = false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableMesh<ELEMENT_DIM, SPACE_DIM>::MutableMesh(unsigned numElements)
{
    mAddedNodes = false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableMesh<ELEMENT_DIM, SPACE_DIM>::MutableMesh(std::vector<Node<SPACE_DIM> *> nodes)
{
    Clear();
    for (unsigned index=0; index<nodes.size(); index++)
    {
        Node<SPACE_DIM>* temp_node = nodes[index];
        this->mNodes.push_back(temp_node);
    }
    mAddedNodes = true;
    NodeMap node_map(nodes.size());
    ReMesh(node_map);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableMesh<ELEMENT_DIM, SPACE_DIM>::~MutableMesh()
{
    Clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM> *pNewNode)
{

    if (mDeletedNodeIndices.empty())
    {
        pNewNode->SetIndex(this->mNodes.size());
        this->mNodes.push_back(pNewNode);
    }
    else
    {
        unsigned index=mDeletedNodeIndices.back();
        pNewNode->SetIndex(index);
        mDeletedNodeIndices.pop_back();
        delete this->mNodes[index];
        this->mNodes[index] = pNewNode;
    }
    mAddedNodes = true;
    return pNewNode->GetIndex();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    mDeletedElementIndices.clear();
    mDeletedBoundaryElementIndices.clear();
    mDeletedNodeIndices.clear();
    mAddedNodes = false;

    TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Clear();        
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements() const
{
    return this->mBoundaryElements.size() - mDeletedBoundaryElementIndices.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return this->mElements.size() - mDeletedElementIndices.size();
}


/// Returns the number of nodes that are actually in use
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size() - mDeletedNodeIndices.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::RescaleMeshFromBoundaryNode(ChastePoint<1> updatedPoint, unsigned boundaryNodeIndex)
{
    assert(this->GetNode(boundaryNodeIndex)->IsBoundaryNode());
    double scaleFactor = updatedPoint[0] / this->GetNode(boundaryNodeIndex)->GetPoint()[0];
    double temp;
    for (unsigned i=0; i < boundaryNodeIndex+1; i++)
    {
        temp = scaleFactor * this->mNodes[i]->GetPoint()[0];
        ChastePoint<1> newPoint(temp);
        this->mNodes[i]->SetPoint(newPoint);
    }
    this->RefreshMesh();
}

/**
 *  SetNode moves the node with a particular index to a new point in space and
  * verifies that the signed areas of the supporting Elements are positive
  * @param index is the index of the node to be moved
  * @param point is the new target location of the node
  * @param concreteMove is set to false if we want to skip the signed area tests
  */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::SetNode(unsigned index,
        ChastePoint<SPACE_DIM> point,
        bool concreteMove)
{
    this->mNodes[index]->SetPoint(point);
    if (concreteMove)
    {
        for (typename Node<SPACE_DIM>::ContainingElementIterator it = this->mNodes[index]->ContainingElementsBegin();
             it != this->mNodes[index]->ContainingElementsEnd();
             ++it)
        {
            try
            {
                GetElement(*it)->RefreshJacobianDeterminant();
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
        for (typename Node<SPACE_DIM>::ContainingBoundaryElementIterator it = this->mNodes[index]->ContainingBoundaryElementsBegin();
             it != this->mNodes[index]->ContainingBoundaryElementsEnd();
             ++it)
        {
            try
            {
                GetBoundaryElement(*it)->RefreshJacobianDeterminant();
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
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::DeleteNode(unsigned index)
{
    if (this->mNodes[index]->IsDeleted())
    {
        EXCEPTION("Trying to delete a deleted node");
    }
    unsigned target_index;
    bool found_target=false;
    for (typename Node<SPACE_DIM>::ContainingElementIterator it = this->mNodes[index]->ContainingElementsBegin();
         !found_target && it != this->mNodes[index]->ContainingElementsEnd();
         ++it)
    {
        Element <ELEMENT_DIM,SPACE_DIM> *p_element = GetElement(*it);
        for (unsigned i=0; i<=ELEMENT_DIM && !found_target; i++)
        {
            target_index = p_element->GetNodeGlobalIndex(i);
            try
            {
                MoveMergeNode(index, target_index, false);
                found_target = true;
            }
            catch (Exception e)
            {
                // Just try the next node
            }
        }
    }
    if (!found_target)
    {
        EXCEPTION("Failure to delete node");
    }

    MoveMergeNode(index, target_index);
}

/**
 * This marks a node as deleted. Note that it DOES NOT deal with the
 * associated elements and therefore should only be called immediately prior
 * to a ReMesh() being called. (Thus saves work compared to DeleteNode()
 * function and does not MoveMerge the node and elements).
 *
 * @param index The index of the node to delete
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::DeleteNodePriorToReMesh(unsigned index)
{
#define COVERAGE_IGNORE
    // A ReMesh can only happen in 2D or 3D so
    assert(SPACE_DIM==2 || SPACE_DIM==3);
#undef COVERAGE_IGNORE
    this->mNodes[index]->MarkAsDeleted();
    mDeletedNodeIndices.push_back(index);
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
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::MoveMergeNode(unsigned index,
        unsigned targetIndex,
        bool concreteMove)
{

    if (this->mNodes[index]->IsDeleted())
    {
        EXCEPTION("Trying to move a deleted node");
    }

    if (index == targetIndex)
    {
        EXCEPTION("Trying to merge a node with itself");
    }
    if (this->mNodes[index]->IsBoundaryNode())
    {
        if (!this->mNodes[targetIndex]->IsBoundaryNode())
        {
            EXCEPTION("A boundary node can only be moved on to another boundary node");
        }
    }
    std::set<unsigned> unshared_element_indices;
    std::set_difference(this->mNodes[index]->rGetContainingElementIndices().begin(),
                        this->mNodes[index]->rGetContainingElementIndices().end(),
                        this->mNodes[targetIndex]->rGetContainingElementIndices().begin(),
                        this->mNodes[targetIndex]->rGetContainingElementIndices().end(),
                        std::inserter(unshared_element_indices, unshared_element_indices.begin()));


    if (unshared_element_indices.size() == this->mNodes[index]->rGetContainingElementIndices().size())
    {
        EXCEPTION("These nodes cannot be merged since they are not neighbours");
    }

    std::set<unsigned> unshared_boundary_element_indices;
    std::set_difference(this->mNodes[index]->rGetContainingBoundaryElementIndices().begin(),
                        this->mNodes[index]->rGetContainingBoundaryElementIndices().end(),
                        this->mNodes[targetIndex]->rGetContainingBoundaryElementIndices().begin(),
                        this->mNodes[targetIndex]->rGetContainingBoundaryElementIndices().end(),
                        std::inserter(unshared_boundary_element_indices, unshared_boundary_element_indices.begin()));


    if (this->mNodes[index]->IsBoundaryNode())
    {
        if (unshared_boundary_element_indices.size()
            == this->mNodes[index]->rGetContainingBoundaryElementIndices().size())
        {
            //May be redundant (only thrown in 1D tests)
            EXCEPTION("These nodes cannot be merged since they are not neighbours on the boundary");
        }
    }

    this->mNodes[index]->rGetModifiableLocation() = this->mNodes[targetIndex]->rGetLocation();

    for (std::set<unsigned>::const_iterator element_iter=unshared_element_indices.begin();
             element_iter != unshared_element_indices.end();
             element_iter++)
    {
        try
        {

            this->GetElement(*element_iter)->RefreshJacobianDeterminant(concreteMove);
            if (concreteMove)
            {
                this->GetElement(*element_iter)->ReplaceNode(this->mNodes[index], this->mNodes[targetIndex]);
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

        this->GetBoundaryElement(*boundary_element_iter)->RefreshJacobianDeterminant(concreteMove);
        if (concreteMove)
        {
            this->GetBoundaryElement(*boundary_element_iter)->ReplaceNode(this->mNodes[index], this->mNodes[targetIndex]);
        }
    }

    std::set<unsigned> shared_element_indices;
    std::set_intersection(this->mNodes[index]->rGetContainingElementIndices().begin(),
                          this->mNodes[index]->rGetContainingElementIndices().end(),
                          this->mNodes[targetIndex]->rGetContainingElementIndices().begin(),
                          this->mNodes[targetIndex]->rGetContainingElementIndices().end(),
                          std::inserter(shared_element_indices, shared_element_indices.begin()));

    for (std::set<unsigned>::const_iterator element_iter=shared_element_indices.begin();
             element_iter != shared_element_indices.end();
             element_iter++)
    {
        if (concreteMove)
        {
            this->GetElement(*element_iter)->MarkAsDeleted();
            mDeletedElementIndices.push_back(*element_iter);
        }
        else
        {
            this->GetElement(*element_iter)->ZeroJacobianDeterminant();
        }
    }


    std::set<unsigned> shared_boundary_element_indices;
    std::set_intersection(this->mNodes[index]->rGetContainingBoundaryElementIndices().begin(),
                          this->mNodes[index]->rGetContainingBoundaryElementIndices().end(),
                          this->mNodes[targetIndex]->rGetContainingBoundaryElementIndices().begin(),
                          this->mNodes[targetIndex]->rGetContainingBoundaryElementIndices().end(),
                          std::inserter(shared_boundary_element_indices, shared_boundary_element_indices.begin()));

    for (std::set<unsigned>::const_iterator boundary_element_iter=shared_boundary_element_indices.begin();
             boundary_element_iter != shared_boundary_element_indices.end();
             boundary_element_iter++)
    {
        if (concreteMove)
        {
            this->GetBoundaryElement(*boundary_element_iter)->MarkAsDeleted();
            mDeletedBoundaryElementIndices.push_back(*boundary_element_iter);
        }
        else
        {
            this->GetBoundaryElement(*boundary_element_iter)->ZeroJacobianDeterminant();
            this->GetBoundaryElement(*boundary_element_iter)->ZeroWeightedDirection();
        }
    }

    if (concreteMove)
    {
        this->mNodes[index]->MarkAsDeleted();
        mDeletedNodeIndices.push_back(index);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMesh<ELEMENT_DIM, SPACE_DIM>::RefineElement(
    Element<ELEMENT_DIM,SPACE_DIM>* pElement,
    ChastePoint<SPACE_DIM> point)
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
            new_elt_index = this->mElements.size();
        }
        else
        {
            new_elt_index = mDeletedElementIndices.back();
            mDeletedElementIndices.pop_back();
        }

        Element<ELEMENT_DIM,SPACE_DIM>* p_new_element=
            new Element<ELEMENT_DIM,SPACE_DIM>(*pElement, new_elt_index);

        // Second, update the node in the element with the new one
        p_new_element->UpdateNode(ELEMENT_DIM-1-i, this->mNodes[new_node_index]);

        p_new_element->RefreshJacobianDeterminant();

        // Third, add the new element to the set
        if ((unsigned) new_elt_index == this->mElements.size())
        {
            this->mElements.push_back(p_new_element);
        }
        else
        {
            delete this->mElements[new_elt_index];
            this->mElements[new_elt_index] = p_new_element;
        }

    }

    // Lastly, update the last node in the element to be refined
    pElement->UpdateNode(ELEMENT_DIM, this->mNodes[new_node_index]);
    pElement->RefreshJacobianDeterminant();

    return new_node_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::DeleteBoundaryNodeAt(unsigned index)
{
    if (!this->mNodes[index]->IsBoundaryNode() )
    {
        EXCEPTION(" You may only delete a boundary node ");
    }

    this->mNodes[index]->MarkAsDeleted();
    mDeletedNodeIndices.push_back(index);
    // Update the boundary node vector
    typename std::vector<Node<SPACE_DIM>*>::iterator b_node_iter
    = std::find(this->mBoundaryNodes.begin(), this->mBoundaryNodes.end(), this->mNodes[index]);
    this->mBoundaryNodes.erase(b_node_iter);

    // Remove boundary elements containing this node
    std::set<unsigned> boundary_element_indices = this->mNodes[index]->rGetContainingBoundaryElementIndices();
    std::set<unsigned>::const_iterator boundary_element_indices_iterator = boundary_element_indices.begin();
    while (boundary_element_indices_iterator != boundary_element_indices.end())
    {
        BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* p_boundary_element = this->GetBoundaryElement(*boundary_element_indices_iterator);
        p_boundary_element->MarkAsDeleted();
        mDeletedBoundaryElementIndices.push_back(*boundary_element_indices_iterator);
        boundary_element_indices_iterator++;
    }

    // Remove elements containing this node
    std::set<unsigned> element_indices = this->mNodes[index]->rGetContainingElementIndices();
    std::set<unsigned>::const_iterator element_indices_iterator = element_indices.begin();
    while (element_indices_iterator != element_indices.end())
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(*element_indices_iterator);
        for (unsigned i=0; i< p_element->GetNumNodes();i++)
        {
            Node<SPACE_DIM>* p_node = p_element->GetNode(i);
            if (!p_node->IsDeleted())
            {
                p_node->SetAsBoundaryNode();
                // Update the boundary node vector
                this->mBoundaryNodes.push_back(p_node);
            }
        }
        p_element->MarkAsDeleted();
        mDeletedElementIndices.push_back(p_element->GetIndex());
        element_indices_iterator++;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::ReIndex(NodeMap& map)
{
    assert(!mAddedNodes);
    map.Resize(this->GetNumAllNodes());

    std::vector<Element<ELEMENT_DIM, SPACE_DIM> *> live_elements;
    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        if (this->mElements[i]->IsDeleted())
        {
            delete this->mElements[i];
        }
        else
        {
            live_elements.push_back(this->mElements[i]);
        }
    }

    assert (mDeletedElementIndices.size() == this->mElements.size()-live_elements.size());
    mDeletedElementIndices.clear();
    this->mElements = live_elements;

    std::vector<Node<SPACE_DIM> *> live_nodes;
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        if (this->mNodes[i]->IsDeleted())
        {
            delete this->mNodes[i];
            map.SetDeleted(i);
        }
        else
        {
            live_nodes.push_back(this->mNodes[i]);
            // the nodes will have their index set to be the index into the live_nodes
            // vector further down
            map.SetNewIndex(i, (unsigned)(live_nodes.size()-1));
        }
    }

    assert (mDeletedNodeIndices.size() == this->mNodes.size()-live_nodes.size());
    this->mNodes = live_nodes;
    mDeletedNodeIndices.clear();

    std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *> live_boundary_elements;
    for (unsigned i=0; i<this->mBoundaryElements.size(); i++)
    {
        if (this->mBoundaryElements[i]->IsDeleted())
        {
            delete this->mBoundaryElements[i];
        }
        else
        {
            live_boundary_elements.push_back(this->mBoundaryElements[i]);
        }
    }

    assert (mDeletedBoundaryElementIndices.size() == this->mBoundaryElements.size()-live_boundary_elements.size());
    this->mBoundaryElements = live_boundary_elements;
    mDeletedBoundaryElementIndices.clear();

    for (unsigned i=0; i<this->mNodes.size();i++)
    {
        this->mNodes[i]->SetIndex(i);
    }
    for (unsigned i=0; i<this->mElements.size();i++)
    {
        this->mElements[i]->ResetIndex(i);
    }
    for (unsigned i=0; i<this->mBoundaryElements.size();i++)
    {
        this->mBoundaryElements[i]->ResetIndex(i);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh(NodeMap& map)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert( SPACE_DIM==2 || SPACE_DIM==3 );
    assert( ELEMENT_DIM == SPACE_DIM );
    #undef COVERAGE_IGNORE
    
    // Avoid some triangle/tetgen errors: need at least four 
    // nodes for tetgen, and at least three for triangle
    if (GetNumNodes() <= SPACE_DIM)
    {
        EXCEPTION("The number of nodes must exceed the spatial dimension.");   
    }

// I think the following dramatically slows down the simulations, as with a 
// decent sized mesh we nearly always do need to remesh.

//    // If there are no nodes waiting to be deleted,
//    // and the current mesh is Voronoi, then we don't
//    // need to call triangle/tetgen
//    if (mDeletedNodeIndices.size()==0 && !mAddedNodes)
//    {        
//        if (CheckVoronoi())
//        {
//            map.ResetToIdentity();
//            return;
//        }
//    }
    
    // Make sure the map is big enough
    map.Resize(this->GetNumAllNodes());
    
    if (SPACE_DIM==2)  // In 2D, remesh using triangle via library calls
    {
        struct triangulateio triangle_input;
        triangle_input.pointlist = (double *) malloc(GetNumNodes() * 2 * sizeof(double));
        triangle_input.numberofpoints = GetNumNodes();
        triangle_input.numberofpointattributes = 0;
        triangle_input.pointmarkerlist = NULL;
        triangle_input.numberofsegments = 0;
        triangle_input.numberofholes = 0;
        triangle_input.numberofregions = 0;
    
        unsigned new_index = 0;
        for (unsigned i=0; i<this->GetNumAllNodes(); i++)
        {
            if (this->mNodes[i]->IsDeleted())
            {
                map.SetDeleted(i);
            }
            else
            {
                map.SetNewIndex(i, new_index);
                triangle_input.pointlist[2*new_index] = this->mNodes[i]->rGetLocation()[0];
                triangle_input.pointlist[2*new_index + 1] = this->mNodes[i]->rGetLocation()[1];
                new_index++;
            }
        }
    
        // Make structure for output
        struct triangulateio triangle_output;
        triangle_output.pointlist = NULL;
        triangle_output.pointattributelist = (double *) NULL;
        triangle_output.pointmarkerlist = (int *) NULL;
        triangle_output.trianglelist = (int *) NULL;
        triangle_output.triangleattributelist = (double *) NULL;
        triangle_output.edgelist = (int *) NULL;
        triangle_output.edgemarkerlist = (int *) NULL;
    
        // Library call
        triangulate((char*)"Qze", &triangle_input, &triangle_output, NULL);
    
        assert(triangle_output.numberofcorners == 3);
    
        // Remove current data
        Clear();
    
        // Construct the nodes
        for (unsigned node_index=0; node_index<(unsigned)triangle_output.numberofpoints; node_index++)
        {
            if (triangle_output.pointmarkerlist[node_index] == 1)
            {
                // Boundary node
                Node<SPACE_DIM> *p_node = new Node<SPACE_DIM>(node_index, true,
                  triangle_output.pointlist[node_index * 2],
                  triangle_output.pointlist[node_index * 2+1]);
                this->mNodes.push_back(p_node);
                this->mBoundaryNodes.push_back(p_node);
            }
            else
            {
                this->mNodes.push_back(new Node<SPACE_DIM>(node_index, false,
                  triangle_output.pointlist[node_index * 2],
                  triangle_output.pointlist[node_index * 2+1]));
            }
        }
    
        // Construct the elements
        this->mElements.reserve(triangle_output.numberoftriangles);
        for (unsigned element_index=0; element_index<(unsigned)triangle_output.numberoftriangles; element_index++)
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            for (unsigned j=0; j<3; j++)
            {
                unsigned global_node_index = triangle_output.trianglelist[element_index*3 + j];
                assert(global_node_index < this->mNodes.size());
                nodes.push_back(this->mNodes[global_node_index]);
            }
            this->mElements.push_back(new Element<ELEMENT_DIM, SPACE_DIM>(element_index, nodes));
        }
    
        // Construct the edges
        // too big mBoundaryElements.reserve(triangle_output.numberoftriangles);
        unsigned next_boundary_element_index = 0;
        for (unsigned boundary_element_index=0; boundary_element_index<(unsigned)triangle_output.numberofedges; boundary_element_index++)
        {
            if (triangle_output.edgemarkerlist[boundary_element_index] == 1)
            {
                std::vector<Node<SPACE_DIM>*> nodes;
                for (unsigned j=0; j<2; j++)
                {
                    unsigned global_node_index=triangle_output.edgelist[boundary_element_index*2 + j];
                    assert(global_node_index < this->mNodes.size());
                    nodes.push_back(this->mNodes[global_node_index]);
                }
                this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>(next_boundary_element_index++, nodes));
            }
        }
    
        free(triangle_input.pointlist);
    
        free(triangle_output.pointlist);
        free(triangle_output.pointattributelist);
        free(triangle_output.pointmarkerlist);
        free(triangle_output.trianglelist);
        free(triangle_output.triangleattributelist);
        free(triangle_output.edgelist);
        free(triangle_output.edgemarkerlist);
    }
    else // in 3D, remesh using tetgen
    {
        std::stringstream pid;
        pid << getpid();
    
        OutputFileHandler handler("");
        std::string full_name = handler.GetOutputDirectoryFullPath("") + "temp_" + pid.str() + ".";
    
        // Only the master process should do IO and call the mesher
        if (handler.IsMaster())
        {
            std::string node_file_name = "temp_" + pid.str() + ".node";
            {
                out_stream node_file = handler.OpenOutputFile(node_file_name);
    
                (*node_file) << GetNumNodes() << "\t" << SPACE_DIM << "\t0\t0\n";
    
                unsigned new_index = 0;
    
                for (unsigned i=0; i<this->GetNumAllNodes(); i++)
                {
                    if (this->mNodes[i]->IsDeleted())
                    {
                        map.SetDeleted(i);
                    }
                    else
                    {
                        map.SetNewIndex(i, new_index);
                        new_index++;
                        const c_vector<double, SPACE_DIM> node_loc = this->mNodes[i]->rGetLocation();
                        (*node_file) << i << "\t" << node_loc[0] << "\t" << node_loc[1] << "\t" << node_loc[2] << "\n";
                    }
                }
                node_file->close();
            }
    
            std::string binary_name = "tetgen";
            std::string command = "./bin/tetgen -Qe " + full_name + "node";
    
            // Tetgen's quiet mode isn't as quiet as Triangle's
            command += " > /dev/null";
            
            int return_value = system(command.c_str());
            if (return_value != 0)
            {
                NEVER_REACHED;
                //EXCEPTION("The tetgen mesher did not succeed in remeshing.");
            }
        }
        // Wait for the new mesh to be available and communicate its name
    #ifndef SPECIAL_SERIAL
        if (!PetscTools::IsSequential())
        {
            char full_name_comm[200];
            strcpy(full_name_comm, full_name.c_str());
            MPI_Bcast(full_name_comm, 200, MPI_CHAR, 0, MPI_COMM_WORLD);
            full_name = full_name_comm;
        }
    #endif //SPECIAL_SERIAL
    
        // Clear all current data
        Clear();
    
        // Read the new mesh back from file
        TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(full_name+"1");
        ConstructFromMeshReader(mesh_reader);
    
        // Make sure the file is not deleted before all the processors have read it
    #ifndef SPECIAL_SERIAL
        if (!PetscTools::IsSequential())
        {
            MPI_Barrier(PETSC_COMM_WORLD);
        }
    #endif //SPECIAL_SERIAL
    
        if (handler.IsMaster())
        {
            std::string remove_command = "rm " + full_name + "*";
            system(remove_command.c_str());
        }
    }    
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh()
{
    NodeMap map(GetNumNodes());
    ReMesh(map);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableMesh<ELEMENT_DIM, SPACE_DIM>::CheckVoronoi(Element<ELEMENT_DIM, SPACE_DIM> *pElement, double maxPenetration)
{
    assert (ELEMENT_DIM == SPACE_DIM);
    unsigned num_nodes = pElement->GetNumNodes();
    std::set<unsigned> neighbouring_elements_indices;
    std::set< Element<ELEMENT_DIM,SPACE_DIM> *> neighbouring_elements;
    std::set<unsigned> neighbouring_nodes_indices;

    //Form a set of neighbouring elements via the nodes
    for (unsigned i = 0; i < num_nodes; i++)
    {
        Node<SPACE_DIM>* node = pElement->GetNode(i);
        neighbouring_elements_indices = node->rGetContainingElementIndices();
        ///\todo Should use a set union operation here
        for (std::set<unsigned>::const_iterator it = neighbouring_elements_indices.begin();
                 it != neighbouring_elements_indices.end(); ++it)
            {
                neighbouring_elements.insert(this->GetElement(*it));
            }
    }
    neighbouring_elements.erase(pElement);

    //For each neighbouring element find the supporting nodes
    typedef typename std::set<Element<ELEMENT_DIM,SPACE_DIM> *>::const_iterator ElementIterator;

    for (ElementIterator it = neighbouring_elements.begin();
         it != neighbouring_elements.end(); ++it)
    {
        for (unsigned i = 0; i < num_nodes; i++)
        {
            neighbouring_nodes_indices.insert((*it)->GetNodeGlobalIndex(i));
        }
    }
    //Remove the nodes that support this element
    for (unsigned i = 0; i < num_nodes; i++)
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
        c_vector <double, ELEMENT_DIM> node_location = this->GetNode(*it)->rGetLocation();

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
bool MutableMesh<ELEMENT_DIM, SPACE_DIM>::CheckVoronoi(double maxPenetration)
{
    // Looping through all the elements in the mesh
    for (unsigned i=0; i < this->mElements.size();i++)
    {
        // Check if the element is not deleted
        if (!this->mElements[i]->IsDeleted())
        {
            // Checking the Voronoi of the Element
            if (CheckVoronoi(this->mElements[i], maxPenetration) == false)
            {
                return false;
            }
        }
    }
    return true;
}



#endif /*MUTABLEMESH_HPP_*/
