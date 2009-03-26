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

/**
 * A concrete mutable mesh class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MutableMesh : public TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>
{
protected:

    /**
     * Indices of elements that have been marked as deleted.
     * These indices can be reused when adding new elements.
     */
    std::vector<unsigned> mDeletedElementIndices;

    /**
     * Indices of boundary elements that have been marked as deleted.
     * These indices can be reused when adding new boundary elements.
     */
    std::vector<unsigned> mDeletedBoundaryElementIndices;

    /**
     * Indices of nodes that have been marked as deleted.
     * These indices can be reused when adding new nodes.
     */
    std::vector<unsigned> mDeletedNodeIndices;

    /** Whether any nodes have been added to the mesh. */
    bool mAddedNodes;

private:

    /**
     * Check whether any neighbouring node is inside the circumsphere of this element.
     * 
     * @param pointer to an element
     * @param maxPenetration is the maximum distance a node is allowed to be inside the
     * circumsphere of the element, as a proportion of the circumsphere radius.
     */
#define COVERAGE_IGNORE
    bool CheckVoronoi(Element<ELEMENT_DIM, SPACE_DIM>  *pElement, double maxPenetration);
#undef COVERAGE_IGNORE

public:

    /**
     * Constructor.
     */
    MutableMesh();

    /**
     * Constructor which takes in a vector of nodes.
     */
    MutableMesh(std::vector<Node<SPACE_DIM> *> nodes);

    /**
     * Destructor.
     */
    virtual ~MutableMesh();

    /**
     * Clear all the data in the mesh.
     */
    void Clear();

    /** 
     * Get the number of nodes that are actually in use.
     */
    unsigned GetNumNodes() const;

    /** 
     * Get the number of elements that are actually in use.
     */
    unsigned GetNumElements() const;

    /** 
     * Get the number of boundary elements that are actually in use.
     */
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

    /**
     * Move the node with a particular index to a new point in space and
     * verifies that the signed areas of the supporting Elements are positive.
     * 
     * @param index is the index of the node to be moved
     * @param point is the new target location of the node
     * @param concreteMove is set to false if we want to skip the signed area tests
     */
    virtual void SetNode(unsigned index, ChastePoint<SPACE_DIM> point, bool concreteMove=true);

    /**
     * Move one node to another (i.e. merges the nodes), refreshing/deleting 
     * elements as appropriate.
     *
     * @param index is the index of the node to be moved
     * @param targetIndex is the index of the node to move to
     * @param concreteMove can be set to false if you just want to check whether this will work.
     *     Set it to true if you're doing the merger for real, in order to do all the bookkeeping.
     */
    void MoveMergeNode(unsigned index, unsigned targetIndex, bool concreteMove=true);

#define COVERAGE_IGNORE
    /**
     * Delete a node from the mesh by finding an appropriate neighbour node
     * to merge it with.
     *
     * @param index is the index of the node to be deleted
     */
    void DeleteNode(unsigned index);
#undef COVERAGE_IGNORE

#define COVERAGE_IGNORE
    /**
     * Mark a node as deleted. Note that this method DOES NOT deal with the
     * associated elements and therefore should only be called immediately prior
     * to a ReMesh() being called. (Thus saves work compared to DeleteNode()
     * function and does not MoveMerge the node and elements).
     *
     * @param index The index of the node to delete
     */
    void DeleteNodePriorToReMesh(unsigned index);
#undef COVERAGE_IGNORE

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
#define COVERAGE_IGNORE     
    void ReIndex(NodeMap& map);
#undef COVERAGE_IGNORE

    /**
     * Re-mesh a mesh using triangle (via library calls) or tetgen
     * @param map is a NodeMap which associates the indices of nodes in the old mesh
     * with indices of nodes in the new mesh.  This should be created with the correct size (NumAllNodes)
     */
#define COVERAGE_IGNORE
    virtual void ReMesh(NodeMap& map);
#undef COVERAGE_IGNORE

    /**
     * Alternative version of remesh which takes no parameters does not require a NodeMap. Note: inherited
     * classes should overload ReMesh(NodeMap&)
     */
#define COVERAGE_IGNORE
    void ReMesh();
#undef COVERAGE_IGNORE
    
    /**
     * Checks the entire mesh element by element and checks whether any neighbouring node
     * is inside the circumsphere of this element.
     * @param maxPenetration is the maximum distance a node is allowed to be inside the
     * circumsphere of an element that it is not a member of, as a proportion of the
     * circumsphere radius.
     */
#define COVERAGE_IGNORE     
    bool CheckVoronoi(double maxPenetration=0.0);    
#undef COVERAGE_IGNORE
};


#endif /*MUTABLEMESH_HPP_*/
