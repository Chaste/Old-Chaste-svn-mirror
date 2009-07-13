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
#ifndef VERTEXMESH_HPP_
#define VERTEXMESH_HPP_

#include <iostream>
#include <map>
#include <algorithm>

#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractTetrahedralMesh.hpp"
#include "VertexMeshReader.hpp"
#include "VertexElement.hpp"
#include "VertexElementMap.hpp"

/**
 * A vertex-based mesh class, for use in vertex-based tissue simulations.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMesh : public AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestVertexMesh;

protected:

    /** Vector of pointers to VertexElements. */
    std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> mElements;

    /** The minimum distance apart that two nodes in the mesh can be without causing element rearrangment. */
    double mCellRearrangementThreshold;

    /** The maximum distance apart that neighbouring nodes in the mesh can be without the edge being divided. */
    double mEdgeDivisionThreshold;
    
    /** The area threshold at which T2 swaps occur in an apoptotic, triangular cell/element */
    double mT2Threshold;

    /** Indices of nodes that have been deleted. These indices can be reused when adding new elements/nodes. */
    std::vector<unsigned> mDeletedNodeIndices;

    /** Indices of elements that have been deleted. These indices can be reused when adding new elements. */
    std::vector<unsigned> mDeletedElementIndices;

    /** Whether nodes have been added to the mesh. */
    bool mAddedNodes;

    /** Whether Elements have been added to the mesh. */
    bool mAddedElements;

    /** Create correspondences between VertexElements and Nodes in the mesh. */
    void SetupVertexElementsOwnedByNodes();

    /**
     * Helper method for ReMesh to Identify the type of swap
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param pNodeA one of the nodes to perform the swap with
     * @param pNodeB the other node to perform the swap
     */
    void IdentifySwapType(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);

    /**
     * Helper method for ReMesh to merge nodes when needed.
     * Replaces the node contained in the least number of elements with the other node.
     *
     * @param pNodeA one of the nodes to perform the merge with
     * @param pNodeB the other node to perform the merge with 
     */
    void PerformNodeMerge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);

    /**
     * Helper method for ReMesh to perform the T1 Swap
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param pNodeA one of the nodes to perform the swap with
     * @param pNodeB the other node to perform the swap
     * @param elementsContainingNodes set of common elements
     */
    void PerformT1Swap(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, std::set<unsigned> elementsContainingNodes);

    /**
     * Method to divide an element given 2 nodes in which to divide the element with
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param pElement the element to divide
     * @param nodeAIndex the local index of node where to divide
     * @param nodeBIndex the local index of node where to divide
     *
     * @return the index of the new element
     */
    unsigned DivideElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned nodeAIndex, unsigned nodeBIndex);

    /**
     * Test whether a given point lies inside a given element.
     *
     * We use a ray-casting algorithm, which relies on the following result:
     * if the point in question is not on the boundary of the element, then
     * the number of intersections is an even number if the point is outside,
     * and it is odd if inside.
     *
     * Currently the method is coded 'strictly', such that points lying on
     * an edge or at a vertex are considered to lie outside the element.
     *
     * @param testPoint the point to test
     * @param elementIndex global index of the element in the mesh
     *
     * @return if the point is included in the element.
     */
    bool ElementIncludesPoint(const c_vector<double, SPACE_DIM>& testPoint, unsigned elementIndex);

    /**
     * Get the local index of a given element which is the start vertex of the edge
     * of the element that the overlapping point testPoint is closest to.
     *
     * @param testPoint the point to test
     * @param elementIndex global index of the element in the mesh
     *
     * @return the local index
     */
    unsigned GetLocalIndexForElementEdgeClosestToPoint(const c_vector<double, SPACE_DIM>& testPoint, unsigned elementIndex);

    /**
     * Called by ReMesh(). Moves a node, which has been found to overlap an element,
     * back onto the edge of that element and associates it with the element.
     *
     * @param pNode pointer to the node
     * @param elementIndex global index of the element in the mesh
     */
    void MoveOverlappingNodeOntoEdgeOfElement(Node<SPACE_DIM>* pNode, unsigned elementIndex);

    /**
     * Solve node mapping method. This overridden method is required 
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the node
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Solve element mapping method. This overridden method is required 
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the element
     */
    unsigned SolveElementMapping(unsigned index) const;

    /**
     * Solve boundary element mapping method. This overridden method is required 
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the boundary element
     */
    unsigned SolveBoundaryElementMapping(unsigned index) const;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archives the member variables of the object which
     * have to be preserved during its lifetime.
     *
     * The remaining member variables are re-initialised before being used
     * by each ReMesh() call so they do not need to be archived.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mCellRearrangementThreshold;
        archive & mEdgeDivisionThreshold;
    }

public:

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param vertexElements vector of pointers to VertexElements
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangment (defaults to 0.01)
     * @param edgeDivisionThreshold the maximum threshold distance for edge division (defaults to 1.5)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     */
    VertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
               std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements,
               double cellRearrangementThreshold=0.01,
               double edgeDivisionThreshold=1.5,
               double t2Threshold=0.001);

    /**
     * Helper constructor, creates a rectangular vertex-based mesh.
     *
     * @param numAcross number of VertexElements across
     * @param numUp number of VertexElements up
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangment
     * @param edgeDivisionThreshold the maximum threshold distance for edge division
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.01)
     */
    VertexMesh(unsigned numAcross, unsigned numUp, double cellRearrangementThreshold, double edgeDivisionThreshold, double t2Threshold = 0.01);

    /**
     * Constructor for use by serializer.
     *
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangment (defaults to 0.01)
     * @param edgeDivisionThreshold the maximum threshold distance for edge division (defaults to 1.5)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.01)
     */
    VertexMesh(double cellRearrangementThreshold=0.01, double edgeDivisionThreshold=1.5, double t2Threshold=0.01);

    /**
     * Destructor.
     */
    virtual ~VertexMesh();

    /**
     * Set method for mCellRearrangementThreshold.
     *
     * @param cellRearrangementThreshold
     */
    void SetCellRearrangementThreshold(double cellRearrangementThreshold);

    /**
     * Set method for mEdgeDivisionThreshold.
     *
     * @param edgeDivisionThreshold
     */
    void SetEdgeDivisionThreshold(double edgeDivisionThreshold);

    /**
     * Set method for mT2Threshold.
     * 
     * @param t2Threshold
     */
    void SetT2Threshold(double t2Threshold);

    /**
     *  Move the node with a particular index to a new point in space.
     *
      * @param nodeIndex the index of the node to be moved
      * @param point the new target location of the node
      */
    void SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point);

    /**
     * @return mCellRearrangementThreshold
     */
    double GetCellRearrangementThreshold() const;

    /**
     * @return mEdgeDivisionThreshold
     */
    double GetEdgeDivisionThreshold() const;

    /**
     * @return mT2Threshold
     */
    double GetT2Threshold() const;

    /**
     * @return the number of Nodes in the mesh.
     */
    unsigned GetNumNodes() const;

    /**
     * @return the number of VertexElements in the mesh.
     */
    unsigned GetNumElements() const;

    /**
     * @return the number of VertexElements in the mesh, including those marked as deleted.
     */
    unsigned GetNumAllElements();

    /**
     * @param index  the global index of a specified vertex element
     *
     * @return a pointer to the vertex element
     */
    VertexElement<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

    /**
     * Compute the area of an element.
     *
     * This needs to be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the area of the element
     */
    virtual double GetAreaOfElement(unsigned index);

    /**
     * Compute the perimeter of an element.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the perimeter of the element
     */
    double GetPerimeterOfElement(unsigned index);

    /**
     * Compute the centroid of an element.
     *
     * This needs to be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (centroid_x,centroid_y).
     */
    virtual c_vector<double, SPACE_DIM> GetCentroidOfElement(unsigned index);

    /**
     * Compute the area gradient of an element at one of its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the area of the element, evaluated at this node.
     */
    c_vector<double, SPACE_DIM> GetAreaGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the gradient of the edge of an element ending at its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the edge of the element that ends at this node.
     */
    c_vector<double, SPACE_DIM> GetPreviousEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the gradient of the edge of an element starting at its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the edge of the element that starts at this node.
     */
    c_vector<double, SPACE_DIM> GetNextEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the gradient of the perimeter of an element at its nodes.
     * This returns the sum of GetPreviousEdgeGradientAtNode() and GetNextEdgeGradientAtNode().
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the perimeter of the element, evaluated at this node.
     */
    c_vector<double, SPACE_DIM> GetPerimeterGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the second moments of area of a given (polygonal) element.
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (Ixx,Iyy,Ixy).
     */
    virtual c_vector<double, 3> CalculateMomentsOfElement(unsigned index);

    /**
     * Calculate the vector of the shortest axis of a given element.
     * This is the eigenvector associated with the largest eigenvalue
     * of the inertial tensor. If the polygon is regular then the
     * eigenvalues are the same, so we return a random unit vector.
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (short_axis_x, short_axis_y).
     */
    c_vector<double, SPACE_DIM> GetShortAxisOfElement(unsigned index);

    /**
     * Given a node, find a set containing the indices of its neighbouring nodes.
     *
     * @param nodeIndex global index of the node
     * @return its neighbouring node indices
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);

    /**
     * Given a node and one of its containing elements, find a set containing
     * the indices of those neighbouring node(s) that are NOT also in the element.
     *
     * Note that we allow for more than one such index, since there is no reason
     * a priori to assume that each node is contained by exactly three elements.
     *
     * @param nodeIndex global index of the node
     * @param elemIndex global index of the element
     *
     * @return its neighbouring nodes that are not in the element
     */
    std::set<unsigned> GetNeighbouringNodeNotAlsoInElement(unsigned nodeIndex, unsigned elemIndex);

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     * @param cullInternalFaces whether to cull internal faces (defaults to false)
     */
    void ConstructFromMeshReader(AbstractTetrahedralMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader,
                                 bool cullInternalFaces=false);

    /**
     * Add a node to the mesh.
     *
     * Note: After calling this one or more times, you must then call ReMesh.
     *
     * @param pNewNode pointer to the new node
     * @return the global index of the new node in the mesh.
     */
    unsigned AddNode(Node<SPACE_DIM>* pNewNode);

    /**
     * Helper method for ReMesh to perform the T2 Swap
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param pElement is the element to remove
     */
    void PerformT2Swap(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement);          

    /**
     * Perform a T2 swap on an apoptotic cell. Checks that it is both a triangle
     * and small enough first.
     * 
     * @param pElement Pointer to the element to be checked
     */
    void PerformT2SwapIfNecessary(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement);

    /**
     * Mark an element as deleted. Note that it DOES NOT deal with the associated
     * nodes and therefore should only be called immediately prior to a ReMesh()
     * being called.
     *
     * @param index  the global index of a specified vertex element
     */
    void DeleteElementPriorToReMesh(unsigned index);

    /**
     * Method to divide an element in half
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param pElement the element to divide
     *
     * @return the index of the new element
     */
    unsigned DivideElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement);

    /**
     * Method to divide an element in half using a specific axis
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param pElement the element to divide
     * @param AxisOfDivision axis to divide the element by 
     * 
     * @return the index of the new element
     */
    unsigned DivideElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, c_vector<double, SPACE_DIM> AxisOfDivision);


    /**
     * Add an element to the mesh.
     *
     * @param pNewElement the new element
     *
     * @return the index of the new element in the mesh
     */
    unsigned AddElement(VertexElement<ELEMENT_DIM, SPACE_DIM>* pNewElement);

    /**
     * Delete mNodes and mElements.
     */
    void Clear();

    /**
     * Add a node on the edge between two nodes.
     *
     * @param pNodeA a pointer to one node
     * @param pNodeB a pointer to the other nodes
     */
    void DivideEdge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);

    /**
     * Re-mesh the mesh.
     *
     * @param elementMap a VertexElementMap which associates the indices of VertexElements in the old mesh
     *                   with indices of VertexElements in the new mesh.  This should be created
     *                   with the correct size, GetNumElements()
     */
    void ReMesh(VertexElementMap& elementMap);

    /**
     * Alternative version of remesh which takes no parameters does not require a VertexElementMap.
     * Note: inherited classes should overload ReMesh(VertexElementMap&).
     */
    void ReMesh();

};

#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS2(VertexMesh, 1, 1)
EXPORT_TEMPLATE_CLASS2(VertexMesh, 1, 2)
EXPORT_TEMPLATE_CLASS2(VertexMesh, 2, 2)
EXPORT_TEMPLATE_CLASS2(VertexMesh, 2, 3)
EXPORT_TEMPLATE_CLASS2(VertexMesh, 3, 3)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a VertexMesh.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const VertexMesh<ELEMENT_DIM, SPACE_DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const double cell_rearrangement_threshold = t->GetCellRearrangementThreshold();
    ar << cell_rearrangement_threshold;

    const double edge_division_threshold = t->GetEdgeDivisionThreshold();
    ar << edge_division_threshold;
}

/**
 * De-serialize constructor parameters and initialise a VertexMesh.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, VertexMesh<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    double cell_rearrangement_threshold;
    ar >> cell_rearrangement_threshold;

    double edge_division_threshold;
    ar >> edge_division_threshold;

    // Invoke inplace constructor to initialise instance
    ::new(t)VertexMesh<ELEMENT_DIM, SPACE_DIM>(cell_rearrangement_threshold, edge_division_threshold);
}
}
} // namespace ...

#endif /*VERTEXMESH_HPP_*/
