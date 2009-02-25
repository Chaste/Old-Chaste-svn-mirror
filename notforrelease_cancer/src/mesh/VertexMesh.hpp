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

#include "AbstractMesh.hpp"
#include "VertexMeshReader2d.hpp"
#include "VertexElement.hpp"
#include "NodeMap.hpp"

/**
 * A vertex-based mesh class, for use in vertex-based tissue simulations.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMesh
{
    friend class TestVertexMesh;
    
protected:

    /** Vector of pointers to nodes. */
    std::vector<Node<SPACE_DIM>*> mNodes;

    /** Vector of pointers to VertexElements. */    
    std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> mElements;
   
    /** The minimum distance apart that two nodes in the mesh can be without causing element rearrangment. */
    double mCellRearrangementThreshold;

    /** The maximum distance apart that neighbouring nodes in the mesh can be without the edge being divided. */ 
    double mEdgeDivisionThreshold;

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
     * @param rNodeA one of the nodes to perform the swap with 
     * @param rNodeB the other node to perform the swap
     */  
    void IdentifySwapType(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);
    
    /**
     * Helper method for ReMesh to merge nodes when needed 
     * Move node with smallest global index to center and remove other node
     * 
     * @param rNodeA one of the nodes to perform the swap with 
     * @param rNodeB the other node to perform the swap
     */  
    void PerformNodeMerge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB,
                          std::set<unsigned> ElementsContainingNodes);
    
    /**
     * Helper method for ReMesh to perform the T1 Swap
     * 
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     * 
     * @param rNodeA one of the nodes to perform the swap with 
     * @param rNodeB the other node to perform the swap 
     */  
    void PerformT1Swap(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, 
                       std::set<unsigned> ElementsContainingNodes);
    
    /**
     * Method to divide an element given 2 nodes in which to divide the element with 
     * 
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     * 
     * @param rElement the element to divide
     * @param nodeAIndex the local index of node where to divide
     * @param nodeBindex the local index of node where to divide
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

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
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
     */
    VertexMesh(std::vector<Node<SPACE_DIM>*> nodes, 
               std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements,
               double cellRearrangementThreshold=0.01,
               double edgeDivisionThreshold=1.5);

    /**
     * Helper constructor, creates a rectangular vertex-based mesh.
     * 
     * @param numAcross number of VertexElements across
     * @param numUp number of VertexElements up
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangment
     * @param edgeDivisionThreshold the maximum threshold distance for edge division
     */
    VertexMesh(unsigned numAcross, unsigned numUp, double cellRearrangementThreshold, double edgeDivisionThreshold);

    /**
     * Constructor for use by serializer.
     *
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangment (defaults to 0.01)
     * @param edgeDivisionThreshold the maximum threshold distance for edge division (defaults to 1.5)
     */
    VertexMesh(double cellRearrangementThreshold=0.01, double edgeDivisionThreshold=1.5);

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
     * Calculates the `width' of any dimension of the mesh.
     *
     * @param rDimension a dimension (0,1 or 2)
     * @return The maximum distance between any nodes in this dimension.
     *
     * N.B. Overwritten in Cylindrical2dVertexMesh.
     */
    virtual double GetWidth(const unsigned& rDimension);

    /**
     * Calculates the `width extremes' of any dimension of the mesh.
     * 
     * @param rDimension a dimension (0,1 or 2)
     * @return The minimum and maximum co-ordinates of any node in this dimension.
     */
    c_vector<double,2> GetWidthExtremes(const unsigned& rDimension);

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
     * @param index the global index of a specified node
     * 
     * @return a pointer to the node
     */
    Node<SPACE_DIM>* GetNode(unsigned index) const;
    
    /**
     * @param index  the global index of a specified vertex element
     * 
     * @return a pointer to the vertex element
     */    
    VertexElement<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

    /*
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

    /*
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
     * Returns distance between two nodes
     *
     * @param indexA a node index
     * @param indexB a node index
     *
     * @return straight line distance between two nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden 
     * in daughter classes for non-Euclidean metrics.
     */
    double GetDistanceBetweenNodes(unsigned indexA, unsigned indexB);

    /**
     * Returns a vector between two points in space
     *
     * @param rLocationA a c_vector of co-ordinates
     * @param rLocationB a c_vector of co-ordinates
     *
     * @return vector from location A to location B.
     *
     * N.B. This can be overridden in daughter classes for non-Euclidean metrics.
     */
    virtual c_vector<double, SPACE_DIM> GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocationA, const c_vector<double, SPACE_DIM>& rLocationB);

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
     * Construct the mesh using a mesh reader.
     * 
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(VertexMeshReader2d& rMeshReader);

    /**
     * Scale the mesh.
     * 
     * @param xFactor is the scale in the x-direction,
     * @param yFactor is the scale in the y-direction,
     * @param zFactor is the scale in the z-direction
     **/
    void Scale(const double xFactor=1.0, const double yFactor=1.0, const double zFactor=1.0);

    /**
     * Add a node to the mesh.
     *
     * Note: After calling this one or more times, you must then call ReMesh.
     * 
     * @param pNewNode pointer to the new node
     * @return the global index of the new node in the mesh.
     */
    unsigned AddNode(Node<SPACE_DIM> *pNewNode);

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
     * Add an element to the mesh.
     */
    unsigned AddElement(VertexElement<ELEMENT_DIM, SPACE_DIM> *pNewElement);

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
     * @param elementMap a NodeMap which associates the indices of VertexElements in the old mesh
     *                   with indices of VertexElements in the new mesh.  This should be created 
     *                   with the correct size, GetNumElements()
     */
    void ReMesh(NodeMap& elementMap);

    /**
     * Alternative version of remesh which takes no parameters does not require a NodeMap. 
     * Note: inherited classes should overload ReMesh(NodeMap&).
     */
    void ReMesh();

};

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a VertexMesh
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
 * De-serialize constructor parameters and initialise VertexMesh.
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
