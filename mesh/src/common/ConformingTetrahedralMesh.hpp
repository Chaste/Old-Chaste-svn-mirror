#ifndef _CONFORMINGTETRAHEDRALMESH_HPP_
#define _CONFORMINGTETRAHEDRALMESH_HPP_

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "AbstractMeshReader.hpp"
#include "TrianglesMeshReader.hpp"
#include "Element.hpp"
#include "BoundaryElement.hpp"
#include "Node.hpp"
#include "NodeMap.hpp"
#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
/**
 * \todo
 * Work still needs to be done with boundary nodes & elements?
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ConformingTetrahedralMesh
{
public:
    typedef typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> *>::const_iterator ElementIterator;
    typedef typename std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *>::const_iterator BoundaryElementIterator;
    typedef typename std::vector<Node<SPACE_DIM> *>::const_iterator BoundaryNodeIterator;
private:
    unsigned mNumCornerNodes;
    // Note that since these are vectors of objects, not pointers, push_back
    // will copy the objects.
    std::vector<Element<ELEMENT_DIM, SPACE_DIM> *> mElements;
    std::vector<Node<SPACE_DIM> *> mNodes;
    std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *> mBoundaryElements;
    
    /// Indices of elements/nodes that have been deleted - these indices can be reused when adding
    /// new elements/nodes
    std::vector<unsigned> mDeletedElementIndices;
    std::vector<unsigned> mDeletedBoundaryElementIndices;
    std::vector<unsigned> mDeletedNodeIndices;
    
    std::vector< Node<SPACE_DIM> *> mBoundaryNodes;
    
    //ElementIterator mpElementIter;
    //BoundaryElementIterator mpBoundaryElementIter;
    //BoundaryNodeIterator mpBoundaryNodeIter;
    
    unsigned AddNode(Node<SPACE_DIM> *pNewNode);
    
    /**
     * Check whether any neighbouring node is inside the circumsphere of this element.
     * @param pointer to an element
     * @param maxPenetration is the maximum distance a node is allowed to be inside the
     * circumsphere of the element, as a proportion of the circumsphere radius.
     */
      bool CheckVoronoi(Element<ELEMENT_DIM, SPACE_DIM>  *pElement, double maxPenetration);
    
public:

    ConformingTetrahedralMesh();
    ConformingTetrahedralMesh(unsigned numElements);
    
    ~ConformingTetrahedralMesh();
    
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader, unsigned orderOfBasisFunctions=1);
    
    void RescaleMeshFromBoundaryNode(Point<1> updatedPoint, unsigned boundaryNodeIndex);
    
    Node<SPACE_DIM> *GetNode(unsigned index);
    
    unsigned GetNumNodes();
    unsigned GetNumElements();
    unsigned GetNumBoundaryElements();
    unsigned GetNumAllNodes();
    unsigned GetNumAllElements();
    unsigned GetNumAllBoundaryElements();
    unsigned GetNumBoundaryNodes();
    unsigned GetNumCornerNodes();
    
    /**
     * Return a pointer to the first element in the mesh.
     */
    ElementIterator GetElementIteratorBegin() const
    {
        return mElements.begin();
    }
    /**
     * Return a pointer to *one past* the last element in the mesh
     * (for consistency with STL iterators).
     */
    ElementIterator GetElementIteratorEnd() const
    {
        return mElements.end();
    }
    
    /**
     * Return a pointer to the first boundary element in the mesh.
     */
    
    BoundaryElementIterator GetBoundaryElementIteratorBegin() const
    {
        return mBoundaryElements.begin();
    }
    /**
     * Return a pointer to *one past* the last boundary element in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryElementIterator GetBoundaryElementIteratorEnd() const
    {
        return mBoundaryElements.end();
    }
    
    /**
     * Return a pointer to the first boundary node in the mesh.
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorBegin() const
    {
        return mBoundaryNodes.begin();
    }
    /**
     * Return a pointer to *one past* the last boundary node in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorEnd() const
    {
        return mBoundaryNodes.end();
    }
    
    Element<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index)
    {
        return (mElements[index]);
    }
    
    BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* GetBoundaryElement(unsigned index)
    {
        return (mBoundaryElements[index]);
    }
    
    void SetNode(unsigned index, Point<SPACE_DIM> point, bool concreteMove=true);
    void MoveMergeNode(unsigned index, unsigned targetIndex, bool concreteMove=true);
    void DeleteNode(unsigned index);
    
    unsigned RefineElement(Element<ELEMENT_DIM,SPACE_DIM>* pElement, Point<SPACE_DIM> Point);
    void RefreshMesh(void);
    
    /**
     * Return the volume of a mesh, calculated by adding the determinant of each element 
     * and dividing by n!, where n is the element dimension.
     */
    double CalculateMeshVolume();
    double CalculateMeshSurface();
    
    void Translate(c_vector<double, SPACE_DIM> displacement);
    void Translate(const double xMovement=0.0, const double yMovement=0.0, const double zMovement=0.0);
    void Scale(const double xFactor=1.0, const double yFactor=1.0, const double zFactor=1.0);
    void Rotate(c_matrix<double , SPACE_DIM, SPACE_DIM> rotation_matrix);
    void Rotate(c_vector<double,3> axis, double angle);
    void RotateX(const double theta);
    void RotateY(const double theta);
    void RotateZ(const double theta);
    /**Rotating a 2D mesh equates that rotation around the z-axis*/
    void Rotate(double theta)
    {
        RotateZ(theta);
    }
    
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
    void ReIndex(void);
    /**
     * Re-mesh a mesh using triangle or tetgen
     * @param map is a NodeMap which associates the indices of nodes in the old mesh
     * with indices of nodes in the new mesh.  This should be created with the correct size (NumAllNodes)
     */
    void ReMesh(NodeMap &map);
    
    /**
     * Permute the nodes so that they appear in a different order in mNodes
     * (and their mIndex's are altered accordingly).
     * @param rRng is a random number generators instance (so that the caller can seed
     * the permutation
     */
    void PermuteNodes(RandomNumberGenerator &rRng);
   
   /**
     * Permute the nodes so that they appear in a different order in mNodes
     * (and their mIndex's are altered accordingly) using Metis binaries.
     * 
     */
    void PermuteNodesWithMetisBinaries();
   
    /**
      * Permute the nodes so that they appear in a different order in mNodes
      * (and their mIndex's are altered accordingly).
     * @param perm is a vector containing the new indices
     */
     void PermuteNodes(std::vector<unsigned> perm);
    
    /**    
     * Checks the entire mesh element by element and checkes whether any neighbouring node
     * inside the circumsphere of this element.
     * @param maxPenetration is the maximum distance a node is allowed to be inside the
     * circumsphere of an element that it is not a member of, as a proportion of the
     * circumsphere radius.
     */
    bool CheckVoronoi(double maxPenetration=0.0);
    
     /**
      * Construct a rectangular grid on [0,width]x[0,height]
      * diagonals can be staggered so that there is no prefered diffusion propagation
      * direction.
      */
     void ConstructRectangularMesh(unsigned width, unsigned height, bool stagger=true);
 
     /**
      * Construct a cuboid grid on [0,width]x[0,height]x[0,depth]
      * diagonals can be staggered so that there is no prefered diffusion propagation
      * direction.
      */
     void ConstructCuboid(unsigned width, unsigned height, unsigned depth);
     
     /**
      *  Returns the element index for the first element that is known to contain a test point
      *  @param testPoint
      *  @param strict Should the element returned contain the point in the interior and
      *  not on an edge/face/vertex (default = not strict)
      */
     unsigned GetContainingElementIndex(Point<SPACE_DIM> testPoint, bool strict=false);

     /**
      *  Returns the element index for an element is closest to the testPoint
      * "Closest" means that the minimum interpolation weights for the testPoint are 
      * maximised for this element
      *  @param testPoint
      * 
      */
     unsigned GetNearestElementIndex(Point<SPACE_DIM> testPoint);

     /**
      *  Returns all element indices for elements that are known to contain a test point
      *  @param testPoint 
      */
     std::vector<unsigned> GetContainingElementIndices(Point<SPACE_DIM> testPoint);


	/**
	 * Sets the ownership of each element according to which nodes are owned by the
	 * process.
	 * @param lo is the lowest node number owned by the process
	 * @param hi is one higher than the highest node number owned by the process
	 * ie. this process owns nodes [lo..hi)
	 * and element is "owned" if one or more of its nodes are owned
	 */
	 void SetElementOwnerships(unsigned lo, unsigned hi);
	 
};

#endif //_CONFORMINGTETRAHEDRALMESH_HPP_
