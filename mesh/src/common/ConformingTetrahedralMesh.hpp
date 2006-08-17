#ifndef _CONFORMINGTETRAHEDRALMESH_HPP_
#define _CONFORMINGTETRAHEDRALMESH_HPP_

#include <vector>
#include "AbstractMeshReader.hpp"
#include "Element.cpp"
#include "BoundaryElement.cpp"
#include "Node.hpp"

/**
 * \todo
 * Work still needs to be done with boundary nodes & elements?
 */

template<int ELEMENT_DIM, int SPACE_DIM>
class ConformingTetrahedralMesh
{
public:
    typedef typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> *>::const_iterator ElementIterator;
    typedef typename std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *>::const_iterator BoundaryElementIterator;
    typedef typename std::vector<const Node<SPACE_DIM> *>::const_iterator BoundaryNodeIterator;
private:
    int mNumCornerNodes;
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
    
    // \todo change to indices?
    std::vector<const Node<SPACE_DIM> *> mBoundaryNodes;
    
    ElementIterator mpElementIter;
    BoundaryElementIterator mpBoundaryElementIter;
    BoundaryNodeIterator mpBoundaryNodeIter;
    
    int AddNode(Node<SPACE_DIM> *pNewNode);
    
    
public:

    ConformingTetrahedralMesh();
    ConformingTetrahedralMesh(long numElements);
    
    ~ConformingTetrahedralMesh();
    
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader, int orderOfBasisFunctions=1);
    
    void RescaleMeshFromBoundaryNode(Point<1> updatedPoint, int boundaryNodeIndex);
    
    Node<SPACE_DIM> *GetNodeAt(long index);
    
    long GetNumNodes();
    long GetNumElements();
    long GetNumBoundaryElements();
    long GetNumAllNodes();
    long GetNumAllElements();
    long GetNumAllBoundaryElements();
    long GetNumBoundaryNodes();
    long GetNumCornerNodes();
    
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
    
    Element<ELEMENT_DIM, SPACE_DIM>* GetElement(int index)
    {
        return (mElements[index]);
    }
    
    BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* GetBoundaryElement(int index)
    {
        return (mBoundaryElements[index]);
    }
    
    void SetNode(unsigned index, Point<SPACE_DIM> point, bool verify=true);
    
    int RefineElement(Element<ELEMENT_DIM,SPACE_DIM>* pElement, Point<SPACE_DIM> Point);
    void RefreshMesh(void);
    
    /**
     * Return the volume of a mesh, calculated by adding the determinant of each element 
     * and dividing by n!, where n is the element dimension.
     */
    double CalculateMeshVolume();
    double CalculateMeshSurface();
    
    void Translate(c_vector<double, SPACE_DIM> displacement);
    void Translate(const double xMovement=0.0, const double yMovement=0.0, const double zMovement=0.0);
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
    void DeleteBoundaryNodeAt(long index);
};

#endif //_CONFORMINGTETRAHEDRALMESH_HPP_
