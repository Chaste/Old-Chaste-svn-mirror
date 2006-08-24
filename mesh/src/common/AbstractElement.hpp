#ifndef _ABSTRACTELEMENT_HPP_
#define _ABSTRACTELEMENT_HPP_
#include <float.h>
/**
 * This class defines an Element for use in FEM.
 */

#include "Node.hpp"
#include "Point.hpp"
#include "UblasCustomFunctions.hpp"

#include "Exception.hpp"

#include <vector>
#include <cmath>

// When creating an element within a mesh one needs to specify its global index
// If the element is not used within a mesh the following
// constant is used instead.
const unsigned INDEX_IS_NOT_USED=0;


template <int ELEMENT_DIM, int SPACE_DIM>
class AbstractElement
{
protected:
   unsigned mIndex;
    std::vector<Node<SPACE_DIM>*> mNodes;
    int mOrderOfBasisFunctions;
    bool mIsDeleted;
    c_matrix<double, SPACE_DIM, SPACE_DIM> *mpJacobian;
    c_matrix<double, SPACE_DIM, SPACE_DIM> *mpInverseJacobian;
    double mJacobianDeterminant;
    bool mOwnership, mOwnershipSet;
     
    void RefreshJacobian(void);
    
    
    /**
     * Method that constructs the element. This is required because of having
     * two different copy constructors, one with a new index and another without
     */
    void CommonConstructor(const AbstractElement &element)
    {
        //Note that the index must be already set by the calling constructor
        mIsDeleted = element.mIsDeleted;
        mNodes = element.mNodes;
        // Allow nodes to keep track of containing elements (but not surface/boundary elements)
        // Only done in copy constructor, since that is what is called to put elements
        // in the vector contained in ConformingTetrahedralMesh.
        
        mJacobianDeterminant = element.mJacobianDeterminant;
        mpJacobian = NULL;
        if (element.mpJacobian != NULL)
        {
            mpJacobian = new c_matrix<double, SPACE_DIM, SPACE_DIM>;
            *mpJacobian = *(element.mpJacobian);
        }
        mpInverseJacobian = NULL;
        if (element.mpInverseJacobian != NULL)
        {
            mpInverseJacobian = new c_matrix<double, SPACE_DIM, SPACE_DIM>;
            *mpInverseJacobian = *(element.mpInverseJacobian);
        }
    }
    
    /**
     * Used by the destructor and the assignment operator
     */
    void FreeMemory()
    {
    
        if (mpJacobian != NULL)
        {
            delete mpJacobian;
        }
        if (mpInverseJacobian != NULL)
        {
            delete mpInverseJacobian;
        }
    }
    
public:
    static const int NUM_CORNER_NODES = ELEMENT_DIM+1;
    
    virtual void RegisterWithNodes()=0;
    
    ///Main constructor
    AbstractElement(unsigned index, std::vector<Node<SPACE_DIM>*> nodes, int orderOfBasisFunctions=1);
    
    /**
     * Copy constructor. This is needed so that copies of an element don't
     * share pointers to the same matrices, which causes problems when copies
     * get destroyed.
     */
    AbstractElement(const AbstractElement &element)
    {
        mIndex=element.mIndex;
        CommonConstructor(element);
    }
    
    
    AbstractElement()
    {
    }
    
    /**
     * Element assignment - make this element equal to the other one.
     */
    virtual void operator=(const AbstractElement &element)
    {
        // First free memory held by us
        FreeMemory();
        
        // Now copy stuff
        mIndex=element.mIndex;
        CommonConstructor(element);
    }
    
    /**
     * Free memory potentially allocated in the constructor (or elsewhere) for
     * holding Jacobian and its inverse.
     */
    virtual ~AbstractElement()
    {
        FreeMemory();
    }
    
    void RefreshJacobianDeterminant(void);
    
    void AddInternalNode(const Node<SPACE_DIM>* internalNodeToAdd)
    {
        assert(mOrderOfBasisFunctions > 1);
        assert(mNodes.size() - NUM_CORNER_NODES < 0.5*SPACE_DIM*(SPACE_DIM+1));
        
        mNodes.push_back(internalNodeToAdd);
    }
    
    double GetNodeLocation(int localIndex, int dimension) const
    {
        assert(dimension < SPACE_DIM);
        assert((unsigned)localIndex < mNodes.size());
        return mNodes[localIndex]->rGetPoint()[dimension];
    }
    
    // note: this used to return a reference to a c_vector, in which case a
    // weird error arose where it compiled, ran and passed on some machines
    // but failed the tests (bad_size errors) on another machine.
    c_vector<double, SPACE_DIM> GetNodeLocation(int localIndex) const
    {
        assert((unsigned)localIndex < mNodes.size());
        Point<SPACE_DIM> point=mNodes[localIndex]->rGetPoint();
        return point.rGetLocation();
    }
    
    long GetNodeGlobalIndex(int localIndex) const
    {
        assert((unsigned)localIndex < mNodes.size());
        return mNodes[localIndex]->GetIndex();
    }
    
    Node<SPACE_DIM>* GetNode(int localIndex) const
    {
        assert((unsigned)localIndex < mNodes.size());
        return mNodes[localIndex];
    }
    
    int GetNumNodes() const
    {
        return mNodes.size();
    }
    
    
    void AddNode(Node<SPACE_DIM>* node)
    {
        mNodes.push_back(node);
    }
    
    const c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> *GetJacobian(void) const
    {
        return mpJacobian;
    }
    const c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> *GetInverseJacobian(void) const
    {
        return mpInverseJacobian;
    }
    double GetJacobianDeterminant(void) const
    {
        return mJacobianDeterminant;
    }
    
    /** Get the index of this element
     */
    const unsigned& GetIndex(void) const
    {
        return mIndex;
    }
    
    /** Update node at the given index
     *  @param rIndex is an local index to which node to change
     *  @param pNode is a pointer to the replacement node
     */
    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)
    {
        assert(rIndex < mNodes.size());
        
        // Remove it from the node at this location
        mNodes[rIndex]->RemoveElement(mIndex);
        
        // Update the node at this location
        mNodes[rIndex] = pNode;
        
        // Add element to this node
        mNodes[rIndex]->AddElement(mIndex);
    }
    
    /**
    * Mark an element as having been removed from the mesh.
    * Also notify nodes in the element that it has been removed.
    */
    virtual void MarkAsDeleted()=0;
    
    bool IsDeleted()
    {
        return mIsDeleted;
    }
    
    bool GetOwnership()
    {
        return mOwnership;
    }
        
    void SetOwnership(bool ownership)
    {
        mOwnership=ownership;
        mOwnershipSet=true;
    }
    
    bool GetOwnershipSet()
    {
        return mOwnershipSet;
    }
    
    
    
    
    
    
};


#endif //_ABSTRACTELEMENT_HPP_
