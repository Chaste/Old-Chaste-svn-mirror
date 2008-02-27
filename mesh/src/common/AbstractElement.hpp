#ifndef _ABSTRACTELEMENT_HPP_
#define _ABSTRACTELEMENT_HPP_
/**
 * This class defines an Element for use in FEM.
 */

#include "Node.hpp"
#include "ChastePoint.hpp"
#include "UblasCustomFunctions.hpp"

#include "Exception.hpp"

#include <vector>
#include <cmath>

// When creating an element within a mesh one needs to specify its global index
// If the element is not used within a mesh the following
// constant is used instead.
const unsigned INDEX_IS_NOT_USED=0;


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractElement
{
protected:
    unsigned mIndex;
    std::vector<Node<SPACE_DIM>*> mNodes;
    c_matrix<double, SPACE_DIM, SPACE_DIM> mJacobian;
    c_matrix<double, SPACE_DIM, SPACE_DIM> mInverseJacobian;
    c_vector<double, SPACE_DIM> mWeightedDirection; //Holds an area-weighted normal or direction.  Only used when ELEMENT_DIM < SPACE_DIM
    double mJacobianDeterminant;
    bool mIsDeleted;
    bool mOwnership;
    bool mFlag;
    
    
    
    /**
     * Method that constructs the element. This is required because of having
     * two different copy constructors, one with a new index and another without
     */
    void CommonConstructor(const AbstractElement &element)
    {
        //Note that the index must be already set by the calling constructor
        mNodes = element.mNodes;
        // Allow nodes to keep track of containing elements (but not surface/boundary elements)
        // Only done in copy constructor, since that is what is called to put elements
        // in the vector contained in ConformingTetrahedralMesh.
        
        mJacobianDeterminant = element.mJacobianDeterminant;
        mJacobian = element.mJacobian;
        mInverseJacobian = element.mInverseJacobian;
        mWeightedDirection = element.mWeightedDirection;
        
        // Copy various flags
        mIsDeleted = element.mIsDeleted;
        mOwnership = element.mOwnership;
        mFlag = element.mFlag;
    }
    
    
public:
    static const unsigned NUM_CORNER_NODES = ELEMENT_DIM+1;
    
    virtual void RegisterWithNodes()=0;
    
    ///Main constructor
    AbstractElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes);
    
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
    
    /**
     * \todo Why does the default constructor not do anything?
     */
    AbstractElement()
    {}
    
    /**
     * Element assignment - make this element equal to the other one.
     */
    virtual AbstractElement& operator=(const AbstractElement &element)
    {
        // Now copy stuff
        mIndex=element.mIndex;
        CommonConstructor(element);
        return *this;
    }
    
    virtual ~AbstractElement()
    {}
    
    void RefreshJacobianDeterminant(bool concreteMove=true);
    void ZeroJacobianDeterminant(void);
    void ZeroWeightedDirection(void);
    
    double GetNodeLocation(unsigned localIndex, unsigned dimension) const
    {
        assert(dimension < SPACE_DIM);
        assert((unsigned)localIndex < mNodes.size());
        return mNodes[localIndex]->rGetLocation()[dimension];
    }
    
    /**
     * \todo this used to return a reference to a c_vector, in which case a
     * weird error arose where it compiled, ran and passed on some machines
     * but failed the tests (bad_size errors) on another machine.
     */
    c_vector<double, SPACE_DIM> GetNodeLocation(unsigned localIndex) const
    {
        assert((unsigned)localIndex < mNodes.size());
        return mNodes[localIndex]->rGetLocation();
    }
    
    c_vector<double, SPACE_DIM> CalculateCentroid() const
    {
        c_vector<double, SPACE_DIM> centroid=zero_vector<double>(SPACE_DIM);
        for (unsigned i=0; i<=ELEMENT_DIM; i++)
        {
            centroid += mNodes[i]->rGetLocation();
        }
        return centroid/((double)(ELEMENT_DIM + 1));
    }
    
    unsigned GetNodeGlobalIndex(unsigned localIndex) const
    {
        assert((unsigned)localIndex < mNodes.size());
        return mNodes[localIndex]->GetIndex();
    }
    
    Node<SPACE_DIM>* GetNode(unsigned localIndex) const
    {
        assert((unsigned)localIndex < mNodes.size());
        return mNodes[localIndex];
    }
    
    unsigned GetNumNodes() const
    {
        return mNodes.size();
    }
    
    
    void AddNode(Node<SPACE_DIM>* node)
    {
        mNodes.push_back(node);
    }
    
    const c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> *GetJacobian(void) const
    {
        return &mJacobian;
    }
    const c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> *GetInverseJacobian(void) const
    {
        return &mInverseJacobian;
    }
    double GetJacobianDeterminant(void) const
    {
        return mJacobianDeterminant;
    }
    
    /** Get the volume of an element (or area in 2d, or length in 1d) */
    double GetVolume(void) const
    {
        assert(SPACE_DIM == ELEMENT_DIM);
        double scale_factor = 1.0;
        
        if (ELEMENT_DIM == 2)
        {
            scale_factor = 2.0;  // both the volume of the canonical triangle is 0.5 
        }
        else if (ELEMENT_DIM == 3)
        {
            scale_factor=  6.0; // both the volume of the canonical triangle is 1/6
        }
        return mJacobianDeterminant/scale_factor;
    }
    
    c_vector<double, SPACE_DIM> *pGetWeightedDirection(void)
    {
        if (ELEMENT_DIM >= SPACE_DIM)
        {
            assert(ELEMENT_DIM == SPACE_DIM);
            EXCEPTION("WeightedDirection undefined for fully dimensional element");
            
        }
        return &mWeightedDirection;
    }
    
    /** Get the index of this element
     */
    unsigned GetIndex(void) const
    {
        return mIndex;
    }
    
    /** Update node at the given index
     *  @param rIndex is an local index to which node to change
     *  @param pNode is a pointer to the replacement node
     */
    virtual void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)=0;
    
    
    void ReplaceNode(Node <SPACE_DIM>* pOldNode, Node <SPACE_DIM>* pNewNode)
    {
        for (unsigned i=0; i<mNodes.size(); i++)
        {
            if (mNodes[i]==pOldNode)
            {
                UpdateNode(i,pNewNode);
                return;
            }
        }
        EXCEPTION("You didn't have that node to start with.");
    }
    
    /**
    * Mark an element as having been removed from the mesh.
    * Also notify nodes in the element that it has been removed.
    */
    virtual void MarkAsDeleted()=0;
    
    bool IsDeleted() const
    {
        return mIsDeleted;
    }
    
    void SetIndex(unsigned index)
    {
        mIndex=index;
    }
    
    bool GetOwnership() const
    {
        return mOwnership;
    }
    
    void SetOwnership(bool ownership)
    {
        mOwnership=ownership;
    }
    
    void Flag()
    {
        mFlag = true;
    }
    
    void Unflag()
    {
        mFlag = false;
    }
    
    bool IsFlagged() const
    {
        return mFlag;
    }
    
    
    /***
     * Place in the pIndices array, the global indices (within the stiffness matrix)
     * of the degrees of freedom associated with this element.
     * 
     * @param problemDim the problem dimension e.g. 2 for Bidomain.
     * @param pIndices where to store results: an unsigned array with ELEMENT_DIM+1 entries.
     * 
     */
    void GetStiffnessMatrixGlobalIndices(unsigned problemDim, unsigned* pIndices) const
    {
        
        for (unsigned local_index=0; local_index<ELEMENT_DIM+1; local_index++)
        {
            unsigned node = GetNodeGlobalIndex(local_index);
            
            for (unsigned problem_index=0; problem_index<problemDim; problem_index++)
            {
                //std::cout << local_index*problemDim + problem_index << std::endl;
                pIndices[local_index*problemDim + problem_index] = node*problemDim + problem_index;            
            }   
        } 
    }
};


#endif //_ABSTRACTELEMENT_HPP_
