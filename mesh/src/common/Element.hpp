#ifndef _ELEMENT_HPP_
#define _ELEMENT_HPP_
#include <float.h>
/**
 * This class defines an Element for use in FEM.
 * 
 * Note that there are two definitions, one for general ELEMENT_DIM and one
 * for ELEMENT_DIM == 0, to avoid infinite loops in the compiler.
 */

#include "Node.hpp"
#include "Point.hpp"
#include "UblasCustomFunctions.hpp"
#include "AbstractMaterial.hpp"

#include "Exception.hpp"

#include <vector>
#include <cmath>

// When creating an element within a mesh one needs to specify its global index
// If the element is not used within a mesh the following
// constant is used instead.
const unsigned INDEX_IS_NOT_USED=0;


template <int ELEMENT_DIM, int SPACE_DIM>
class Element
{
private:
    unsigned mIndex;
    std::vector<Node<SPACE_DIM>*> mNodes;
    int mOrderOfBasisFunctions;
    const Element<ELEMENT_DIM-1,SPACE_DIM>* mLowerOrderElements[ELEMENT_DIM+1];
    bool mHasLowerOrderElements;

    c_matrix<double, SPACE_DIM, SPACE_DIM> *mpJacobian;
    c_matrix<double, SPACE_DIM, SPACE_DIM> *mpInverseJacobian;

    double mJacobianDeterminant;
	
	AbstractMaterial<SPACE_DIM>* mpMaterial;
    void RefreshJacobian(void);
	
    
public:
    static const int NUM_CORNER_NODES = ELEMENT_DIM+1;

	Element(unsigned index,
            std::vector<Node<SPACE_DIM>*> nodes, 
            int orderOfBasisFunctions=1,
            bool createLowerOrderElements=false, 
            bool createJacobian=true);
      
         
            
    
    /**
     * Copy constructor. This is needed so that copies of an element don't
     * share pointers to the same matrices, which causes problems when copies
     * get destroyed.
     */
    Element(const Element &element)
    {
		mNodes = element.mNodes;
		mIndex = element.mIndex;
        // Allow nodes to keep track of containing elements (but not surface/boundary elements)
        // Only done in copy constructor, since that is what is called to put elements
        // in the vector contained in ConformingTetrahedralMesh.
        
        for (unsigned i=0; i<mNodes.size(); i++)
        {
            mNodes[i]->AddElement(mIndex);
        }

		mHasLowerOrderElements = element.mHasLowerOrderElements;
    	if (mHasLowerOrderElements)
    	{
    		// Copy lower order elements, rather than pointers
    		for (int i=0; i<ELEMENT_DIM+1; i++)
    		{
    			mLowerOrderElements[i] =
    				new Element<ELEMENT_DIM-1, SPACE_DIM>(*(element.mLowerOrderElements[i]));
    		}
    	}

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
     * Free memory potentially allocated in the constructor (or elsewhere) for
     * holding lower order elements and the Jacobian and its inverse.
     */
    ~Element()
    {
    	if (mHasLowerOrderElements)
    	{
    		for (int i=0; i<ELEMENT_DIM+1; i++)
    		{
    			delete mLowerOrderElements[i];
    		}
    	}
     	
    	if (mpJacobian != NULL)
    	{
    		delete mpJacobian;
    	}
    	if (mpInverseJacobian != NULL)
    	{
    		delete mpInverseJacobian;
    	}
    }
    
    void RefreshJacobianDeterminant(void);
    
    void AddInternalNode(const Node<SPACE_DIM>* internalNodeToAdd)
    {
    	assert(mOrderOfBasisFunctions > 1);
    	assert(mNodes.size() - NUM_CORNER_NODES < 0.5*SPACE_DIM*(SPACE_DIM+1));
    	
    	mNodes.push_back(internalNodeToAdd);
    }
    
    void CreateLowerOrderElements()
    {
    	for (int i = 0; i < ELEMENT_DIM + 1; i++)
        {
	        std::vector<Node<SPACE_DIM>*> somenodes;
			for (int j = 0; j < ELEMENT_DIM; j++)
           	{
               	int nodeIndex = ((i + j + 1) % (ELEMENT_DIM + 1));
               	somenodes.push_back(mNodes[nodeIndex]);
           	}
           	//WARNING: Lower Order Elements are not constructed with internal nodes at present
           	mLowerOrderElements[i] = new Element<ELEMENT_DIM-1, SPACE_DIM>(0,somenodes,1,true,false);   
        }
        mHasLowerOrderElements = true;
    }
    
    double GetNodeLocation(int localIndex, int dimension) const
    {
    	assert(dimension < SPACE_DIM);
    	assert((unsigned)localIndex < mNodes.size());
    	return mNodes[localIndex]->rGetPoint()[dimension];
    }
    
    long GetNodeGlobalIndex(int localIndex) const
    {
    	assert((unsigned)localIndex < mNodes.size());
    	return mNodes[localIndex]->GetIndex();
    }
    
    const Node<SPACE_DIM>* GetNode(int localIndex) const
    {
    	assert((unsigned)localIndex < mNodes.size());
    	return mNodes[localIndex];
    }
    
    int GetNumNodes() const
    {
    	return mNodes.size();
    }
    
    const Element<ELEMENT_DIM-1,SPACE_DIM>* GetLowerOrderElement(int localIndex) const
    {
    	assert(mHasLowerOrderElements);
    	assert(localIndex < ELEMENT_DIM+1);
		return mLowerOrderElements[localIndex];   
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
	
	/** Set the material details for this element
	 */
	void SetMaterial(AbstractMaterial<SPACE_DIM>* pMaterial)
	{
		mpMaterial = pMaterial;
	}

	/** Get the material details for this element
	 */
	AbstractMaterial<SPACE_DIM>* GetMaterial() const
	{

		if(!mpMaterial)
		{
			assert(0);
		}
		
		return mpMaterial;
	}
};


// This is necessary to stop the compiler recursing into negative dimensions...
template <int SPACE_DIM>
class Element<0, SPACE_DIM>
{
private:
    unsigned mIndex;
    std::vector<Node<SPACE_DIM>*> mNodes;

    c_matrix<double, SPACE_DIM, SPACE_DIM> *mpJacobian;
    c_matrix<double, SPACE_DIM, SPACE_DIM> *mpInverseJacobian;

    double mJacobianDeterminant;

public:
    static const int NUM_CORNER_NODES = 1;

	/**
	 * Create a new Element.
	 * 
	 * @param nodes A vector of pointers to Node objects, of length at least
	 *     ELEMENT_DIM+1, specifying the nodes associated with this element.
	 *     The corner nodes must come first.
     * @param index - Global index to this element (not used)
	 * @param createLowerOrderElements Whether to create elements for e.g.
	 *     the faces of a tetrahedral element.
	 *     These only inherit corner nodes.
	 * @param createJacobian Whether to create the Jacobian matrix for mapping
	 *     the element into the appropriate canonical space, e.g. [0,1] in 1D.
	 *     Currently only works for non-sub-elements with straight edges.
	 */
    Element(unsigned index,
    		std::vector<Node<SPACE_DIM>*> nodes,
            int notUsed=1,
    	    bool createLowerOrderElements=false, bool createJacobian=true)
    {
    	// Store Node pointers
    	assert(nodes.size() == 1);
    	mNodes = nodes;

     	// Create Jacobian?
    	mpJacobian = NULL;
    	mpInverseJacobian = NULL;
    	if (createJacobian)
    	{
    		mpJacobian = new c_matrix<double, SPACE_DIM, SPACE_DIM>;
    		mpInverseJacobian = new c_matrix<double, SPACE_DIM, SPACE_DIM>;
    		for (int i=0; i<SPACE_DIM; i++)
    		{
    			for (int j=0; j<SPACE_DIM; j++)
    			{
					(*mpJacobian)(i,j) = 0.0;
					(*mpInverseJacobian)(i,j) = 0.0;
    			}
    			(*mpJacobian)(i,i) = 1.0;
				(*mpInverseJacobian)(i,i) = 1.0;
    		}
			mJacobianDeterminant = 1.0;
    	}
    }
    
    /**
     * Create a new element from a Node
     * The element has ELEMENT_DIM=0 and
     * SPACE_DIM identical to that of the node from which it is constructed
     * 
     */
    Element(unsigned index,
            Node<SPACE_DIM> *node, 
            bool createLowerOrderElements=false, 
            bool createJacobian=true)
            : mIndex(index)
    {
        // Store Node pointer
        mNodes.push_back(node);
        
        // Create Jacobian?
        mpJacobian = NULL;
        mpInverseJacobian = NULL;
        if (createJacobian)
        {
            // Note: this might be a source of errors for simulations
            // with space dim not equal to elem dim
            mpJacobian = new c_matrix<double, SPACE_DIM, SPACE_DIM>;
            mpInverseJacobian = new c_matrix<double, SPACE_DIM, SPACE_DIM>;
            (*mpJacobian)(0,0) = 1.0;
            (*mpInverseJacobian)(0,0) = 1.0;
            mJacobianDeterminant = 1.0;
        }
    }
   
   /**
     * Copy constructor. This is needed so that copies of an element don't
     * share pointers to the same matrices, which causes problems when copies
     * get destroyed.
     */
    Element(const Element &element)
    {
		mNodes = element.mNodes;
		mIndex = element.mIndex;
        
        // Allow nodes to keep track of containing elements (but not surface/boundary elements)
        // Only done in copy constructor, since that is what is called to put elements
        // in the vector contained in ConformingTetrahedralMesh.
        if (0 == SPACE_DIM)
        {
            for (unsigned i=0; i<mNodes.size(); i++)
            {
                mNodes[i]->AddElement(mIndex);
            }
        }
 
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
     * Free memory potentially allocated in the constructor for
     * holding the Jacobian and its inverse.
     */
    ~Element()
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

    double GetNodeLocation(int localIndex, int dimension) const
    {
    	assert(dimension < SPACE_DIM);
    	assert((unsigned)localIndex < mNodes.size());
    	return mNodes[localIndex]->GetPoint()[dimension];
    }
    
    long GetNodeGlobalIndex(int localIndex) const
    {
    	assert((unsigned)localIndex < mNodes.size());
    	return mNodes[localIndex]->GetIndex();
    }
    
    const Node<SPACE_DIM>* GetNode(int localIndex) const
    {
    	assert((unsigned)localIndex < mNodes.size());
    	return mNodes[localIndex];
    }

    int GetNumNodes() const
    {
    	return mNodes.size(); // Will be 1
    }
    
	const c_matrix<double, SPACE_DIM, SPACE_DIM> *GetJacobian(void) const
	{
		return mpJacobian;
	}
	const c_matrix<double, SPACE_DIM, SPACE_DIM> *GetInverseJacobian(void) const
	{
		return mpInverseJacobian;
	}
	double GetJacobianDeterminant(void) const
	{
		return mJacobianDeterminant;
	}
    
};


#endif //_ELEMENT_HPP_
