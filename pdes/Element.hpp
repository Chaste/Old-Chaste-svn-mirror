#ifndef _ELEMENT_HPP_
#define _ELEMENT_HPP_

/**
 * This class defines an Element for use in FEM.
 * 
 * Note that there are two definitions, one for general ELEMENT_DIM and one
 * for ELEMENT_DIM == 0, to avoid infinite loops in the compiler.
 */

#include "Node.hpp"
#include "Point.hpp"
#include "MatrixDouble.hpp"
#include "VectorDouble.hpp"
#include "LinearBasisFunction.cpp"

#include "Exception.hpp"

#include <vector>
#include <iostream>

template <int ELEMENT_DIM, int SPACE_DIM>
class Element
{
private:
    std::vector<Node<SPACE_DIM>*> mNodes;
    const Element<ELEMENT_DIM-1,SPACE_DIM>* mLowerOrderElements[ELEMENT_DIM+1];
    bool mHasLowerOrderElements;

    MatrixDouble *mpJacobian;
    MatrixDouble *mpInverseJacobian;

    double mJacobianDeterminant;

public:
    static const int NUM_CORNER_NODES = ELEMENT_DIM+1;

	/**
	 * Create a new Element.
	 * 
	 * @param nodes A vector of pointers to Node objects, of length at least
	 *     ELEMENT_DIM+1, specifying the nodes associated with this element.
	 *     The corner nodes must come first.
	 * @param createLowerOrderElements Whether to create elements for e.g.
	 *     the faces of a tetrahedral element.
	 *     These only inherit corner nodes.
	 * @param createJacobian Whether to create the Jacobian matrix for mapping
	 *     the element into the appropriate canonical space, e.g. [0,1] in 1D.
	 *     Currently only works for non-sub-elements with straight edges.
	 */
    Element(std::vector<Node<SPACE_DIM>*> nodes,
    	    bool createLowerOrderElements=false, bool createJacobian=true)
    {
    	// Sanity checking
    	assert(ELEMENT_DIM <= SPACE_DIM);
    	
    	// Store Node pointers
    	mNodes = nodes;
    	
    	// Create lower order elements?
    	mHasLowerOrderElements = false;
    	if (createLowerOrderElements)
    	{
    		CreateLowerOrderElements();
    	}
    	
    	// Create Jacobian?
    	mpJacobian = NULL;
    	mpInverseJacobian = NULL;
    	if (createJacobian)
    	{
    		if (ELEMENT_DIM ==  SPACE_DIM)
    		{
    			// TODO: Jacobian for sub-elements    
    			mpJacobian = new MatrixDouble(SPACE_DIM,SPACE_DIM);
    			mpInverseJacobian = new MatrixDouble(SPACE_DIM,SPACE_DIM);
    			
		        for(int i=0; i<SPACE_DIM; i++)
		        {
		            for(int j=0; j<SPACE_DIM; j++)
		            {		                
	                    (*mpJacobian)(i,j) = GetNodeLocation(j+1,i) - GetNodeLocation(0,i);
		            }
		        }
		       
		        *mpInverseJacobian   = mpJacobian->Inverse();
		        mJacobianDeterminant = mpJacobian->Determinant();
    		}
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
			mpJacobian = new MatrixDouble(SPACE_DIM, SPACE_DIM);
			*mpJacobian = *(element.mpJacobian);
		}
		mpInverseJacobian = NULL;
		if (element.mpInverseJacobian != NULL)
		{
			mpInverseJacobian = new MatrixDouble(SPACE_DIM, SPACE_DIM);
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
           	mLowerOrderElements[i] = new Element<ELEMENT_DIM-1, SPACE_DIM>(somenodes,true,false);   
        }
        mHasLowerOrderElements = true;
    }
    
    double GetNodeLocation(int localIndex, int dimension) const
    {
    	assert(dimension < SPACE_DIM);
    	assert(localIndex < mNodes.size());
    	return mNodes[localIndex]->GetPoint()[dimension];
    }
    
    long GetNodeGlobalIndex(int localIndex) const
    {
    	assert(localIndex < mNodes.size());
    	return mNodes[localIndex]->GetIndex();
    }
    
    Node<SPACE_DIM>* GetNode(int localIndex) const
    {
    	assert(localIndex < mNodes.size());
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
	
	const MatrixDouble *GetJacobian(void) const
	{
		return mpJacobian;
	}
	const MatrixDouble *GetInverseJacobian(void) const
	{
		return mpInverseJacobian;
	}
	double GetJacobianDeterminant(void) const
	{
		return mJacobianDeterminant;
	}

};


// This is necessary to stop the compiler recursing into negative dimensions...
template <int SPACE_DIM>
class Element<0, SPACE_DIM>
{
private:
    std::vector<Node<SPACE_DIM>*> mNodes;

    MatrixDouble *mpJacobian;
    MatrixDouble *mpInverseJacobian;

    double mJacobianDeterminant;

public:
    static const int NUM_CORNER_NODES = 1;

	/**
	 * Create a new Element.
	 * 
	 * @param nodes A vector of pointers to Node objects, of length at least
	 *     ELEMENT_DIM+1, specifying the nodes associated with this element.
	 *     The corner nodes must come first.
	 * @param createLowerOrderElements Whether to create elements for e.g.
	 *     the faces of a tetrahedral element.
	 *     These only inherit corner nodes.
	 * @param createJacobian Whether to create the Jacobian matrix for mapping
	 *     the element into the appropriate canonical space, e.g. [0,1] in 1D.
	 *     Currently only works for non-sub-elements with straight edges.
	 */
    Element(std::vector<Node<SPACE_DIM>*> nodes,
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
    		mpJacobian = new MatrixDouble(SPACE_DIM,SPACE_DIM);
    		mpInverseJacobian = new MatrixDouble(SPACE_DIM,SPACE_DIM);
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
     * Copy constructor. This is needed so that copies of an element don't
     * share pointers to the same matrices, which causes problems when copies
     * get destroyed.
     */
    Element(const Element &element)
    {
		mNodes = element.mNodes;
		
		mJacobianDeterminant = element.mJacobianDeterminant;
		mpJacobian = NULL;
		if (element.mpJacobian != NULL)
		{
			mpJacobian = new MatrixDouble(SPACE_DIM, SPACE_DIM);
			*mpJacobian = *(element.mpJacobian);
		}
		mpInverseJacobian = NULL;
		if (element.mpInverseJacobian != NULL)
		{
			mpInverseJacobian = new MatrixDouble(SPACE_DIM, SPACE_DIM);
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
    	assert(localIndex < mNodes.size());
    	return mNodes[localIndex]->GetPoint()[dimension];
    }
    
    long GetNodeGlobalIndex(int localIndex) const
    {
    	assert(localIndex < mNodes.size());
    	return mNodes[localIndex]->GetIndex();
    }
    
    Node<SPACE_DIM>* GetNode(int localIndex) const
    {
    	assert(localIndex < mNodes.size());
    	return mNodes[localIndex];
    }

    int GetNumNodes() const
    {
    	return mNodes.size(); // Will be 1
    }
    
	const MatrixDouble *GetJacobian(void) const
	{
		return mpJacobian;
	}
	const MatrixDouble *GetInverseJacobian(void) const
	{
		return mpInverseJacobian;
	}
	double GetJacobianDeterminant(void) const
	{
		return mJacobianDeterminant;
	}
    
};


#endif //_ELEMENT_HPP_
