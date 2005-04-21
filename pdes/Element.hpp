#ifndef _ELEMENT_HPP_
#define _ELEMENT_HPP_

/**
 * This class defines an Element for use in FEM.
 */

#include "Node.hpp"
#include "Point.hpp"
#include "MatrixDouble.hpp"
#include "VectorDouble.hpp"
#include "LinearBasisFunction.cpp"

#include <vector>

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
	 *     Also only implemented in 1 element and space dimension.
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
    	if (createJacobian)
    	{
    		if (ELEMENT_DIM == 1 && SPACE_DIM == 1)
    		{
    			// TODO: Jacobian for sub-elements    
    			mpJacobian = new MatrixDouble(1,1);
    			mpInverseJacobian = new MatrixDouble(1,1);
    			
		        for(int i=0; i<1; i++)
		        {
		            for(int j=0; j<1; j++)
		            {    
		                //loop over corner nodes
		                for(int node=0; node<1+1; node++)
		                {
		                    (*mpJacobian)(i,j) = GetNodeLocation(1,0) - GetNodeLocation(0,0);
		                }
		            }            
		        }
		        *mpInverseJacobian   = mpJacobian->Inverse();
		        mJacobianDeterminant = mpJacobian->Determinant();
    		}
	    }
    }
    
    /**
     * Free memory potentially allocated in the constructor (or elsewhere) for
     * holding lower order elements.
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
    	if (createJacobian)
    	{
    		mpJacobian = new MatrixDouble(1,1);
    		mpInverseJacobian = new MatrixDouble(1,1);
			(*mpJacobian)(0,0) = 1.0;
			(*mpInverseJacobian)(0,0) = 1.0;
			mJacobianDeterminant = 1.0;
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
