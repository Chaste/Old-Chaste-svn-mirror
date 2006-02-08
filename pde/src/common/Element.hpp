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
#include "AbstractMaterial.hpp"

#include "global/src/Exception.hpp"

#include <vector>
#include <iostream>
#include <cmath>

template <int ELEMENT_DIM, int SPACE_DIM>
class Element
{
private:
    std::vector<Node<SPACE_DIM>*> mNodes;
    int mOrderOfBasisFunctions;
    const Element<ELEMENT_DIM-1,SPACE_DIM>* mLowerOrderElements[ELEMENT_DIM+1];
    bool mHasLowerOrderElements;

    MatrixDouble *mpJacobian;
    MatrixDouble *mpInverseJacobian;

    double mJacobianDeterminant;
	
	AbstractMaterial<SPACE_DIM>* mpMaterial;
	
public:
    static const int NUM_CORNER_NODES = ELEMENT_DIM+1;

	/**
	 * Create a new Element.
	 * 
	 * @param nodes A vector of pointers to Node objects, of length at least
	 *     ELEMENT_DIM+1, specifying the nodes associated with this element.
	 *     The corner nodes must come first, and must be listed in anticlockwise
     *     order.
	 * @param createLowerOrderElements Whether to create elements for e.g.
	 *     the faces of a tetrahedral element.
	 *     These only inherit corner nodes.
	 * @param createJacobian Whether to create the Jacobian matrix for mapping
	 *     the element into the appropriate canonical space, e.g. [0,1] in 1D.
	 *     The full matrix is only calculated for non-sub-elements with straight edges.
     *     For surface (boundary) elements we only calculate the determinant,
     *     but this is all that is needed.
	 */
    Element(std::vector<Node<SPACE_DIM>*> nodes,
    		int orderOfBasisFunctions=1,
    	    bool createLowerOrderElements=false, bool createJacobian=true)
    {
    	// Sanity checking
    	assert(ELEMENT_DIM <= SPACE_DIM);
    	//added extra 0.5 to ensure in correct interval for floor() function
    	int total_nodes = (int)floor((ELEMENT_DIM+1)*(1 + 0.5*ELEMENT_DIM*(orderOfBasisFunctions - 1)) + 0.5);

    	assert(nodes.size() == total_nodes);
    	
    	// Store Node pointers
    	mNodes = nodes;
    	
    	// Specify order of basis functions
    	mOrderOfBasisFunctions = orderOfBasisFunctions;
    	
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
    		if (ELEMENT_DIM == SPACE_DIM)
    		{
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
                // If determinant < 0 then element nodes are listed clockwise.
                // We want them anticlockwise.
                assert(mJacobianDeterminant > 0.0);
    		}
            else if (ELEMENT_DIM == SPACE_DIM-1)
            {
                // For boundary elements we only need to know the determinant
                VectorDouble twod_r1_minus_r0(2);
                VectorDouble r1_minus_r0(3);
                VectorDouble r2_minus_r0(3);
                switch (ELEMENT_DIM)
                {
                    case 0:
                        // End point of a line
                        mJacobianDeterminant = 1;
                        break;
                    case 1:
                        // Linear edge in a plane
                        twod_r1_minus_r0(0) = GetNodeLocation(1,0) - GetNodeLocation(0,0); // x1-x0
                        twod_r1_minus_r0(1) = GetNodeLocation(1,1) - GetNodeLocation(0,1); // y1-y0
                        mJacobianDeterminant = twod_r1_minus_r0.L2Norm();
                        break;
                    case 2:
                        // Surface triangle in a 3d mesh
                        r1_minus_r0(0) = GetNodeLocation(1,0) - GetNodeLocation(0,0); // x1-x0
                        r1_minus_r0(1) = GetNodeLocation(1,1) - GetNodeLocation(0,1); // y1-y0
                        r1_minus_r0(2) = GetNodeLocation(1,2) - GetNodeLocation(0,2); // z1-z0
                        r2_minus_r0(0) = GetNodeLocation(2,0) - GetNodeLocation(0,0); // x2-x0
                        r2_minus_r0(1) = GetNodeLocation(2,1) - GetNodeLocation(0,1); // y2-y0
                        r2_minus_r0(2) = GetNodeLocation(2,2) - GetNodeLocation(0,2); // z2-z0
                        mJacobianDeterminant = (r1_minus_r0.VectorProduct(r2_minus_r0)).L2Norm();
                        break;
                    default:
                        assert(0); // TODO? Might want to change this
                }
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
		
        // Allow nodes to keep track of containing elements (but not surface/boundary elements)
        // Only done in copy constructor, since that is what is called to put elements
        // in the vector contained in ConformingTetrahedralMesh.
        if (ELEMENT_DIM == SPACE_DIM)
        {
            for (int i=0; i<mNodes.size(); i++)
            {
                mNodes[i]->AddElement((const void*)this);
            }
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
           	mLowerOrderElements[i] = new Element<ELEMENT_DIM-1, SPACE_DIM>(somenodes,1,true,false);   
        }
        mHasLowerOrderElements = true;
    }
    
    double GetNodeLocation(int localIndex, int dimension) const
    {
    	assert(dimension < SPACE_DIM);
    	assert(localIndex < mNodes.size());
    	return mNodes[localIndex]->rGetPoint()[dimension];
    }
    
    long GetNodeGlobalIndex(int localIndex) const
    {
    	assert(localIndex < mNodes.size());
    	return mNodes[localIndex]->GetIndex();
    }
    
    const Node<SPACE_DIM>* GetNode(int localIndex) const
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
     * Create a new element from a Node
     * The element has ELEMENT_DIM=0 and
     * SPACE_DIM identical to that of the node from which it is constructed
     * 
     */
    Element(Node<SPACE_DIM> *node,
            bool createLowerOrderElements=false, bool createJacobian=true)
    {
        // Store Node pointer
        mNodes.push_back(node);

        // Create Jacobian?
        mpJacobian = NULL;
        mpInverseJacobian = NULL;
        if (createJacobian)
        {
            mpJacobian = new MatrixDouble(1,1);
            mpInverseJacobian = new MatrixDouble(1,1);
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
		
        // Allow nodes to keep track of containing elements (but not surface/boundary elements)
        // Only done in copy constructor, since that is what is called to put elements
        // in the vector contained in ConformingTetrahedralMesh.
        if (0 == SPACE_DIM)
        {
            for (int i=0; i<mNodes.size(); i++)
            {
                mNodes[i]->AddElement((const void*)this);
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
    
    const Node<SPACE_DIM>* GetNode(int localIndex) const
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
