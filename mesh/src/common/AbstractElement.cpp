#ifndef _ABSTRACTELEMENT_CPP_
#define _ABSTRACTELEMENT_CPP_
#include "AbstractElement.hpp"
#include "UblasCustomFunctions.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractElement<ELEMENT_DIM, SPACE_DIM>::AbstractElement(unsigned index,
                                                         const std::vector<Node<SPACE_DIM>*>& rNodes,
                                                         unsigned orderOfBasisFunctions)
        : mIndex(index), mNodes(rNodes)
{
    // Sanity checking
    assert(ELEMENT_DIM <= SPACE_DIM);
    //added extra 0.5 to ensure in correct interval for floor() function
    unsigned total_nodes = (unsigned)floor((ELEMENT_DIM+1)*(1 + 0.5*ELEMENT_DIM*(orderOfBasisFunctions - 1)) + 0.5);
    
    assert(mNodes.size() == total_nodes);
    
    // Specify order of basis functions
    mOrderOfBasisFunctions = orderOfBasisFunctions;
    
    // Initialise flags.
    // This must be done before the Jacobian calculations, or assertions trip.    
    mIsDeleted = false;
    mFlag = false;
    mOwnership = true;

    // This is so we know it's the first time of asking
    mJacobianDeterminant=0.0;
    // Create Jacobian
    try
    {
        RefreshJacobianDeterminant();
    }
    catch (Exception)
    {
        // if the Jacobian is negative the orientation of the element is probably
        // wrong, so swap the last two nodes around.
        
        mNodes[total_nodes-1] = rNodes[total_nodes-2];
        mNodes[total_nodes-2] = rNodes[total_nodes-1];
        RefreshJacobianDeterminant();
    }
    
    // If determinant < 0 then element nodes are listed clockwise.
    // We want them anticlockwise.
    assert(mJacobianDeterminant > 0.0);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::ZeroJacobianDeterminant(void)
{
    mJacobianDeterminant=0.0;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::ZeroWeightedDirection(void)
{
    mWeightedDirection=zero_vector<double>(SPACE_DIM);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::RefreshJacobianDeterminant(bool concreteMove)
{
    if (mIsDeleted)
    {
        EXCEPTION("Attempting to Refresh a deleted element");
    }
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        for (unsigned j=0; j<ELEMENT_DIM; j++)
        {
            mJacobian(i,j) = GetNodeLocation(j+1,i) - GetNodeLocation(0,i);
        }
    }
    if (ELEMENT_DIM == SPACE_DIM)
    {
        mJacobianDeterminant = Determinant(mJacobian);
        if (mJacobianDeterminant <= DBL_EPSILON)
        {
            EXCEPTION("Jacobian determinant is non-positive");
        }
        mInverseJacobian   = Inverse(mJacobian);
        return;
    }
    
    
    bool refresh=false;
    c_vector<double, SPACE_DIM> weighted_direction;
    
    if (mJacobianDeterminant > 0)
    {
        refresh=true;
    }
    //This code is only used when ELEMENT_DIM<SPACE_DIM
    switch (ELEMENT_DIM)
    {
        case 0:
            // End point of a line
            weighted_direction(0)=1.0;
            if (SPACE_DIM == 2)
            {
                weighted_direction(1)=0.0;
            }
            break;
        case 1:
            // Linear edge in a 2D plane or in 3D
            
            weighted_direction=matrix_column<c_matrix<double,SPACE_DIM,SPACE_DIM> >(mJacobian,0);
            break;
        case 2:
            // Surface triangle in a 3d mesh
            assert(SPACE_DIM == 3);
            weighted_direction(0)=-SubDeterminant(mJacobian,0,2);
            weighted_direction(1)= SubDeterminant(mJacobian,1,2);
            weighted_direction(2)=-SubDeterminant(mJacobian,2,2);
            break;
        default:
            ; // Not going to happen
    }
    mJacobianDeterminant = norm_2(weighted_direction);
    if (mJacobianDeterminant < DBL_EPSILON)
    {
        EXCEPTION("Jacobian determinant is zero");
    }
    if (refresh == true)
    {
        if ( inner_prod(mWeightedDirection,weighted_direction) < 0)
        {
            EXCEPTION("Subspace element has changed direction");
        }
    }
    if (concreteMove)
    {
        mWeightedDirection = weighted_direction;
    }
    
}

#endif //_ABSTRACTELEMENT_CPP_
