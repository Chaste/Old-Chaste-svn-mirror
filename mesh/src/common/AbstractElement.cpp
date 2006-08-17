#ifndef _ABSTRACTELEMENT_CPP_
#define _ABSTRACTELEMENT_CPP_
#include "AbstractElement.hpp"

template<int ELEMENT_DIM, int SPACE_DIM>
AbstractElement<ELEMENT_DIM, SPACE_DIM>::AbstractElement(unsigned index,
                                                         std::vector<Node<SPACE_DIM>*> nodes,
                                                         int orderOfBasisFunctions)
        : mIndex(index)
{
    mIsDeleted = false;
    // Sanity checking
    assert(ELEMENT_DIM <= SPACE_DIM);
    //added extra 0.5 to ensure in correct interval for floor() function
    unsigned total_nodes = (unsigned)floor((ELEMENT_DIM+1)*(1 + 0.5*ELEMENT_DIM*(orderOfBasisFunctions - 1)) + 0.5);
    
    assert(nodes.size() == total_nodes);
    
    // Store Node pointers
    mNodes = nodes;
    
    // Specify order of basis functions
    mOrderOfBasisFunctions = orderOfBasisFunctions;
    
    // Create Jacobian
    mpJacobian = NULL;
    mpInverseJacobian = NULL;
    if (ELEMENT_DIM == SPACE_DIM)
    {
        mpJacobian = new c_matrix<double, SPACE_DIM, SPACE_DIM>;
        mpInverseJacobian = new c_matrix<double, SPACE_DIM, SPACE_DIM>;
        
        try
        {
            RefreshJacobian();
        }
        catch (Exception)
        {
            // if the Jacobian is negative the orientation of the element is probably
            // wrong, so swap the last two nodes around.
            
            mNodes[nodes.size()-1] = nodes[nodes.size()-2];
            mNodes[nodes.size()-2] = nodes[nodes.size()-1];
            RefreshJacobian();
        }
        
        // If determinant < 0 then element nodes are listed clockwise.
        // We want them anticlockwise.
        assert(mJacobianDeterminant > 0.0);
        
    }
    else
    {
        RefreshJacobianDeterminant();
    }
    
    
}


template<int ELEMENT_DIM, int SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::RefreshJacobian(void)
{
    for (int i=0; i<ELEMENT_DIM; i++)
    {
        for (int j=0; j<ELEMENT_DIM; j++)
        {
            (*mpJacobian)(i,j) = GetNodeLocation(j+1,i) - GetNodeLocation(0,i);
        }
    }
    
    mJacobianDeterminant = Determinant(*mpJacobian);
    if (mJacobianDeterminant < DBL_EPSILON)
    {
        EXCEPTION("Jacobian determinant is non-positive");
    }
    *mpInverseJacobian   = Inverse(*mpJacobian);
}


template<int ELEMENT_DIM, int SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::RefreshJacobianDeterminant(void)
{
    if (ELEMENT_DIM == SPACE_DIM)
    {
        RefreshJacobian();
        return;
    }
    
    //For boundary elements we only need to know the determinant
    c_vector<double, SPACE_DIM> line_r1_minus_r0;
    c_vector<double, 3> r1_minus_r0;
    c_vector<double, 3> r2_minus_r0;
    switch (ELEMENT_DIM)
    {
        case 0:
            // End point of a line
            mJacobianDeterminant = 1;
            break;
        case 1:
            // Linear edge in a 2D plane or in 3D
            line_r1_minus_r0 = GetNodeLocation(1) - GetNodeLocation(0);
            
            mJacobianDeterminant = norm_2(line_r1_minus_r0); //r1_minus_r0);
            break;
        case 2:
            // Surface triangle in a 3d mesh
            assert(SPACE_DIM == 3);
            r1_minus_r0 = GetNodeLocation(1) - GetNodeLocation(0);
            r2_minus_r0 = GetNodeLocation(2) - GetNodeLocation(0);
            
            mJacobianDeterminant = norm_2( VectorProduct(r1_minus_r0, r2_minus_r0) );
            break;
        default:
            ; // Not going to happen
    }
}

#endif //_ABSTRACTELEMENT_CPP_
