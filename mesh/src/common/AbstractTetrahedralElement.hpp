/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef _ABSTRACTTETRAHEDRALELEMENT_HPP_
#define _ABSTRACTTETRAHEDRALELEMENT_HPP_

#include "AbstractElement1.hpp"

/**
 * This class defines an Element for use in the Finite Element Method.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractTetrahedralElement : public AbstractElement<ELEMENT_DIM,SPACE_DIM>
{
protected:
    c_matrix<double, SPACE_DIM, SPACE_DIM> mJacobian;
    c_matrix<double, SPACE_DIM, SPACE_DIM> mInverseJacobian;
    
    /*< Holds an area-weighted normal or direction.  Only used when ELEMENT_DIM < SPACE_DIM */
    c_vector<double, SPACE_DIM> mWeightedDirection; 
    double mJacobianDeterminant;

    /**
     * Method that constructs the element. This is required because of having
     * two different copy constructors, one with a new index and another without
     */
    void CommonConstructor(const AbstractTetrahedralElement& rElement)
    {
        mJacobianDeterminant = rElement.mJacobianDeterminant;
        mJacobian = rElement.mJacobian;
        mInverseJacobian = rElement.mInverseJacobian;
        mWeightedDirection = rElement.mWeightedDirection;
        
        AbstractElement<ELEMENT_DIM,SPACE_DIM>::CommonConstructor(rElement);
    }


public:
    ///Main constructor
    AbstractTetrahedralElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Copy constructor. This is needed so that copies of an element don't
     * share pointers to the same matrices, which causes problems when copies
     * get destroyed.
     */
    AbstractTetrahedralElement(const AbstractTetrahedralElement &element)
        : AbstractElement<ELEMENT_DIM, SPACE_DIM>(element)
    {
//        CommonConstructor(element);
    }

    /**
     * \todo Why does the default constructor not do anything?
     */
    AbstractTetrahedralElement()
    {}

    /**
     * Element assignment - make this element equal to the other one.
     */
    virtual AbstractTetrahedralElement& operator=(const AbstractTetrahedralElement &element)
    {
        CommonConstructor(element);
        return *this;
    }

    virtual ~AbstractTetrahedralElement()
    {}

    void RefreshJacobianDeterminant(bool concreteMove=true);
    void ZeroJacobianDeterminant(void);
    void ZeroWeightedDirection(void);


    c_vector<double, SPACE_DIM> CalculateCentroid() const
    {
        c_vector<double, SPACE_DIM> centroid=zero_vector<double>(SPACE_DIM);
        for (unsigned i=0; i<=ELEMENT_DIM; i++)
        {
            centroid += this->mNodes[i]->rGetLocation();
        }
        return centroid/((double)(ELEMENT_DIM + 1));
    }

///////////////////////////////////
    const c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> *CalculateJacobian(void) const
    {
        return &mJacobian;
    }
    const c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> *CalculateInverseJacobian(void) const
    {
        return &mInverseJacobian;
    }
    double CalculateJacobianDeterminant(void) const
    {
        return mJacobianDeterminant;
    }
///////////////////////////////////

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
            scale_factor= 6.0; // both the volume of the canonical triangle is 1/6
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



    /** Update node at the given index
     *  @param rIndex is an local index to which node to change
     *  @param pNode is a pointer to the replacement node
     */
    virtual void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)=0;


    void ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode)
    {
        for (unsigned i=0; i<this->mNodes.size(); i++)
        {
            if (this->mNodes[i]==pOldNode)
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

    /**
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
            unsigned node = this->GetNodeGlobalIndex(local_index);

            for (unsigned problem_index=0; problem_index<problemDim; problem_index++)
            {
                //std::cout << local_index*problemDim + problem_index << std::endl;
                pIndices[local_index*problemDim + problem_index] = node*problemDim + problem_index;
            }
        }
    }
};



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::AbstractTetrahedralElement(unsigned index,
                                                                               const std::vector<Node<SPACE_DIM>*>& rNodes)
        : AbstractElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes)
{
    // Sanity checking
    unsigned total_nodes = ELEMENT_DIM+1;
    assert(this->mNodes.size() == total_nodes);

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

        this->mNodes[total_nodes-1] = rNodes[total_nodes-2];
        this->mNodes[total_nodes-2] = rNodes[total_nodes-1];
        RefreshJacobianDeterminant();
    }

    // If determinant < 0 then element nodes are listed clockwise.
    // We want them anticlockwise.
    assert(mJacobianDeterminant > 0.0);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::ZeroJacobianDeterminant(void)
{
    mJacobianDeterminant=0.0;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::ZeroWeightedDirection(void)
{
    mWeightedDirection=zero_vector<double>(SPACE_DIM);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::RefreshJacobianDeterminant(bool concreteMove)
{
    if (this->mIsDeleted)
    {
        EXCEPTION("Attempting to Refresh a deleted element");
    }
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        for (unsigned j=0; j!=ELEMENT_DIM; j++) //Does a j<ELEMENT_DIM without ever having to test j<0U (#186: pointless comparison of unsigned integer with zero)
        {
            mJacobian(i,j) = this->GetNodeLocation(j+1,i) - this->GetNodeLocation(0,i);
        }
    }


    if (ELEMENT_DIM == SPACE_DIM)
    {
        mJacobianDeterminant = Determinant(mJacobian);
        if (mJacobianDeterminant <= DBL_EPSILON)
        {
            std::stringstream string_stream;
            string_stream << "Jacobian determinant is non-positive: "
                          << "determinant = " << mJacobianDeterminant
                          << " for element " << this->mIndex;
            EXCEPTION(string_stream.str());
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
    double jacobian_determinant = norm_2(weighted_direction);
    if (jacobian_determinant < DBL_EPSILON)
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
        assert(ELEMENT_DIM < SPACE_DIM);
        mJacobianDeterminant = jacobian_determinant;
        mWeightedDirection = weighted_direction;
    }
}

#endif //_ABSTRACTTETRAHEDRALELEMENT_HPP_
