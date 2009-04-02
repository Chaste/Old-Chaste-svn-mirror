/*

Copyright (C) University of Oxford, 2005-2009

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

#include "AbstractElement.hpp"

/**
 * This abstract class defines a tetrahedral element for use in the Finite Element Method.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractTetrahedralElement : public AbstractElement<ELEMENT_DIM,SPACE_DIM>
{
protected:

    /**
     * Refresh the Jacobian for this element.
     *
     * @param rJacobian  the Jacobian matrix
     */
    void RefreshJacobian(c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian)
    {
        if (this->mIsDeleted)
        {
            EXCEPTION("Attempting to Refresh a deleted element");
        }
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            for (unsigned j=0; j!=ELEMENT_DIM; j++) //Does a j<ELEMENT_DIM without ever having to test j<0U (#186: pointless comparison of unsigned integer with zero)
            {
                rJacobian(i,j) = this->GetNodeLocation(j+1,i) - this->GetNodeLocation(0,i);
            }
        }
    }

public:

    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param index  the index of the element in the mesh
     * @param rNodes  the nodes owned by the element
     */
    AbstractTetrahedralElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Default constructor, which doesn't fill in any nodes.
     * The nodes must be added later.
     *
     * @param index  the index of the element in the mesh (defaults to INDEX_IS_NOT_USED)
     */
    AbstractTetrahedralElement(unsigned index=INDEX_IS_NOT_USED)
        : AbstractElement<ELEMENT_DIM,SPACE_DIM>(index)
    {}

    /**
     * Virtual destructor, since this class has virtual methods.
     */
    virtual ~AbstractTetrahedralElement()
    {}

    /**
     * \todo This method does not appear to be used anywhere - remove it?
     */
    void ZeroJacobianDeterminant(void);

    /**
     * \todo This method does not appear to be used anywhere - remove it?
     */
    void ZeroWeightedDirection(void);

    /**
     * Get the location of the centroid of the element.
     */
    c_vector<double, SPACE_DIM> CalculateCentroid() const
    {
        c_vector<double, SPACE_DIM> centroid=zero_vector<double>(SPACE_DIM);
        for (unsigned i=0; i<=ELEMENT_DIM; i++)
        {
            centroid += this->mNodes[i]->rGetLocation();
        }
        return centroid/((double)(ELEMENT_DIM + 1));
    }

    /**
     * Compute the Jacobian for this element.
     *
     * @param rJacobian  the Jacobian matrix
     * @param rJacobianDeterminant  the determinant of the Jacobian
     * @param concreteMove \todo this argument is not used in the method - should it be removed? (defaults to true)
     */
    void CalculateJacobian(c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian, double &rJacobianDeterminant, bool concreteMove=true);

    /**
     * Compute the weighted direction for this element.
     *
     * @param rWeightedDirection  the weighted direction vector
     * @param rJacobianDeterminant  the determinant of the Jacobian
     * @param concreteMove \todo this argument is not used in the method - should it be removed? (defaults to true)
     */
    void CalculateWeightedDirection(c_vector<double, SPACE_DIM>& rWeightedDirection, double &rJacobianDeterminant, bool concreteMove=true);

    /**
     * Compute the inverse Jacobian for this element.
     *
     * @param rJacobian  the Jacobian matrix
     * @param rJacobianDeterminant  the determinant of the Jacobian
     * @param rInverseJacobian  the inverse Jacobian matrix
     */
    void CalculateInverseJacobian(c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian, double &rJacobianDeterminant, c_matrix<double, SPACE_DIM, SPACE_DIM>& rInverseJacobian) //const
    {
        assert(ELEMENT_DIM==SPACE_DIM);
        CalculateJacobian(rJacobian, rJacobianDeterminant);
        //CalculateJacobian should make sure that the determinant is not close to zero (or, in fact, negative)
        assert(rJacobianDeterminant > 0.0);
        rInverseJacobian = Inverse(rJacobian);
    }


///\todo Re-implement
    /** Get the volume of an element (or area in 2d, or length in 1d) */
    double GetVolume(void) //const?
    {
        assert(SPACE_DIM == ELEMENT_DIM);

        if (this->mIsDeleted)
        {
            return 0.0;
        }

        // Create Jacobian
        ///\todo We don't want to create new data, calculation and throw the answer away
        c_matrix<double, SPACE_DIM, SPACE_DIM> jacobian;
        double determinant;

        CalculateJacobian(jacobian, determinant);
        double scale_factor = 1.0;

        if (ELEMENT_DIM == 2)
        {
            scale_factor = 2.0;  // both the volume of the canonical triangle is 0.5
        }
        else if (ELEMENT_DIM == 3)
        {
            scale_factor= 6.0; // both the volume of the canonical triangle is 1/6
        }
        return determinant/scale_factor;
    }

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
    // Create Jacobian
    ///\todo We don't want to create new data, calculation and throw the answer away
    c_matrix<double, SPACE_DIM, SPACE_DIM> jacobian;
    c_vector<double, SPACE_DIM> weighted_direction;
    double det;

    if (SPACE_DIM == ELEMENT_DIM)
    {
        try
        {
            CalculateJacobian(jacobian, det);
        }
        catch (Exception)
        {
            // if the Jacobian is negative the orientation of the element is probably
            // wrong, so swap the last two nodes around.

            this->mNodes[total_nodes-1] = rNodes[total_nodes-2];
            this->mNodes[total_nodes-2] = rNodes[total_nodes-1];

            CalculateJacobian(jacobian, det);
        }
    }
    else
    {
        CalculateWeightedDirection(weighted_direction, det);
    }

    // If determinant < 0 then element nodes are listed clockwise.
    // We want them anticlockwise.
    assert(det > 0.0);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::CalculateJacobian(c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian, double &rJacobianDeterminant, bool concreteMove)
{

    assert(ELEMENT_DIM == SPACE_DIM);
    RefreshJacobian(rJacobian);

    {
        rJacobianDeterminant = Determinant(rJacobian);
        if (rJacobianDeterminant <= DBL_EPSILON)
        {
            std::stringstream string_stream;
            string_stream << "Jacobian determinant is non-positive: "
                          << "determinant = " << rJacobianDeterminant
                          << " for element " << this->mIndex;
            EXCEPTION(string_stream.str());
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::CalculateWeightedDirection(c_vector<double, SPACE_DIM>& rWeightedDirection, double &rJacobianDeterminant, bool concreteMove)
{
    if(ELEMENT_DIM >= SPACE_DIM)
    {
        assert(ELEMENT_DIM == SPACE_DIM);
        EXCEPTION("WeightedDirection undefined for fully dimensional element");
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> jacobian;
    RefreshJacobian(jacobian);

    //At this point we're only dealing with subspace (ELEMENT_DIM < SPACE_DIM) elem
    //We assume that the rWeightedDirection vector and rJacobianDeterminant (length of vector)
    //are the values from a previous call.
    //rJacobianDeterminant=0.0 signifies that this is the first calculation on this element.
    c_vector<double, SPACE_DIM> weighted_direction;
//    bool refresh=false;
//
//    if (rJacobianDeterminant > 0) // 767 Checking against the reference we are getting?
//    {
//        refresh=true;
//    }


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

            weighted_direction=matrix_column<c_matrix<double,SPACE_DIM,SPACE_DIM> >(jacobian,0);
            break;
        case 2:
            // Surface triangle in a 3d mesh
            assert(SPACE_DIM == 3);
            weighted_direction(0)=-SubDeterminant(jacobian,0,2);
            weighted_direction(1)= SubDeterminant(jacobian,1,2);
            weighted_direction(2)=-SubDeterminant(jacobian,2,2);
            break;
        default:
           ; // Not going to happen
    }
    double jacobian_determinant = norm_2(weighted_direction);
    if (jacobian_determinant < DBL_EPSILON)
    {
        EXCEPTION("Jacobian determinant is zero");
    }
//    if (refresh == true)
//    {
//        if ( inner_prod(rWeightedDirection, weighted_direction) < 0)
//        {
//            EXCEPTION("Subspace element has changed direction");
//        }
//    }
//    if (concreteMove)
//    {
        rJacobianDeterminant = jacobian_determinant;
        rWeightedDirection = weighted_direction;
//    }
}

#endif //_ABSTRACTTETRAHEDRALELEMENT_HPP_
