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
    void RefreshJacobian(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian);

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
    AbstractTetrahedralElement(unsigned index=INDEX_IS_NOT_USED);

    /**
     * Virtual destructor, since this class has virtual methods.
     */
    virtual ~AbstractTetrahedralElement()
    {}

    /**
     * Get the location of the centroid of the element.
     */
    c_vector<double, SPACE_DIM> CalculateCentroid() const;

    /**
     * Compute the Jacobian for this element.
     *
     * @param rJacobian  the Jacobian matrix
     * @param rJacobianDeterminant  the determinant of the Jacobian
     * @param concreteMove \todo this argument is not used in the method - should it be removed? (defaults to true)
     */
    void CalculateJacobian(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian, double &rJacobianDeterminant, bool concreteMove=true);

    /**
     * Compute the weighted direction for this element.
     *
     * @param rWeightedDirection  the weighted direction vector
     * @param rJacobianDeterminant  the determinant of the Jacobian
     */
    void CalculateWeightedDirection(c_vector<double, SPACE_DIM>& rWeightedDirection, double &rJacobianDeterminant);

    /**
     * Compute the inverse Jacobian for this element.
     *
     * @param rJacobian  the Jacobian matrix
     * @param rJacobianDeterminant  the determinant of the Jacobian
     * @param rInverseJacobian  the inverse Jacobian matrix
     */
    void CalculateInverseJacobian(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian, double &rJacobianDeterminant, c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian); //const

    /** Get the volume of an element (or area in 2d, or length in 1d) 
     * @param determinant a pre-calculated Jacobian determinant for this element.
     * @return volume (which is simply the determinant weighted by the SPACE_DIM)
     */
    double GetVolume(double determinant) const;

    /**
     * Place in the pIndices array, the global indices (within the stiffness matrix)
     * of the degrees of freedom associated with this element.
     *
     * @param problemDim the problem dimension e.g. 2 for Bidomain.
     * @param pIndices where to store results: an unsigned array with ELEMENT_DIM+1 entries.
     *
     */
    void GetStiffnessMatrixGlobalIndices(unsigned problemDim, unsigned* pIndices) const;

};


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::RefreshJacobian(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian)
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::AbstractTetrahedralElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
    : AbstractElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes)
{
    // Sanity checking
    unsigned total_nodes = ELEMENT_DIM+1;
    assert(this->mNodes.size() == total_nodes);

    // This is so we know it's the first time of asking
    // Create Jacobian
    ///\todo We don't want to create new data, calculation and throw the answer away
    c_matrix<double, SPACE_DIM, ELEMENT_DIM> jacobian;

    if (SPACE_DIM == ELEMENT_DIM)
    {
        double det;
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
            // If determinant < 0 then element nodes are listed clockwise.
            // We want them anticlockwise.
            assert(det > 0.0);
        }
    }
    // else - don't bother working out the chirality
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::AbstractTetrahedralElement(unsigned index)
    : AbstractElement<ELEMENT_DIM,SPACE_DIM>(index)
{}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::CalculateJacobian(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian, double &rJacobianDeterminant, bool concreteMove)
{

    assert(ELEMENT_DIM <= SPACE_DIM);
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
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::CalculateWeightedDirection(c_vector<double, SPACE_DIM>& rWeightedDirection, double &rJacobianDeterminant)
{
 
    if(ELEMENT_DIM >= SPACE_DIM)
    {
        assert(ELEMENT_DIM == SPACE_DIM);
        EXCEPTION("WeightedDirection undefined for fully dimensional element");
    }

    c_matrix<double, SPACE_DIM, ELEMENT_DIM> jacobian;
    RefreshJacobian(jacobian);

    //At this point we're only dealing with subspace (ELEMENT_DIM < SPACE_DIM) elem
    //We assume that the rWeightedDirection vector and rJacobianDeterminant (length of vector)
    //are the values from a previous call.

    //This code is only used when ELEMENT_DIM<SPACE_DIM
    switch (ELEMENT_DIM)
    {
        case 0:
            NEVER_REACHED; //See specialised template for ELEMENT_DIM==0
            break;
        case 1:
            // Linear edge in a 2D plane or in 3D

            rWeightedDirection=matrix_column<c_matrix<double,SPACE_DIM,ELEMENT_DIM> >(jacobian,0);
            break;
        case 2:
            // Surface triangle in a 3d mesh
            assert(SPACE_DIM == 3);
            rWeightedDirection(0)=-SubDeterminant(jacobian,0,2);
            rWeightedDirection(1)= SubDeterminant(jacobian,1,2);
            rWeightedDirection(2)=-SubDeterminant(jacobian,2,2);
            break;
        default:
           ; // Not going to happen
    }
    rJacobianDeterminant = norm_2(rWeightedDirection);
    
    if (rJacobianDeterminant < DBL_EPSILON)
    {
        EXCEPTION("Jacobian determinant is zero");
    }

//    if (ELEMENT_DIM == 2)
//    {
//        //Check that weighted-direction and the generalised determinant are giving the same thing
//        double absolute_difference=fabs(rJacobianDeterminant - Determinant(jacobian));
//        double relative_difference=absolute_difference/rJacobianDeterminant;
//        if (relative_difference >= 10 * DBL_EPSILON)
//        {
//#include "Debug.hpp"
//            PRINT_4_VARIABLES(absolute_difference, relative_difference/DBL_EPSILON, rJacobianDeterminant, SPACE_DIM);
//        }
//    }

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::CalculateCentroid() const
{
    c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);
    for (unsigned i=0; i<=ELEMENT_DIM; i++)
    {
        centroid += this->mNodes[i]->rGetLocation();
    }
    return centroid/((double)(ELEMENT_DIM + 1));
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::CalculateInverseJacobian(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian, double &rJacobianDeterminant, c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian)
{
    assert(ELEMENT_DIM <= SPACE_DIM);
    CalculateJacobian(rJacobian, rJacobianDeterminant);

    // CalculateJacobian should make sure that the determinant is not close to zero (or, in fact, negative)
    assert(rJacobianDeterminant > 0.0);
    rInverseJacobian = Inverse(rJacobian);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::GetVolume(double determinant) const
{
    assert(SPACE_DIM == ELEMENT_DIM);

    if (this->mIsDeleted)
    {
        return 0.0;
    }

    double scale_factor = 1.0;

    if (ELEMENT_DIM == 2)
    {
        scale_factor = 2.0;  // both the volume of the canonical triangle is 0.5
    }
    else if (ELEMENT_DIM == 3)
    {
        scale_factor = 6.0; // both the volume of the canonical triangle is 1/6
    }
    return determinant/scale_factor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>::GetStiffnessMatrixGlobalIndices(unsigned problemDim, unsigned* pIndices) const
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


//////////////////////////////////////////////////////////////////////
//                  Specialization for 0d elements                  //
//////////////////////////////////////////////////////////////////////

/**
 * Specialization for 0d elements so we don't get errors from Boost on some
 * compilers.
 */
template<unsigned SPACE_DIM>
class AbstractTetrahedralElement<0, SPACE_DIM> : public AbstractElement<0,SPACE_DIM>
{
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
    AbstractTetrahedralElement(unsigned index=INDEX_IS_NOT_USED);

    /**
     * Virtual destructor, since this class has virtual methods.
     */
    virtual ~AbstractTetrahedralElement()
    {}

    /**
     * Get the location of the centroid of the element.
     */
    c_vector<double, SPACE_DIM> CalculateCentroid() const;

    /**
     * Compute the weighted direction for this element.
     *
     * @param rWeightedDirection  the weighted direction vector
     * @param rJacobianDeterminant  the determinant of the Jacobian
     */
    void CalculateWeightedDirection(c_vector<double, SPACE_DIM>& rWeightedDirection, double &rJacobianDeterminant);

    /**
     * Place in the pIndices array, the global indices (within the stiffness matrix)
     * of the degrees of freedom associated with this element.
     *
     * @param problemDim the problem dimension e.g. 2 for Bidomain.
     * @param pIndices where to store results: an unsigned array with ELEMENT_DIM+1 entries.
     *
     */
    void GetStiffnessMatrixGlobalIndices(unsigned problemDim, unsigned* pIndices) const;
};


template<unsigned SPACE_DIM>
AbstractTetrahedralElement<0, SPACE_DIM>::AbstractTetrahedralElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
    : AbstractElement<0, SPACE_DIM>(index, rNodes)
{
    // Sanity checking
    //unsigned total_nodes = 1;
    assert(this->mNodes.size() == 1);
    assert(SPACE_DIM > 0);

    // This is so we know it's the first time of asking
    // Create Jacobian
    ///\todo We don't want to create new data, calculation and throw the answer away
    c_vector<double, SPACE_DIM> weighted_direction;
    double det;
    
    CalculateWeightedDirection(weighted_direction, det);

    // If determinant < 0 then element nodes are listed clockwise.
    // We want them anticlockwise.
    assert(det > 0.0);
}

template<unsigned SPACE_DIM>
AbstractTetrahedralElement<0, SPACE_DIM>::AbstractTetrahedralElement(unsigned index)
    : AbstractElement<0, SPACE_DIM>(index)
{
}

template<unsigned SPACE_DIM>
void AbstractTetrahedralElement<0, SPACE_DIM>::CalculateWeightedDirection(
        c_vector<double, SPACE_DIM>& rWeightedDirection,
        double &rJacobianDeterminant)
{
    assert(SPACE_DIM > 0);

    // End point of a line
    rWeightedDirection = zero_vector<double>(SPACE_DIM);
    rWeightedDirection(0) = 1.0;

    rJacobianDeterminant = 1.0;
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractTetrahedralElement<0, SPACE_DIM>::CalculateCentroid() const
{
    c_vector<double, SPACE_DIM> centroid = this->mNodes[0]->rGetLocation();
    return centroid;
}

template<unsigned SPACE_DIM>
void AbstractTetrahedralElement<0, SPACE_DIM>::GetStiffnessMatrixGlobalIndices(unsigned problemDim, unsigned* pIndices) const
{
    for (unsigned local_index=0; local_index<1; local_index++)
    {
        unsigned node = this->GetNodeGlobalIndex(local_index);

        for (unsigned problem_index=0; problem_index<problemDim; problem_index++)
        {
            //std::cout << local_index*problemDim + problem_index << std::endl;
            pIndices[local_index*problemDim + problem_index] = node*problemDim + problem_index;
        }
    }
}


#endif //_ABSTRACTTETRAHEDRALELEMENT_HPP_
