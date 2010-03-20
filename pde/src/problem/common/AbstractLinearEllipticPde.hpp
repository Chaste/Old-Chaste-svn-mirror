/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef _ABSTRACTLINEARELLIPTICPDE_HPP_
#define _ABSTRACTLINEARELLIPTICPDE_HPP_

#include "UblasCustomFunctions.hpp"
#include "ChastePoint.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include <petscvec.h>


/**
 * AbstractLinearEllipticPde class.
 *
 * A general PDE of the form:
 * 0 =   Grad.(DiffusionTerm(x)*Grad(u))
 *     + ComputeConstantInUSourceTerm(x)
 *     + ComputeLinearInUCoeffInSourceTerm(x, u)
 *
 * Parabolic PDEs are be derived from this (AbstractLinearParabolicPde)
 */

//// OLD NOTE: remember this if AbstractPde is brought back
// IMPORTANT NOTE: the inheritance of AbstractPde has to be 'virtual', ie
// "class AbstractCardiacPde : public virtual AbstractPde"
// because AbstractPde will be the top class in a 'dreaded diamond':
//      A
//     / \     A = AbstractPde, B = AbstractCardiac,
//    B   C    C = AbtractLinearElliptic (and AbstractLinearParabolicPde)
//     \ /     D = MonodomainPde
//      D
//
// B and C must use virtual inheritence of A in order for D to only contain 1 instance
// of the member variables in A


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractLinearEllipticPde
{
public:

    /**
     * Compute the constant in u part of the source term, i.e g(x) in
     * Div(D Grad u)  +  f(x)u + g(x) = 0, at a given point.
     *
     * @param rX The point in space
     */
    virtual double ComputeConstantInUSourceTerm(const ChastePoint<SPACE_DIM>& rX)=0;

    /**
     * Compute the coefficient of u in the linear part of the source term, i.e f(x) in
     * Div(D Grad u)  +  f(x)u + g(x) = 0, at a given point in space.
     *
     * @param rX The point in space
     * @param pElement
     */
    virtual double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<SPACE_DIM>& rX,
                                                     Element<ELEMENT_DIM,SPACE_DIM>* pElement)=0;

    /**
     * Compute the diffusion term at a given point.
     *
     * @param rX The point in space at which the diffusion term is computed.
     * @return A matrix.
     */
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX)=0;

    /**
     * Compute the constant in u part of the source term, i.e g(x) in
     * Div(D Grad u)  +  f(x)u + g(x) = 0, at a given node.
     *
     * @param rNode the node
     */
    virtual double ComputeConstantInUSourceTermAtNode(const Node<SPACE_DIM>& rNode);

    /**
     * Compute the coefficient of u in the linear part of the source term, i.e f(x) in
     * Div(D Grad u)  +  f(x)u + g(x) = 0, at a given node.
     *
     * @param rNode the node
     */
    virtual double ComputeLinearInUCoeffInSourceTermAtNode(const Node<SPACE_DIM>& rNode);

    /**
     * Destructor.
     */
    virtual ~AbstractLinearEllipticPde()
    {}
};


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::ComputeConstantInUSourceTermAtNode(const Node<SPACE_DIM>& rNode)
{
    return ComputeConstantInUSourceTerm(rNode.GetPoint());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::ComputeLinearInUCoeffInSourceTermAtNode(const Node<SPACE_DIM>& rNode)
{
    return ComputeLinearInUCoeffInSourceTerm(rNode.GetPoint(), NULL);
}

#endif //_ABSTRACTLINEARELLIPTICPDE_HPP_
