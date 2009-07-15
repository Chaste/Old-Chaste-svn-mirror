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


#ifndef _FOURTHORDERTENSOR_HPP_
#define _FOURTHORDERTENSOR_HPP_

#include <cassert>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;
#include "Exception.hpp"


/**
 *  FourthOrderTensor
 *
 *  A class of fourth order tensors (i.e. tensors with four indices), over arbitrary dimension.
 *
 */
template<unsigned DIM>
class FourthOrderTensor
{
private:

    std::vector<double> mData;  /**< The components of the tensor. */
    unsigned mDimSqd;           /**< The squared dimension, DIM^2. */
    unsigned mDimCubed;         /**< The cubed dimension, DIM^3. */
    unsigned mDimToFour;        /**< The fourth power of the dimension, DIM^4. */

public:

    /**
     * Constructor.
     */
    FourthOrderTensor();

    /**
     *  Set to be the inner product of another fourth order tensor and a matrix
     *
     *  @param rTensor A fourth order tensor
     *  @param rMatrix A matrix
     *  @param component  The component in the fourth order tensor with which to sum
     *    (indexed from ZERO)
     *
     *  ie. if component=0, X_{RM} T_{MNPQ} is returned
     *  ie. if component=2, X_{RQ} T_{MNPQ} is returned
     *
     */
    void SetAsProduct(FourthOrderTensor<DIM>& rTensor, const c_matrix<double,DIM,DIM>& rMatrix, unsigned component);

    /**
     * Access the MNPQ-component of the tensor.
     *
     * @param M  first index
     * @param N  second index
     * @param P  third index
     * @param Q  fourth index
     */
    double& operator()(unsigned M, unsigned N, unsigned P, unsigned Q);

    /**
     * Set all components of the tensor to zero.
     */
    void Zero();

};

#endif //_FOURTHORDERTENSOR_HPP_
