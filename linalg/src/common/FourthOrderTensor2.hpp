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


#ifndef _FOURTHORDERTENSOR2_HPP_
#define _FOURTHORDERTENSOR2_HPP_

#include <cassert>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;
#include "Exception.hpp"


/**
 *  FourthOrderTensor
 *
 *  A class of fourth order tensors (ie tensors with four indices), over arbitrary dimension
 *
 */
template<unsigned DIM>
class FourthOrderTensor2
{
private:

    std::vector<double> mData;
    unsigned mDimSqd;
    unsigned mDimCubed;
    unsigned mDimToFour;

public:

    /**
     * Constructor.
     */
    FourthOrderTensor2();

    /**
     *  Set to be the inner product of another fourth order tensor and a matrix
     *
     *  @param tensor A fourth order tensor
     *  @param matrix A Deal.II matrix
     *  @param component  The component in the fourth order tensor with which to sum
     *    (indexed from ZERO)
     *
     *  ie. if component=0, X_{RM} T_{MNPQ} is returned
     *  ie. if component=2, X_{RQ} T_{MNPQ} is returned
     *
     */
    void SetAsProduct(FourthOrderTensor2<DIM>& tensor, const c_matrix<double,DIM,DIM>& matrix, unsigned component);

    double& operator()(unsigned M, unsigned N, unsigned P, unsigned Q);

    void Zero();

};

template<unsigned DIM>
FourthOrderTensor2<DIM>::FourthOrderTensor2()
{
    // check dim>0 but <4
    assert(DIM>0);
    assert(DIM<4);

    mDimSqd = DIM*DIM;
    mDimCubed = DIM*DIM*DIM;
    mDimToFour = DIM*DIM*DIM*DIM;

    // allocate memory and zero entries
    mData.resize(mDimToFour, 0.0);
}

template<unsigned DIM>
void FourthOrderTensor2<DIM>::SetAsProduct(FourthOrderTensor2<DIM>& tensor, const c_matrix<double,DIM,DIM>& matrix, unsigned component)
{
    Zero();

    // messy repeated code but not sure how to do this neatly and efficiently..
    switch (component)
    {
        case 0:
        {
            for(unsigned M=0; M<DIM; M++)
            {
                for(unsigned N=0; N<DIM; N++)
                {
                    for(unsigned P=0; P<DIM; P++)
                    {
                        for(unsigned Q=0; Q<DIM; Q++)
                        {
                            unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
                            for(unsigned s=0; s<DIM; s++)
                            {
                                mData[index] += matrix(M,s) * tensor(s,N,P,Q);
                            }
                        }
                    }
                }
            }
            break;
        }
        case 1:
        {
            for(unsigned M=0; M<DIM; M++)
            {
                for(unsigned N=0; N<DIM; N++)
                {
                    for(unsigned P=0; P<DIM; P++)
                    {
                        for(unsigned Q=0; Q<DIM; Q++)
                        {
                            unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
                            for(unsigned s=0; s<DIM; s++)
                            {
                                mData[index] += matrix(N,s) * tensor(M,s,P,Q);
                            }
                        }
                    }
                }
            }
            break;
        }
        case 2:
        {
            for(unsigned M=0; M<DIM; M++)
            {
                for(unsigned N=0; N<DIM; N++)
                {
                    for(unsigned P=0; P<DIM; P++)
                    {
                        for(unsigned Q=0; Q<DIM; Q++)
                        {
                            unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
                            for(unsigned s=0; s<DIM; s++)
                            {
                                mData[index] += matrix(P,s) * tensor(M,N,s,Q);
                            }
                        }
                    }
                }
            }
            break;
        }
        case 3:
        {
            for(unsigned M=0; M<DIM; M++)
            {
                for(unsigned N=0; N<DIM; N++)
                {
                    for(unsigned P=0; P<DIM; P++)
                    {
                        for(unsigned Q=0; Q<DIM; Q++)
                        {
                            unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
                            for(unsigned s=0; s<DIM; s++)
                            {
                                mData[index] += matrix(Q,s) * tensor(M,N,P,s);
                            }
                        }
                    }
                }
            }
            break;
        }
        default:
        {
            EXCEPTION("Component not 0, 1, 2, or 3");
        }
    }
}
    
template<unsigned DIM>
double& FourthOrderTensor2<DIM>::operator()(unsigned M, unsigned N, unsigned P, unsigned Q)
{
    assert(M<DIM);
    assert(N<DIM);
    assert(P<DIM);
    assert(Q<DIM);

    unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
    return mData[index];
}

template<unsigned DIM>
void FourthOrderTensor2<DIM>::Zero()
{
    for(unsigned i=0; i<mDimToFour; i++)
    {
        mData[i] = 0.0;
    }
}

#endif //_FOURTHORDERTENSOR2_HPP_
