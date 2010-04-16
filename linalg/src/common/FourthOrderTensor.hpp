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
template<unsigned DIM1, unsigned DIM2, unsigned DIM3, unsigned DIM4>
class FourthOrderTensor
{
private:
    std::vector<double> mData;  /**< The components of the tensor. */

    /** Get the index into the mData vector corresponding to this set of indices
      * @param M  first index
      * @param N  second index
      * @param P  third index
      * @param Q  fourth index
      */
    unsigned GetVectorIndex(unsigned M, unsigned N, unsigned P, unsigned Q)
    {
        assert(M<DIM1);
        assert(N<DIM2);
        assert(P<DIM3);
        assert(Q<DIM4);
        return M + DIM1*N + DIM1*DIM2*P + DIM1*DIM2*DIM3*Q;
    }

public:

    /**
     * Constructor.
     */
    FourthOrderTensor();

    /**
     *  Set to be the inner product of a matrix another fourth order tensor, contracting on first component,
     *  ie. sets this tensor to be R, where
     *  R_{abcd} = X_{Na} T_{Nbcd}
     *
     *  @param rMatrix A matrix
     *  @param rTensor A fourth order tensor
     *
     */
    template<unsigned CONTRACTED_DIM>
    void SetAsContractionOnFirstDimension(const c_matrix<double,CONTRACTED_DIM,DIM1>& rMatrix, FourthOrderTensor<CONTRACTED_DIM,DIM2,DIM3,DIM4>& rTensor);


    /**
     *  Set to be the inner product of a matrix another fourth order tensor, contracting on second component,
     *  ie. sets this tensor to be R, where
     *  R_{abcd} = X_{Nb} T_{aNcd}
     *
     *  @param rMatrix A matrix
     *  @param rTensor A fourth order tensor
     *
     */
    template<unsigned CONTRACTED_DIM>
    void SetAsContractionOnSecondDimension(const c_matrix<double,CONTRACTED_DIM,DIM2>& rMatrix, FourthOrderTensor<DIM1,CONTRACTED_DIM,DIM3,DIM4>& rTensor);

    /**
     *  Set to be the inner product of a matrix another fourth order tensor, contracting on third component,
     *  ie. sets this tensor to be R, where
     *  R_{abcd} = X_{Nc} T_{abNd}
     *
     *  @param rMatrix A matrix
     *  @param rTensor A fourth order tensor
     *
     */
    template<unsigned CONTRACTED_DIM>
    void SetAsContractionOnThirdDimension(const c_matrix<double,CONTRACTED_DIM,DIM3>& rMatrix, FourthOrderTensor<DIM1,DIM2,CONTRACTED_DIM,DIM4>& rTensor);

    /**
     *  Set to be the inner product of a matrix another fourth order tensor, contracting on fourth component,
     *  ie. sets this tensor to be R, where
     *  R_{abcd} = X_{Nd} T_{abcN}
     *
     *  @param rMatrix A matrix
     *  @param rTensor A fourth order tensor
     *
     */
    template<unsigned CONTRACTED_DIM>
    void SetAsContractionOnFourthDimension(const c_matrix<double,CONTRACTED_DIM,DIM4>& rMatrix, FourthOrderTensor<DIM1,DIM2,DIM3,CONTRACTED_DIM>& rTensor);


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




///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned DIM1, unsigned DIM2, unsigned DIM3, unsigned DIM4>
FourthOrderTensor<DIM1,DIM2,DIM3,DIM4>::FourthOrderTensor()
{
    unsigned size = DIM1*DIM2*DIM3*DIM4;
    mData.resize(size, 0.0);

    // allocate memory and zero entries
    mData.resize(size, 0.0);
}

template<unsigned DIM1, unsigned DIM2, unsigned DIM3, unsigned DIM4>
template<unsigned CONTRACTED_DIM>
void FourthOrderTensor<DIM1,DIM2,DIM3,DIM4>::SetAsContractionOnFirstDimension(const c_matrix<double,CONTRACTED_DIM,DIM1>& rMatrix, FourthOrderTensor<CONTRACTED_DIM,DIM2,DIM3,DIM4>& rTensor)
{
    Zero();

            for (unsigned M=0; M<DIM1; M++)
            {
                for (unsigned N=0; N<DIM1; N++)
                {
                    for (unsigned P=0; P<DIM1; P++)
                    {
                        for (unsigned Q=0; Q<DIM1; Q++)
                        {
                            unsigned index = GetVectorIndex(M,N,P,Q);
                            for (unsigned s=0; s<DIM1; s++)
                            {
                                mData[index] += rMatrix(M,s) * rTensor(s,N,P,Q);
                            }
                        }
                    }
                }
            }
}


template<unsigned DIM1, unsigned DIM2, unsigned DIM3, unsigned DIM4>
template<unsigned CONTRACTED_DIM>
void FourthOrderTensor<DIM1,DIM2,DIM3,DIM4>::SetAsContractionOnSecondDimension(const c_matrix<double,CONTRACTED_DIM,DIM2>& rMatrix, FourthOrderTensor<DIM1,CONTRACTED_DIM,DIM3,DIM4>& rTensor)
{
    Zero();

            for (unsigned M=0; M<DIM1; M++)
            {
                for (unsigned N=0; N<DIM1; N++)
                {
                    for (unsigned P=0; P<DIM1; P++)
                    {
                        for (unsigned Q=0; Q<DIM1; Q++)
                        {
                            unsigned index = GetVectorIndex(M,N,P,Q);
                            for (unsigned s=0; s<DIM1; s++)
                            {
                                mData[index] += rMatrix(N,s) * rTensor(M,s,P,Q);
                            }
                        }
                    }
                }
            }
}


template<unsigned DIM1, unsigned DIM2, unsigned DIM3, unsigned DIM4>
template<unsigned CONTRACTED_DIM>
void FourthOrderTensor<DIM1,DIM2,DIM3,DIM4>::SetAsContractionOnThirdDimension(const c_matrix<double,CONTRACTED_DIM,DIM3>& rMatrix, FourthOrderTensor<DIM1,DIM2,CONTRACTED_DIM,DIM4>& rTensor)
{
    Zero();

            for (unsigned M=0; M<DIM1; M++)
            {
                for (unsigned N=0; N<DIM1; N++)
                {
                    for (unsigned P=0; P<DIM1; P++)
                    {
                        for (unsigned Q=0; Q<DIM1; Q++)
                        {
                            unsigned index = GetVectorIndex(M,N,P,Q);
                            for (unsigned s=0; s<DIM1; s++)
                            {
                                mData[index] += rMatrix(P,s) * rTensor(M,N,s,Q);
                            }
                        }
                    }
                }
            }
}


template<unsigned DIM1, unsigned DIM2, unsigned DIM3, unsigned DIM4>
template<unsigned CONTRACTED_DIM>
void FourthOrderTensor<DIM1,DIM2,DIM3,DIM4>::SetAsContractionOnFourthDimension(const c_matrix<double,CONTRACTED_DIM,DIM4>& rMatrix, FourthOrderTensor<DIM1,DIM2,DIM3,CONTRACTED_DIM>& rTensor)
{
    Zero();
            for (unsigned M=0; M<DIM1; M++)
            {
                for (unsigned N=0; N<DIM1; N++)
                {
                    for (unsigned P=0; P<DIM1; P++)
                    {
                        for (unsigned Q=0; Q<DIM1; Q++)
                        {
                            unsigned index = GetVectorIndex(M,N,P,Q);
                            for (unsigned s=0; s<DIM1; s++)
                            {
                                mData[index] += rMatrix(Q,s) * rTensor(M,N,P,s);
                            }
                        }
                    }
                }
            }
}

template<unsigned DIM1, unsigned DIM2, unsigned DIM3, unsigned DIM4>
double& FourthOrderTensor<DIM1,DIM2,DIM3,DIM4>::operator()(unsigned M, unsigned N, unsigned P, unsigned Q)
{
    assert(M<DIM1);
    assert(N<DIM2);
    assert(P<DIM3);
    assert(Q<DIM4);

    return mData[GetVectorIndex(M,N,P,Q)];
}

template<unsigned DIM1, unsigned DIM2, unsigned DIM3, unsigned DIM4>
void FourthOrderTensor<DIM1,DIM2,DIM3,DIM4>::Zero()
{
    for (unsigned i=0; i<mData.size(); i++)
    {
        mData[i] = 0.0;
    }
}

#endif //_FOURTHORDERTENSOR_HPP_
