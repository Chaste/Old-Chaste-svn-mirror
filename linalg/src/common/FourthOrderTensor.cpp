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

#include "FourthOrderTensor.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
FourthOrderTensor<DIM>::FourthOrderTensor()
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
void FourthOrderTensor<DIM>::SetAsProduct(FourthOrderTensor<DIM>& rTensor, const c_matrix<double,DIM,DIM>& rMatrix, unsigned component)
{
    Zero();

    // messy repeated code but needs to be efficiently over neat..
    switch (component)
    {
        case 0:
        {
            for (unsigned M=0; M<DIM; M++)
            {
                for (unsigned N=0; N<DIM; N++)
                {
                    for (unsigned P=0; P<DIM; P++)
                    {
                        for (unsigned Q=0; Q<DIM; Q++)
                        {
                            unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
                            for (unsigned s=0; s<DIM; s++)
                            {
                                mData[index] += rMatrix(M,s) * rTensor(s,N,P,Q);
                            }
                        }
                    }
                }
            }
            break;
        }
        case 1:
        {
            for (unsigned M=0; M<DIM; M++)
            {
                for (unsigned N=0; N<DIM; N++)
                {
                    for (unsigned P=0; P<DIM; P++)
                    {
                        for (unsigned Q=0; Q<DIM; Q++)
                        {
                            unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
                            for (unsigned s=0; s<DIM; s++)
                            {
                                mData[index] += rMatrix(N,s) * rTensor(M,s,P,Q);
                            }
                        }
                    }
                }
            }
            break;
        }
        case 2:
        {
            for (unsigned M=0; M<DIM; M++)
            {
                for (unsigned N=0; N<DIM; N++)
                {
                    for (unsigned P=0; P<DIM; P++)
                    {
                        for (unsigned Q=0; Q<DIM; Q++)
                        {
                            unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
                            for (unsigned s=0; s<DIM; s++)
                            {
                                mData[index] += rMatrix(P,s) * rTensor(M,N,s,Q);
                            }
                        }
                    }
                }
            }
            break;
        }
        case 3:
        {
            for (unsigned M=0; M<DIM; M++)
            {
                for (unsigned N=0; N<DIM; N++)
                {
                    for (unsigned P=0; P<DIM; P++)
                    {
                        for (unsigned Q=0; Q<DIM; Q++)
                        {
                            unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
                            for (unsigned s=0; s<DIM; s++)
                            {
                                mData[index] += rMatrix(Q,s) * rTensor(M,N,P,s);
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
double& FourthOrderTensor<DIM>::operator()(unsigned M, unsigned N, unsigned P, unsigned Q)
{
    assert(M<DIM);
    assert(N<DIM);
    assert(P<DIM);
    assert(Q<DIM);

    unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
    return mData[index];
}

template<unsigned DIM>
void FourthOrderTensor<DIM>::Zero()
{
    for (unsigned i=0; i<mDimToFour; i++)
    {
        mData[i] = 0.0;
    }
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class FourthOrderTensor<1>;
template class FourthOrderTensor<2>;
template class FourthOrderTensor<3>;
