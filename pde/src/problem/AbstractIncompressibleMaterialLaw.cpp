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

#include "AbstractIncompressibleMaterialLaw.hpp"

template<unsigned DIM>
AbstractIncompressibleMaterialLaw<DIM>::~AbstractIncompressibleMaterialLaw()
{
}

template<unsigned DIM>
void AbstractIncompressibleMaterialLaw<DIM>::ComputeCauchyStress(c_matrix<double,DIM,DIM>& rF,
                                                                 double pressure,
                                                                 c_matrix<double,DIM,DIM>& rSigma)
{
    double detF = Determinant(rF);

    c_matrix<double,DIM,DIM> C = prod(trans(rF), rF);
    c_matrix<double,DIM,DIM> invC = Inverse(C);

    c_matrix<double,DIM,DIM> T;

    static FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE; // not filled in, made static for efficiency

    ComputeStressAndStressDerivative(C, invC, pressure, T, dTdE, false);

    // looping it probably more eficient then doing rSigma = (1/detF)*rF*T*transpose(rF)
    // which doesn't seem to compile anyway, as rF is a Tensor<2,DIM> and T is a
    // SymmetricTensor<2,DIM>
    for (unsigned i=0; i<DIM; i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            rSigma(i,j) = 0.0;
            for (unsigned M=0; M<DIM; M++)
            {
                for (unsigned N=0; N<DIM; N++)
                {
                    rSigma(i,j) += rF(i,M)*T(M,N)*rF(j,N);
                }
            }
            rSigma(i,j) /= detF;
        }
    }
}

template<unsigned DIM>
void AbstractIncompressibleMaterialLaw<DIM>::Compute1stPiolaKirchoffStress(c_matrix<double,DIM,DIM>& rF,
                                                                           double pressure,
                                                                           c_matrix<double,DIM,DIM>& rS)
{
    c_matrix<double,DIM,DIM> C = prod(trans(rF), rF);
    c_matrix<double,DIM,DIM> invC = Inverse(C);

    c_matrix<double,DIM,DIM> T;

    static FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE; // not filled in, made static for efficiency

    ComputeStressAndStressDerivative(C, invC, pressure, T, dTdE, false);

    rS = prod(T, trans(rF));
}

template<unsigned DIM>
void AbstractIncompressibleMaterialLaw<DIM>::Compute2ndPiolaKirchoffStress(c_matrix<double,DIM,DIM>& rC,
                                                                           double pressure,
                                                                           c_matrix<double,DIM,DIM>& rT)
{
    c_matrix<double,DIM,DIM> invC = Inverse(rC);

    static FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE; // not filled in, made static for efficiency

    ComputeStressAndStressDerivative(rC, invC, pressure, rT, dTdE, false);
}

template<unsigned DIM>
void AbstractIncompressibleMaterialLaw<DIM>::ScaleMaterialParameters(double scaleFactor)
{
    #define COVERAGE_IGNORE
    EXCEPTION("[the material law you are using]::ScaleMaterialParameters() has not been implemented\n");
    #undef COVERAGE_IGNORE
}


////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

//template class AbstractIncompressibleMaterialLaw<1>;
template class AbstractIncompressibleMaterialLaw<2>;
template class AbstractIncompressibleMaterialLaw<3>;
