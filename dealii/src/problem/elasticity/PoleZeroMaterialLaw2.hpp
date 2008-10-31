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


#ifndef POLEZEROMATERIALLAW2_HPP_
#define POLEZEROMATERIALLAW2_HPP_


#include "AbstractIncompressibleMaterialLaw2.hpp"
#include "Exception.hpp"



/**
 *  Pole-zero material law
 *
 *  NOTE: CHANGED THE CODE AT THE MOMENT TO NOT USE THE POSITIVE PART BIT BELOW
 *  AS THEN IT IS NOT TWICE DIFFERENTIABLE
 *
 *  W = \Sum_{M,N=1..3} k_{MN} ([ E_{MN} ]_+)^2 / (a_{MN} - E_{MN})^b_{MN}
 *
 *  Note the positive part operator in the numerator, so that
 *  the term in W corresponding to M,N is zero if E_MN < 0. (This differs
 *  from the original pole-zero paper but seems to be what they meant..)
 *
 *  Note that is the parameters k4,k5,k6,a4,a5,a6 etc are known, then
 *  k01=k10=0.5*k4 and similarly with k5,k6, but a01=a10=a4 etc.
 *
 *  Not isotropic, so inherits directly from AbstractIncompressibleMaterialLaw2
 */
template<unsigned DIM>
class PoleZeroMaterialLaw2 : public AbstractIncompressibleMaterialLaw2<DIM>
{
friend class TestMaterialLaws2;

private :
    std::vector<std::vector<double> > mK;
    std::vector<std::vector<double> > mA;
    std::vector<std::vector<double> > mB;

    Tensor<2,DIM> mIdentity;

protected :
    /**
     *  Protected default constructor doing nothing. Just saw inherited classes
     *  can be instantiated and THEN set up the parameters
     */
    PoleZeroMaterialLaw2()
    {
    }

    /**
     *  Set k, a, and b. To be called by the constuctor or a child class
     *  Set comments for constructor.
     */
    void SetParameters(std::vector<std::vector<double> > k,
                       std::vector<std::vector<double> > a,
                       std::vector<std::vector<double> > b)
    {
        if (DIM!=2 && DIM !=3)
        {
            EXCEPTION("Can only have a 2 or 3d incompressible pole-zero law");
        }

        assert(k.size()==DIM);
        assert(a.size()==DIM);
        assert(b.size()==DIM);

        for(unsigned i=0; i<DIM; i++)
        {
            assert(k[i].size()==DIM);
            assert(a[i].size()==DIM);
            assert(b[i].size()==DIM);

            for(unsigned j=0; j<DIM; j++)
            {
                assert( k[i][j] = k[j][i] );
                assert( a[i][j] = a[j][i] );
                assert( b[i][j] = b[j][i] );
            }
        }

        mK = k;
        mA = a;
        mB = b;

        for(unsigned M=0; M<DIM; M++)
        {
            for(unsigned N=0; N<DIM; N++)
            {
                mIdentity[M][N] = M==N ? 1.0 : 0.0;
            }
        }
    }

public :
    /**
     *  Constructor, taking in parameters k_i, a_i, b_i as matrices.
     *  These matrices must be of size DIM-by-DIM and must be symmetric
     *
     *  Note: using the k_1..k_6 convention,  k_4 = 2*k[0][1] = 2*k[1][0], etc
     */
     PoleZeroMaterialLaw2(std::vector<std::vector<double> > k,
                          std::vector<std::vector<double> > a,
                          std::vector<std::vector<double> > b)
    {
        SetParameters(k,a,b);
    }

    void ComputeStressAndStressDerivative(Tensor<2,DIM>&          C,
                                          Tensor<2,DIM>&          invC,
                                          double                  pressure,
                                          SymmetricTensor<2,DIM>& T,
                                          FourthOrderTensor<DIM>& dTdE,
                                          bool                    computeDTdE)
    {
        assert(fabs(C[0][1]-C[1][0]) < 1e-6);

        Tensor<2,DIM> E = 0.5*(C-mIdentity);

        for(unsigned M=0; M<DIM; M++)
        {
            for(unsigned N=0; N<DIM; N++)
            {
                double e = E[M][N];
              //  if(e > 0)
                {
                    double b = mB[M][N];
                    double a = mA[M][N];
                    double k = mK[M][N];

                    //if this fails one of the strain values got too large for the law
                    assert(e < a);

                    T[M][N] =   k
                              * e
                              * (2*(a-e) + b*e)
                              * pow(a-e,-b-1)
                              - pressure*invC[M][N];
                }
//                else
//                {
//                    T[M][N] = 0.0;
//                }
            }
        }

        if(computeDTdE)
        {
            for(unsigned M=0; M<DIM; M++)
            {
                for(unsigned N=0; N<DIM; N++)
                {
                    for(unsigned P=0; P<DIM; P++)
                    {
                        for(unsigned Q=0; Q<DIM; Q++)
                        {
                            dTdE(M,N,P,Q) = 2 * pressure * invC[M][P] * invC[Q][N];
                        }
                    }

                    double e = E[M][N];
                 //   if(e > 0)
                    {
                        double b = mB[M][N];
                        double a = mA[M][N];
                        double k = mK[M][N];

                        dTdE(M,N,M,N) +=   k
                                         * pow(a-e, -b-2)
                                         * (
                                              2*(a-e)*(a-e)
                                            + 4*b*e*(a-e)
                                            + b*(b+1)*e*e
                                           );
                    }
                }
            }
        }
    }

    double GetZeroStrainPressure()
    {
        return 0.0;
    }

    /** Scale the dimensional material parameters (ie the K's) */
    void ScaleMaterialParameters(double scaleFactor)
    {
        assert(scaleFactor > 0.0);
        for(unsigned i=0; i<mK.size(); i++)
        {
            for(unsigned j=0; j<mK[i].size(); j++)
            {
                mK[i][j] /= scaleFactor;
            }
        }
    }
};


#endif /*POLEZEROMATERIALLAW2_HPP_*/
