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


#ifndef POLYNOMIALMATERIALLAW3D_HPP_
#define POLYNOMIALMATERIALLAW3D_HPP_

#include "AbstractIsotropicIncompressibleMaterialLaw.hpp"


/**
 *  PolynomialMaterialLaw3d
 *
 *  An incompressible, isotropic, hyperelastic material law with a polynomial form
 *
 *  W(I_1,I_2)  =  \Sigma_{0<p+q<=N}  alpha_{pq} (I_1-3)^p (I_2-3)^q   -  (pressure/2) C^{-1}
 *
 *  For example, if N=1, this reduces to the Mooney Rivlin law
 *     W(I_1,I_2)  =  alpha_{10} (I_1-3) +  alpha_{01} (I_2-3)   -  (pressure/2) C^{-1}
 *  ie the matrix alpha has the form
 *  [ 0  c1 ]
 *  [ c2  0 ]
 *  where c1 and c2 is the usual notation for the Mooney-Rivlin constants
 *
 *  The polynomial is specified by passing in N and the matrix (actually a std::vector
 *  of std::vector<double>s) alpha. alpha should be of size N+1 by N+1, with the bottom
 *  right hand block (ie the components such that p+q>N) all zero. alpha[0][0] should
 *  really also be 0, but, being since alpha[0][0] (I1_3)^0 (I2-3)^0 is a constant
 *  and disappears when the strain energy W is differentiated to obtain the stress, it is
 *  not used. An exception is thrown if alpha[p][q]!=0 for p+q > N though.
 */
class PolynomialMaterialLaw3d : public AbstractIsotropicIncompressibleMaterialLaw<3>
{
private :
    unsigned mN;
    std::vector< std::vector<double> > mAlpha;

public :
    double Get_dW_dI1(double I1, double I2)
    {
        double ret = 0.0;
        // notes: use ints not unsigned as doing p-1
        // (except indexing from p=1 because multiplying by p, but
        // still safer to use ints)
        for (int p=1; p<=(int)mN; p++)
        {
            for (int q=0; q<=(int)mN-p; q++)
            {
                ret += mAlpha[p][q] * p * pow(I1-3,p-1) * pow(I2-3,q);
            }
        }

        return ret;
    }


    double Get_dW_dI2(double I1, double I2)
    {
        double ret = 0.0;
        // notes: use ints not unsigned as doing q-1
        // (except indexing from q=1 because multiplying by q, but
        // still safer to use ints)
        for (int p=0; p<=(int)mN; p++)
        {
            for (int q=1; q<=(int)mN-p; q++)
            {
                ret += mAlpha[p][q] * q * pow(I1-3,p) * pow(I2-3,q-1);
            }
        }
        return ret;
    }


    double Get_d2W_dI1(double I1, double I2)
    {
        double ret = 0.0;

        // notes: use ints not unsigned as doing p-1
        // (except indexing from p=2 because multiplying by p(p-1), but
        // still safer to use ints)
        for (int p=2; p<=(int)mN; p++)
        {
            for (int q=0; q<=(int)mN-p; q++)
            {
                ret += mAlpha[p][q] * p * (p-1) * pow(I1-3,(p-1)*(p-2)) * pow(I2-3,q);
            }
        }
        return ret;
    }


    double Get_d2W_dI2(double I1, double I2)
    {
        double ret = 0.0;

        // notes: use ints not unsigned as doing q-1
        // (except indexing from q=2 because multiplying by q(q-1), but
        // still safer to use ints)
        for (int p=0; p<=(int)mN; p++)
        {
            for (int q=2; q<=(int)mN-p; q++)
            {
                ret += mAlpha[p][q] * q * (q-1) * pow(I1-3,p) * pow(I2-3,(q-1)*(q-2));
            }
        }
        return ret;
    }

    double Get_d2W_dI1I2(double I1, double I2)
    {
        double ret = 0.0;

        // notes: use ints not unsigned as doing p-1
        // (except indexing from p=1,q=1 because multiplying by pq, but
        // still safer to use ints)
        for (int p=1; p<=(int)mN; p++)
        {
            for (int q=1; q<=(int)mN-p; q++)
            {
                ret += mAlpha[p][q] * p * q * pow(I1-3,p-1) * pow(I2-3,q-1);
            }
        }
        return ret;
    }

    double GetAlpha(unsigned i, unsigned j)
    {
        assert(i+j > 0);
        assert(i+j <= mN);

        return mAlpha[i][j];
    }

public :
    PolynomialMaterialLaw3d(unsigned N, std::vector<std::vector<double> > alpha)
    {
        if (N==0)
        {
            EXCEPTION("N must be positive");
        }

        mN = N;

        // error checking: must have alpha[p][q]=0 if p+q>N
        for (unsigned p=0; p<=mN; p++)
        {
            if (alpha[p].size() < mN+1-p)
            {
                EXCEPTION("alpha not big enough");
            }

            for (unsigned q=0; q<alpha[p].size(); q++)
            {
                if ((p+q>mN) && (fabs(alpha[p][q]) > 1e-12))
                {
                    std::stringstream err_mess;
                    err_mess << "alpha[" << p << "][" << q << "] should be zero, as p+q > " << N;
                    EXCEPTION(err_mess.str());
                }
            }
        }

        mAlpha = alpha;
    }

    static std::vector<std::vector<double> > GetZeroedAlpha(unsigned N)
    {
        std::vector<std::vector<double> > alpha(N+1);

        for (unsigned i=0; i<N+1; i++)
        {
            alpha[i].resize(N+1);
            for (unsigned j=0; j<N+1; j++)
            {
                alpha[i][j] = 0.0;
            }
        }

        return alpha;
    }
};



#endif /*POLYNOMIALMATERIALLAW3D_HPP_*/
