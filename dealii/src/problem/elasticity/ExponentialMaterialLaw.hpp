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


#ifndef EXPONENTIALMATERIALLAW_HPP_
#define EXPONENTIALMATERIALLAW_HPP_

#include "AbstractIsotropicIncompressibleMaterialLaw.hpp"
#include "Exception.hpp"

/**
 *  ExponentialMaterialLaw
 *
 *  An exponential isotropic incompressible hyperelastic material law for finite
 *  elasticity
 *
 *  The law is given by a strain energy function
 *      W(I_1,I_2,I_3) = a exp( b(I_1-3) ) - p/2 C^{-1}
 *  in 3d, or
 *      W(I_1,I_2,I_3) = a exp( b(I_1-2) ) - p/2 C^{-1}
 *  in 2d.
 *
 *  Here I_i are the principal invariants of C, the Lagrangian deformation tensor.
 *  (I1=trace(C), I2=trace(C)^2-trace(C^2), I3=det(C)).

 *  Note: only dimension equals 2 or 3 is permitted.
 */

template<unsigned DIM>
class ExponentialMaterialLaw : public AbstractIsotropicIncompressibleMaterialLaw<DIM>
{
private :
    double mA;
    double mB;

public :
    double Get_dW_dI1(double I1, double I2)
    {
        return mA * mB * exp(mB*(I1-DIM));
    }
    double Get_dW_dI2(double I1, double I2)
    {
        // this is covered, but gcov doesn't see this as being covered
        // for some reason, maybe because of optimisations
        #define COVERAGE_IGNORE
        assert(DIM==3);
        #undef COVERAGE_IGNORE

        return 0.0;
    }
    double Get_d2W_dI1(double I1, double I2)
    {
        return mA * mB * mB * exp(mB*(I1-DIM));
    }
    double Get_d2W_dI2(double I1, double I2)
    {
        // this is covered, but gcov doesn't see this as being covered
        // for some reason, maybe because of optimisations
        #define COVERAGE_IGNORE
        assert(DIM==3);
        #undef COVERAGE_IGNORE

        return 0.0;
    }
    double Get_d2W_dI1I2(double I1, double I2)
    {
        // this is covered, but gcov doesn't see this as being covered
        // for some reason, maybe because of optimisations
        #define COVERAGE_IGNORE
        assert(DIM==3);
        #undef COVERAGE_IGNORE

        return 0.0;
    }

    double GetA()
    {
        return mA;
    }
    double GetB()
    {
        return mB;
    }

public :
    /**
     *  Constructor, Taking in the parameters a and b. a must be positive.
     */
    ExponentialMaterialLaw(double a, double b)
    {
        if (DIM!=2 && DIM !=3)
        {
            EXCEPTION("Can only have 2 or 3d incompressible Mooney-Rivlin laws");
        }
        if (a<0.0)
        {
            EXCEPTION("a must be positive");
        }

        mA = a;
        mB = b;
    }
};

#endif /*EXPONENTIALMATERIALLAW_HPP_*/
