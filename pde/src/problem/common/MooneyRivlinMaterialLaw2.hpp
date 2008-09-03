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


#ifndef MOONEYRIVLINMATERIALLAW2_HPP_
#define MOONEYRIVLINMATERIALLAW2_HPP_

#include "AbstractIsotropicIncompressibleMaterialLaw2.hpp"
#include "Exception.hpp"

#define MINUS_LARGE -1e6


/**
 *  MooneyRivlinMaterialLaw
 *
 *  A Mooney-Rivlin isotropic incompressible hyperelastic material law for finite
 *  elasticity
 *
 *  The law is given by a strain energy function
 *      W(I_1,I_2,I_3) = c1(I_1-3) + c2(I_2-3) - p/2 C^{-1}
 *  in 3d, or
 *      W(I_1,I_3) = c1(I_1-2) - p/2 C^{-1}
 *  in 2d.
 *
 *  Here I_i are the principal invariants of C, the Lagrangian deformation tensor.
 *  (I1=trace(C), I2=trace(C)^2-trace(C^2), I3=det(C)).

 *  Note: only dimension equals 2 or 3 is permitted.
 */
template<unsigned DIM>
class MooneyRivlinMaterialLaw2 : public AbstractIsotropicIncompressibleMaterialLaw2<DIM>
{
private :
    double mC1;
    double mC2;

public :
    double Get_dW_dI1(double I1, double I2)
    {
        return mC1;
    }
    double Get_dW_dI2(double I1, double I2)
    {
        // this is covered, but gcov doesn't see this as being covered
        // for some reason, maybe because of optimisations
        #define COVERAGE_IGNORE
        assert(DIM==3);
        #undef COVERAGE_IGNORE
        return mC2;
    }
    double Get_d2W_dI1(double I1, double I2)
    {
        return 0.0;
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

    double GetC1()
    {
        return mC1;
    }

    double GetC2()
    {
        assert(DIM==3);
        return mC2;
    }

public :
    /**
     *  Constructor, Taking in mooney-rivlin parameters c1 and c2.
     *  Note: c2 is not used if the dimension is 2. Just pass in c1 if 2d.
     */
    MooneyRivlinMaterialLaw2(double c1, double c2 = MINUS_LARGE)
    {
        assert(DIM==2 || DIM ==3);

        // if dim==3, check that c2 was passed in, ie c2 isn't the default value
        if ((DIM==3) && (c2<MINUS_LARGE+1))
        {
            EXCEPTION("Two parameters needed for 3d Mooney-Rivlin");
        }

        if (c1 < 0.0)
        {
            EXCEPTION("c1 must be positive in mooney-rivlin"); // is this correct?
        }

        mC1 = c1;
        mC2 = c2;
    }

    /** Scale the dimensional material parameters */
    void ScaleMaterialParameters(double scaleFactor)
    {
        assert(scaleFactor > 0.0);
        mC1 /= scaleFactor;
        mC2 /= scaleFactor;
    }
};


#endif /*MOONEYRIVLINMATERIALLAW2_HPP_*/
