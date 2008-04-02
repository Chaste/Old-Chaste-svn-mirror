/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MOONEYRIVLINMATERIALLAW_HPP_
#define MOONEYRIVLINMATERIALLAW_HPP_

#include "AbstractIsotropicIncompressibleMaterialLaw.hpp"
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
class MooneyRivlinMaterialLaw : public AbstractIsotropicIncompressibleMaterialLaw<DIM>
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
    MooneyRivlinMaterialLaw(double c1, double c2 = MINUS_LARGE)
    {
        if (DIM!=2 && DIM !=3)
        {
            EXCEPTION("Can only have 2 or 3d incompressible Mooney-Rivlin laws");
        }
        
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


#endif /*MOONEYRIVLINMATERIALLAW_HPP_*/
