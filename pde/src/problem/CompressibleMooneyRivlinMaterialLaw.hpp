/*

Copyright (C) University of Oxford, 2005-2011

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


#ifndef COMPRESSIBLEMOONEYRIVLINMATERIALLAW_HPP_
#define COMPRESSIBLEMOONEYRIVLINMATERIALLAW_HPP_

#include "AbstractIsotropicSimpleCompressibleMaterialLaw.hpp"
#include "Exception.hpp"



/**
 *  CompressibleMooneyRivlinMaterialLaw
 *
 *  A Mooney-Rivlin isotropic compressible hyperelastic material law for finite
 *  elasticity
 *
 *  The law is given by a strain energy function
 *      W(I_1,I_2,I_3) = c1(I_1-3) + c2(I_2-3) + c3(I3-1)
 *  in 3d, or
 *      W(I_1,I_3) = c1(I_1-2) + c3(I3-1)
 *  in 2d.
 *
 *  Here I_i are the principal invariants of C, the Lagrangian deformation tensor.
 *  (I1=trace(C), I2=trace(C)^2-trace(C^2), I3=det(C)).
 *
 *  Note c1+c2+c3 must be zero (so zero strain => zero stress).
 */
template<unsigned DIM>
class CompressibleMooneyRivlinMaterialLaw : public AbstractIsotropicSimpleCompressibleMaterialLaw<DIM>
{
private :

    /** Parameter c1. */
    double mC1;

    /** Parameter c2. */
    double mC2;

    /** Parameter c3 */
    double mC3;

public :

    /**
     * Get the first derivative dW/dI1.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_dW_dI1(double I1, double I2)
    {
        return mC1;
    }

    /**
     * Get the first derivative dW/dI2.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_dW_dI2(double I1, double I2)
    {
        return mC2;
    }

    /**
     * Get the second derivative d^2W/dI1^2.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_d2W_dI1(double I1, double I2)
    {
        return 0.0;
    }


    /**
     * Get the second derivative d^2W/dI2^2.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_d2W_dI2(double I1, double I2)
    {
        return 0.0;
    }

    /**
     * Get the second derivative d^2W/dI1dI2.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_d2W_dI1I2(double I1, double I2)
    {
        return 0.0;
    }


    /**
     * Get the first derivative dW/dI3.
     *
     * @param I3 first principal invariant of C (ie det(C))
     */
    double Get_dW_dI3(double I3)
    {
        return mC3;
    }

    /**
     * Get the second derivative d^2W/dI3^2.
     * @param I3 first principal invariant of C (ie det(C))
     */
    double Get_d2W_dI3(double I3)
    {
        return 0.0;
    }

    /** Get method for mC1. */
    double GetC1()
    {
        return mC1;
    }

    /** Get method for mC2. */
    double GetC2()
    {
        return mC2;
    }

    /** Get method for mC3. */
    double GetC3()
    {
        return mC3;
    }

    /**
     * Constructor, taking in Mooney-Rivlin parameters c1, c2 and c3.
     * Note: c2 is not used if the dimension is 2. Just pass in 0.0.
     * c1+c2+c3 must be equal to zero.
     *
     * @param c1 parameter c1
     * @param c2 parameter c2 (should be 0.0 if 2D)
     * @param c3 parameter c3
     */
    CompressibleMooneyRivlinMaterialLaw(double c1, double c2, double c3)
    {
        assert(c1 > 0.0);
        assert(DIM!=2 || c2==0.0);
        mC1 = c1;
        mC2 = c2;
        mC3 = c3;
        if(fabs(c1+c2+c3)>1e-8)
        {
            EXCEPTION("c1+c2+c3 should be equal to zero");
        }
    }

    /**
     * Scale the dimensional material parameters.
     *
     * @param scaleFactor
     */
    void ScaleMaterialParameters(double scaleFactor)
    {
        assert(scaleFactor > 0.0);
        mC1 /= scaleFactor;
        mC2 /= scaleFactor;
        mC3 /= scaleFactor;
    }
};


#endif /*COMPRESSIBLEMOONEYRIVLINMATERIALLAW_HPP_*/
