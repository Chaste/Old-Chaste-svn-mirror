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
    double Get_dW_dI1(double I1, double I2);

    double Get_dW_dI2(double I1, double I2);

    double Get_d2W_dI1(double I1, double I2);

    double Get_d2W_dI2(double I1, double I2);

    double Get_d2W_dI1I2(double I1, double I2);

    double GetC1();

    double GetC2();

public :
    /**
     *  Constructor, Taking in mooney-rivlin parameters c1 and c2.
     *  Note: c2 is not used if the dimension is 2. Just pass in c1 if 2d.
     */
    MooneyRivlinMaterialLaw(double c1, double c2 = MINUS_LARGE);

    /** Scale the dimensional material parameters */
    void ScaleMaterialParameters(double scaleFactor);
};


#endif /*MOONEYRIVLINMATERIALLAW_HPP_*/
