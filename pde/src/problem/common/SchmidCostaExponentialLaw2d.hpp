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


#ifndef SCHMIDCOSTAEXPONENTIALLAW2D_HPP_
#define SCHMIDCOSTAEXPONENTIALLAW2D_HPP_

#include "AbstractIncompressibleMaterialLaw.hpp"

/**
 *  A 2d version of the material law in Costa, Holmes, McCulloch "Modelling Cardiac
 *  Mechanical Properties in Three Dimensions" Philo. Trans. R. Soc.
 *
 *  W = aexp(Q-1)
 *  where Q = bff*Eff^2 + bfs*Efs^2 + bsf*Esf^2 + bss*Ess^2
 *
 *  where the parameters are taken from the fitting in Schmid,Nash,Young,Hunter
 *  "Myocardial Material Parameter Estimation - A Comparative Study for Simple
 *  Shear" Transactions of the ASME.
 */
class SchmidCostaExponentialLaw2d : public AbstractIncompressibleMaterialLaw<2>
{
private :

    /** Parameter a. */
    double mA;         // should be kPa as the assembler gets the active tension in kPa

    /** Matrix of parameters b. */
    std::vector<std::vector<double> > mB;

    /** 2D identity matrix. */
    c_matrix<double,2,2> mIdentity;

public :

    /** Constructor. */
    SchmidCostaExponentialLaw2d();

    /**
     *  Compute the (2nd Piola Kirchoff) stress T and the stress derivative dT/dE for
     *  a given strain.
     *
     *  NOTE: the strain E is not expected to be passed in, instead the Lagrangian
     *  deformation tensor C is required (recall, E = 0.5(C-I))
     *
     *  dT/dE is a fourth-order tensor, where dT/dE[M][N][P][Q] = dT^{MN}/dE_{PQ}
     *
     *  @param rC The Lagrangian deformation tensor (F^T F)
     *  @param rInvC The inverse of C. Should be computed by the user. (Change this?)
     *  @param pressure the current pressure
     *  @param rT the stress will be returned in this parameter
     *  @param rDTdE the stress derivative will be returned in this parameter, assuming
     *    the final parameter is true
     *  @param computeDTdE a boolean flag saying whether the stress derivative is
     *    required or not.
     */
    void ComputeStressAndStressDerivative(c_matrix<double,2,2>&  rC,
                                          c_matrix<double,2,2>&  rInvC,
                                          double                 pressure,
                                          c_matrix<double,2,2>&  rT,
                                          FourthOrderTensor<2>&  rDTdE,
                                          bool                   computeDTdE);

    /**
     * Get method for mA.
     */
    double GetA();

    /**
     * Get method for mB.
     */
    std::vector<std::vector<double> > GetB();

    /**
     * Get the pressure corresponding to E=0, ie C=identity.
     */
    double GetZeroStrainPressure();
};

#endif /* SCHMIDCOSTAEXPONENTIALLAW2D_HPP_*/
