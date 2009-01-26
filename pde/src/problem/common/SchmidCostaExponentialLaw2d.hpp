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
    double mA;         // should be kPa as the assembler gets the active tension in kPa
    std::vector<std::vector<double> > mB;
    c_matrix<double,2,2> mIdentity;

public :
    SchmidCostaExponentialLaw2d();

    void ComputeStressAndStressDerivative(c_matrix<double,2,2>&  C,
                                          c_matrix<double,2,2>&  invC,
                                          double                 pressure,
                                          c_matrix<double,2,2>&  T,
                                          FourthOrderTensor2<2>& dTdE,
                                          bool                   computeDTdE);

    double GetA();

    std::vector<std::vector<double> > GetB();

    double GetZeroStrainPressure();
};

#endif /* SCHMIDCOSTAEXPONENTIALLAW2D_HPP_*/
