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


#ifndef _RKFTESTODE_HPP_
#define _RKFTESTODE_HPP_

#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"

/**
 * Analytic solution to this, for y(0) = 0.5, is y = (t+1)^2 - 0.5*exp(t).
 */
class RkfTestOde : public AbstractOdeSystem
{
public :

    RkfTestOde() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<RkfTestOde>::Instance();
        SetStateVariables(GetInitialConditions());
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] = rY[0] - time*time + 1.0;
    }
};

template<>
void OdeSystemInformation<RkfTestOde>::Initialise()
{
    this->mVariableNames.push_back("Variable 1");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.5);

    this->mInitialised = true;
}

#endif //_RKFTESTODE_HPP_
