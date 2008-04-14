/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VANDERPOLODE_HPP_
#define VANDERPOLODE_HPP_

#include "AbstractOdeSystem.hpp"


class VanDerPolOde : public AbstractOdeSystem
{
public :
    VanDerPolOde() : AbstractOdeSystem(2)  // 2 here is the number of unknowns
    {
        mInitialConditions.push_back(10.0);
        mVariableNames.push_back("x");
        mVariableUnits.push_back("m");
        
        mInitialConditions.push_back(10.0);
        mVariableNames.push_back("v");
        mVariableUnits.push_back("m/s");

        SetStateVariables(mInitialConditions);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        double mu = 1;
        rDY[0]= rY[1] + mu*(rY[0] - rY[0]*rY[0]*rY[0]);
        rDY[1] = -rY[0];
    }
    
};


#endif /*VANDERPOLODE_HPP_*/
