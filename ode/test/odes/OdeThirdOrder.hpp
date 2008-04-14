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

/**
 * Concrete OdeThirdOrder class
 */
#ifndef _ODETHIRDORDER_HPP
#define _ODETHIRDORDER_HPP
#include "AbstractOdeSystem.hpp"


class OdeThirdOrder : public AbstractOdeSystem
{
public :
    OdeThirdOrder()
            : AbstractOdeSystem(3) // 3 here is the number of unknowns
    {
        mInitialConditions.push_back(0.0);
        mInitialConditions.push_back(1.0);
        mInitialConditions.push_back(0.0);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0]=rY[0]-rY[1]+rY[2];
        rDY[1]=rY[1]-rY[2];
        rDY[2]=2*rY[1]-rY[2];
    }
};

#endif //_ODETHIRDORDER_HPP
