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

#ifndef ODE5_HPP_
#define ODE5_HPP_

#include "AbstractOdeSystem.hpp"


class Ode5 : public AbstractOdeSystem
{
public :
    Ode5() : AbstractOdeSystem(1)  // 1 here is the number of unknowns
    {
        mInitialConditions.push_back(0.2);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        double alpha = 100;
        rDY[0]=alpha*rY[0]*(1-rY[0]);
    }
    
};

#endif /*ODE5_HPP_*/
