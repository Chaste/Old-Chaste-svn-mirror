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

/*
 * Concrete Ode2 class
 */
#ifndef _ODE2_HPP_
#define _ODE2_HPP_
#include "AbstractOdeSystem.hpp"


class Ode2 : public AbstractOdeSystem
{
public :

    Ode2() : AbstractOdeSystem(1) // 1 here is the number of unknowns
    {
        mInitialConditions.push_back(4.0);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0]=rY[0]*time;
    }
};

#endif //_ODE2_HPP_
