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

/*
 * Concrete class
 */
#ifndef _RKFTESTODE_HPP_
#define _RKFTESTODE_HPP_
#include "AbstractOdeSystem.hpp"

/*
 * Analytic solution to this, for y(0) = 0.5, is y = (t+1)^2 - 0.5*exp(t).
 */
class RkfTestOde : public AbstractOdeSystem
{
public :

    RkfTestOde() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mInitialConditions.push_back(0.5);
        mStateVariables = mInitialConditions;
    }
    
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0] = rY[0] - time*time + 1.0;
    }
};

#endif //_RKFTESTODE_HPP_
