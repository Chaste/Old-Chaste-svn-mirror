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
 * Concrete Jacobian class
 */
#ifndef _ODEWITHJACOBIAN1_HPP_
#define _ODEWITHJACOBIAN1_HPP_
#include "AbstractOdeSystemWithAnalyticJacobian.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "PetscException.hpp"

class OdeWithJacobian1 : public AbstractOdeSystemWithAnalyticJacobian
{
public :

    OdeWithJacobian1()
            : AbstractOdeSystemWithAnalyticJacobian(1) // 1 here is the number of variables
    {
        mInitialConditions.push_back(0.0);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0]=rY[0]*rY[0];
    }
    
    void AnalyticJacobian(const std::vector<double> &solutionGuess, double** jacobian, double time, double timeStep)
    {
        jacobian[0][0] = 1 - 2.0*timeStep*solutionGuess[0];
    }
    
};


#endif //_ODEWITHJACOBIAN1_HPP_
