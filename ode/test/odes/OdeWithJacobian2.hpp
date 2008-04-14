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
 * Concrete Jacobian class
 */
#ifndef _ODEWITHJACOBIAN2_HPP_
#define _ODEWITHJACOBIAN2_HPP_
#include "AbstractOdeSystemWithAnalyticJacobian.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "PetscException.hpp"

class OdeWithJacobian2 : public AbstractOdeSystemWithAnalyticJacobian
{
public :

    OdeWithJacobian2()
            : AbstractOdeSystemWithAnalyticJacobian(2) // 1 here is the number of variables
    {
        mInitialConditions.push_back(0.0);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0]=rY[0]*rY[0] + rY[1]*rY[1];
        rDY[1]=rY[0]*rY[0] + 2*rY[1]*rY[1];
    }
    
    void AnalyticJacobian(const std::vector<double> &solutionGuess, double** jacobian, double time, double timeStep)
    {
        jacobian[0][0] = 1 - 2.0*timeStep*solutionGuess[0];
        jacobian[0][1] =   - 2.0*timeStep*solutionGuess[1];
        jacobian[1][0] =   - 2.0*timeStep*solutionGuess[0];
        jacobian[1][1] = 1 - 4.0*timeStep*solutionGuess[1];
    }
    
};


#endif //_ODEWITHJACOBIAN1_HPP_
