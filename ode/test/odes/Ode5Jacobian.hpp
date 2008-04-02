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

#ifndef ODE5JACOBIAN_HPP_
#define ODE5JACONIAN_HPP_

#include "AbstractOdeSystemWithAnalyticJacobian.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "PetscException.hpp"


class Ode5Jacobian : public AbstractOdeSystemWithAnalyticJacobian
{
private :
    double mAlpha;
    
public :
    Ode5Jacobian() : AbstractOdeSystemWithAnalyticJacobian(1)  // 1 here is the number of unknowns
    {
        mInitialConditions.push_back(0.2);
        mAlpha = 100;
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0]=mAlpha*rY[0]*(1-rY[0]);
    }
    
    void AnalyticJacobian(const std::vector<double> &solutionGuess, double** jacobian, double time, double timeStep)
    {
        jacobian[0][0] = 1.0 - timeStep*mAlpha + 2.0*timeStep*mAlpha*solutionGuess[0];
    }
    
};

#endif /*ODE5JACOBIAN_HPP_*/
