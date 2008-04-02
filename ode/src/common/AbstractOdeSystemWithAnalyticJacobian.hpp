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

#ifndef _ABSTRACTODESYSTEMWITHANALYTICJACOBIAN_HPP_
#define _ABSTRACTODESYSTEMWITHANALYTICJACOBIAN_HPP_

#include <vector>

#include "Exception.hpp"
#include "AbstractOdeSystem.hpp"

/**
 * Abstract Analytic Jacobian
 *
 * Represents an ODE system with an analytic Jacobian available,
 * which can be computed using the method AnalyticJacobian.
 */
class AbstractOdeSystemWithAnalyticJacobian : public AbstractOdeSystem
{
protected:

public:
    AbstractOdeSystemWithAnalyticJacobian(unsigned numberOfStateVariables = 0)
        : AbstractOdeSystem(numberOfStateVariables)
    {
        mUseAnalytic = true;
    }
    
    /**
     * Compute the analytic Jacobian matrix of the ODE system.
     *
     * @param rSolutionGuess  the current guess at the solution for this time step
     * @param jacobian  will be filled in with the Jacobian matrix entries
     * @param time  the current simulation time
     * @param timeStep  the time step in use by the integrator at present
     */
    virtual void AnalyticJacobian(const std::vector<double>& rSolutionGuess,
                                  double** jacobian,
                                  double time,
                                  double timeStep) = 0;
    
};

#endif //_ABSTRACTODESYSTEMWITHANALYTICJACOBIAN_HPP_
