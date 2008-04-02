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

/**
 * Concrete RungeKutta2IvpOdeSolver class.
 */
#ifndef _RUNGEKUTTA4IVPODESOLVER_HPP_
#define _RUNGEKUTTA4IVPODESOLVER_HPP_

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class RungeKutta4IvpOdeSolver : public AbstractOneStepIvpOdeSolver
{
protected:
    void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                             double timeStep,
                             double time,
                             std::vector<double>& currentYValues,
                             std::vector<double>& nextYValues);

private:
    // Working memory
    std::vector<double> k1;
    std::vector<double> k2;
    std::vector<double> k3;
    std::vector<double> k4;
    std::vector<double> yki;
};

#endif //_RUNGEKUTTA4IVPODESOLVER_HPP_
