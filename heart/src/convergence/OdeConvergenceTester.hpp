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

#ifndef ODECONVERGENCETESTER_HPP_
#define ODECONVERGENCETESTER_HPP_

#include "AbstractConvergenceTester.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM, unsigned PROBLEM_DIM>
class OdeConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM, PROBLEM_DIM>
{
public:
    void SetInitialConvergenceParameters()
    {
        this->PdeTimeStep = 2e-2;
        this->OdeTimeStep = 1e-2;
    }
    void UpdateConvergenceParameters()
    {
        this->OdeTimeStep *= 0.5;
    
    }
    bool GiveUpConvergence()
    {
        return this->OdeTimeStep<=1e-8;
    }
    
    double Abscissa()
    {
        return this->OdeTimeStep;
    }
    
};
#endif /*ODECONVERGENCETESTER_HPP_*/
