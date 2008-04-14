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

#ifndef _MOCKEULERIVPODESOLVER_HPP_
#define _MOCKEULERIVPODESOLVER_HPP_

#include "EulerIvpOdeSolver.hpp"

/**
 * This 'mock' class is only used in testing. It is the same
 * as the EulerIvpOdeSolver, but also keeps a count of how many
 * times it has been called. This is useful to check ode solving
 * has been parallelised.
 */


class MockEulerIvpOdeSolver : public EulerIvpOdeSolver
{
private:
    unsigned mCallCount;
    
protected:
    virtual void InternalSolve(AbstractOdeSystem* pAbstractOdeSystem,
                               std::vector<double>& rCurrentYValues,
                               std::vector<double>& rWorkingMemory,
                               double startTime,
                               double endTime,
                               double timeStep);
                               
public:
    MockEulerIvpOdeSolver();
    
    unsigned GetCallCount();
    
    virtual ~MockEulerIvpOdeSolver()
    {}
    
};

#endif //_MOCKEULERIVPODESOLVER_HPP_
