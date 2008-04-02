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

#ifndef ZEROSTIMULUSCELLFACTORY_HPP_
#define ZEROSTIMULUSCELLFACTORY_HPP_

#include "AbstractCardiacCellFactory.hpp"

template <class CELL, unsigned DIM>
class ZeroStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{

public:
    ZeroStimulusCellFactory(double timeStep=0.01) : AbstractCardiacCellFactory<DIM>(timeStep)
    {
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        return new CELL(this->mpSolver, this->mTimeStep, this->mpZeroStimulus, this->mpZeroStimulus);
    }
    
    ~ZeroStimulusCellFactory(void)
    {
    }
};

#endif /*ZEROSTIMULUSCELLFACTORY_HPP_*/
