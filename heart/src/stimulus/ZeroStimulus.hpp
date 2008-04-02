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

#ifndef ZEROSTIMULUS_HPP_
#define ZEROSTIMULUS_HPP_

#include "AbstractStimulusFunction.hpp"

/**
 *  Stimulus which is always zero. More efficient than using an InitialStimulus
 *  with magnitude zero
 */
class ZeroStimulus : public AbstractStimulusFunction
{
public:
    ZeroStimulus()
    {
    }
    
    virtual ~ZeroStimulus()
    {
    }

    double GetStimulus(double time)
    {
        return 0.0;
    }
};


#endif /*ZEROSTIMULUS_HPP_*/
