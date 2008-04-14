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

#ifndef _ABSTRACTSTIMULUSFUNCTION_HPP_
#define _ABSTRACTSTIMULUSFUNCTION_HPP_
#include <float.h>
/**
 * Represents an abstract stimulus function. Sub-classes will implement the
 * GetStimulus() function to represent the various type of stimuli to the cardiac
 * cell.
 */
class AbstractStimulusFunction
{
public:
    //Returns stimulus at time 'time'
    virtual double GetStimulus(double time) = 0;
    virtual ~AbstractStimulusFunction()
    {}
};

#endif //_ABSTRACTSTIMULUSFUNCTION_HPP_

