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

#ifndef _REGULARSTIMULUS_HPP_
#define _REGULARSTIMULUS_HPP_

#include "AbstractStimulusFunction.hpp"

/**
 * Provides a regular stimulus of magnitude 'magnitudeOfStimulus' at
 * frequency 'frequency' for duration 'duration' starting at time
 * 'startTime'.
 */
class RegularStimulus : public AbstractStimulusFunction
{
private:
    double mMagnitudeOfStimulus;
    double mDuration;
    double mPeriod;
    double mStartTime;
    
public:
    RegularStimulus(double magnitudeOfStimulus, double duration, double period, double startTime);
    ~RegularStimulus();
    double GetStimulus(double time);
    
};
#endif //_REGULARSTIMULUS_HPP_
