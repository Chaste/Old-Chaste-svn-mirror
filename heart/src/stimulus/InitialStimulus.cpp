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

#include "InitialStimulus.hpp"
#include <cmath>

/**
 * Constructor
 *
 */
InitialStimulus::InitialStimulus(double magnitudeOfStimulus, double duration, double timeOfStimulus)
{
    mMagnitudeOfStimulus = magnitudeOfStimulus;
    mDuration = duration;
    mTimeOfStimulus = timeOfStimulus;
    
    mDuration += (mDuration+mTimeOfStimulus)*DBL_EPSILON; 
}
/**
* Destructor
*/

InitialStimulus::~InitialStimulus()
{
    // Do nothing
}

/**
 * Get the magnitude of stimulus at time 'time'
 *
 * @return  Magnitude of stimulus at time 'time'
 */
double InitialStimulus::GetStimulus(double time)
{
    if (mTimeOfStimulus <= time && time <= mDuration+mTimeOfStimulus)
    {
        return mMagnitudeOfStimulus;
    }
    else
    {
        return 0.0;
    }
}
