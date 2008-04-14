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

#include "RegularStimulus.hpp"
#include <cmath>
#include <cassert>

/**
 * Constructor
 */
RegularStimulus::RegularStimulus(double magnitudeOfStimulus, double duration, double period, double startTime)
{
    assert(period > 0);
    assert(period >= duration);

    mMagnitudeOfStimulus = magnitudeOfStimulus;
    mDuration = duration;
    mPeriod = period;
    mStartTime = startTime;
  
    //Swell duration to avoid rounding issues   
    mDuration += period*DBL_EPSILON;
}

/**
 * Destructor
 */
RegularStimulus::~RegularStimulus()
{
    // Do nothing
}

/**
 * Get the magnitude of stimulus at time 'time'
 *
 * @return  Magnitude of stimulus at time 'time'
 */
double RegularStimulus::GetStimulus(double time)
{
    double beatTime = fmod(time-mStartTime,mPeriod);
    
    if (beatTime >=0 && beatTime <= mDuration)
    {
        return mMagnitudeOfStimulus;
    }
    else
    {
        return 0.0;
    }
}
