/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#include "RegularStimulus.hpp"
#include <cmath>
#include <cassert>

RegularStimulus::RegularStimulus(double magnitudeOfStimulus, double duration, double period, double startTime, double stopTime)
{
    assert(period > 0);
    assert(period >= duration);

    mMagnitudeOfStimulus = magnitudeOfStimulus;
    mDuration = duration;
    mPeriod = period;
    mStartTime = startTime;
    mStopTime = stopTime;

    //Swell duration to avoid rounding issues
    mDuration += period*DBL_EPSILON;
}

double RegularStimulus::GetStimulus(double time)
{
    double beatTime = fmod(time-mStartTime,mPeriod);

    if (beatTime >=0 && beatTime <= mDuration && time <= mStopTime)
    {
        return mMagnitudeOfStimulus;
    }
    else
    {
        return 0.0;
    }
}

double RegularStimulus::GetStimulusPeriod()
{
    return mPeriod;
}

void RegularStimulus::SetStimulusStartTime(double startTime)
{
    mStartTime = startTime;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RegularStimulus)
