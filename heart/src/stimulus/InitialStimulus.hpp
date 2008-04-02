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

#ifndef _INITIALSTIMULUS_HPP_
#define _INITIALSTIMULUS_HPP_

#include "AbstractStimulusFunction.hpp"

/**
 * Provides an initial stimulus of magnitude 'magnitudeOfStimulus'
 * from time 'timeOfStimulus' for duration 'duration'.
 */
class InitialStimulus : public AbstractStimulusFunction
{
private:
    /** The stimulus magnitude, typically in microA/cm^2 */
    double mMagnitudeOfStimulus;
    /** Duration of initial stimulus, typically in milliseconds */
    double mDuration;
    /** The time at which the stimulus starts, typically in milliseconds */
    double mTimeOfStimulus;

public:
    InitialStimulus(double magnitudeOfStimulus, double duration, double timeOfStimulus=0.0);
    virtual ~InitialStimulus();
    double GetStimulus(double time);
};

#endif //_INITIALSTIMULUS_HPP_

