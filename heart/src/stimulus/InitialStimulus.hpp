/*

Copyright (C) University of Oxford, 2008

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

