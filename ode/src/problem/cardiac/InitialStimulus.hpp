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

