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
    double mFrequency;
    double mStartTime;
    
public:
    RegularStimulus(double magnitudeOfStimulus, double duration, double frequency, double startTime);
    ~RegularStimulus();
    double GetStimulus(double time);
    
};
#endif //_REGULARSTIMULUS_HPP_
