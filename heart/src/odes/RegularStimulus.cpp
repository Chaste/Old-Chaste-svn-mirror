#include "RegularStimulus.hpp"
#include <cmath>
#include <cassert>

/**
 * Constructor
 */
RegularStimulus::RegularStimulus(double magnitudeOfStimulus, double duration, double frequency, double startTime)
{
    assert(frequency > 0);
    mMagnitudeOfStimulus = magnitudeOfStimulus;
    mDuration = duration;
    mFrequency = frequency;
    mStartTime = startTime;
  
    double period = 1.0/mFrequency;
    assert(period >= mDuration);
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
    double period = 1.0/mFrequency;
    
    double beatTime = fmod(time-mStartTime,period);
    
    if (beatTime >=0 && beatTime <= mDuration)
    {
        return mMagnitudeOfStimulus;
    }
    else
    {
        return 0.0;
    }
}
