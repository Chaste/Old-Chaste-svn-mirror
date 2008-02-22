#include "RegularStimulus.hpp"
#include <cmath>
#include <cassert>

/**
 * Constructor
 */
RegularStimulus::RegularStimulus(double magnitudeOfStimulus, double duration, double period, double startTime)
{
    assert(period > 0);
    assert(period >= mDuration);

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
