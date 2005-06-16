#include "RegularStimulus.hpp"
#include <cmath>
#include <cassert>

/**
 * Constructor
 * 
 */
RegularStimulus::RegularStimulus(double magnitudeOfStimulus, double duration, double frequency, double startTime)
{  
    mMagnitudeOfStimulus = magnitudeOfStimulus;
    mDuration = duration;
    mFrequency = frequency;
    mStartTime = startTime;
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
 * @return double Magnitude of stimulus at time 'time'
 */
double RegularStimulus::GetStimulus(double time)
{
    
    assert(1.0/mFrequency >= mDuration);
    
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
