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
 * @return double Magnitude of stimulus at time 'time'
 */
double InitialStimulus::GetStimulus(double time)
{
    if (mTimeOfStimulus <= time  && time <= mDuration+mTimeOfStimulus) 
    {
    	return mMagnitudeOfStimulus;
    } 
    else 
    {
    	return 0.0;
    }
}
