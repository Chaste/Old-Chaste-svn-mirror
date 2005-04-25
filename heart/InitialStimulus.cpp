#include "InitialStimulus.hpp"
#include <cmath>

/**
 * Constructor
 * 
 */
InitialStimulus::InitialStimulus(double magnitudeOfStimulus, double duration)
{  
    mMagnitudeOfStimulus = magnitudeOfStimulus;
    mDuration = duration;
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
    if (time <= mDuration) 
    {
    	return mMagnitudeOfStimulus;
    } 
    else 
    {
    	return 0.0;
    }
}
