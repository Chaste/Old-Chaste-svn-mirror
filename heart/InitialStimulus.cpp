/**
 * InitialStimulusFunction.cpp
 * 
 */
#include "InitialStimulus.hpp"
#include <cmath>

/**
 * Constructor
 * 
 */
InitialStimulus::InitialStimulus(double magnitudeOfStimulus)
{  
    mMagnitudeOfStimulus = magnitudeOfStimulus;
}
 /**
 * Destructor
 */

InitialStimulus::~InitialStimulus()
{   
    // Do nothing
}

/**
 * Gets the size of stimulus at a given time
 */
double InitialStimulus::GetStimulus(double time)
{
    if (time <=0.5) 
        {
            return mMagnitudeOfStimulus;
        } 
     else 
        {
            return 0.0;
        }
}
