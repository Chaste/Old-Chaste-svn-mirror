#ifndef _INITIALSTIMULUS_HPP_
#define _INITIALSTIMULUS_HPP_

#include "AbstractStimulusFunction.hpp"

/**
 * Provides an initial stimulus of magnitude 'magnitudeOfStimulus' from time 0 for duration
 * 'duration'
 */
class InitialStimulus : public AbstractStimulusFunction 
{
     private:
        double mMagnitudeOfStimulus;
        // Duration of initial stimulus
        double mDuration ;
        
     public: 
         InitialStimulus(double magnitudeOfStimulus, double duration);
        ~InitialStimulus();
         double GetStimulus(double time);       
};

#endif //_INITIALSTIMULUS_HPP_

