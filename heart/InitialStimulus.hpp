#ifndef _INITIALSTIMULUS_HPP_
#define _INITIALSTIMULUS_HPP_

#include "AbstractStimulusFunction.hpp"

class InitialStimulus : public AbstractStimulusFunction 
{
     private:
        double mMagnitudeOfStimulus;
        
     public: 
         InitialStimulus(double magnitudeOfStimulus);
        ~InitialStimulus();
         double GetStimulus(double time);
        
};

#endif //_INITIALSTIMULUS_HPP_

