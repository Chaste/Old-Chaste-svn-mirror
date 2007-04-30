#ifndef SUMSTIMULUS_HPP_
#define SUMSTIMULUS_HPP_

#include "AbstractStimulusFunction.hpp"

/**
 * Provides a stimulus function which is the sum of two such functions.
 */
class SumStimulus : public AbstractStimulusFunction
{
private:
    AbstractStimulusFunction* mpStimulus1;
    AbstractStimulusFunction* mpStimulus2;
    
public:
    SumStimulus(AbstractStimulusFunction *mStimulus1, AbstractStimulusFunction *mStimulus2);
    double GetStimulus(double time);
    
};
#endif /*SUMSTIMULUS_HPP_*/
