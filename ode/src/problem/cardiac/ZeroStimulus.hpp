#ifndef ZEROSTIMULUS_HPP_
#define ZEROSTIMULUS_HPP_

#include "AbstractStimulusFunction.hpp"

/**
 *  Stimulus which is always zero. More efficient than using an InitialStimulus
 *  with magnitude zero
 */
class ZeroStimulus : public AbstractStimulusFunction
{
public:
    ZeroStimulus()
    {
    }
    
    virtual ~ZeroStimulus()
    {
    }

    double GetStimulus(double time)
    {
        return 0.0;
    }
};


#endif /*ZEROSTIMULUS_HPP_*/
