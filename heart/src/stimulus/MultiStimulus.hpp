#ifndef MULTISTIMULUS_HPP_
#define MULTISTIMULUS_HPP_

#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class provides a stimulus function which is the 
 * sum of an arbitrary number of stimuli.
 * 
 * After creation it behaves like a ZeroStimulus until
 * any number of stimuli are added.
 */

class MultiStimulus : public AbstractStimulusFunction
{
private:
    std::vector<AbstractStimulusFunction*> mStimuli;
        
public:
    /**
     * Combine a stimulus with the existing ones.
     * 
     * @param pStimulus pointer to the stimulus to be added
     */   
     void AddStimulus(AbstractStimulusFunction* pStimulus);

    /**
     * Get the magnitude of the multiple stimuli at time 'time'
     *
     * @return  Magnitude of stimulus at time 'time'
     */
     double GetStimulus(double time);
};

#endif /*MULTISTIMULUS_HPP_*/
