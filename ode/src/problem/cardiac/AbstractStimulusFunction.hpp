#ifndef _ABSTRACTSTIMULUSFUNCTION_HPP_
#define _ABSTRACTSTIMULUSFUNCTION_HPP_

/**
 * Represents an abstract stimulus function. Sub-classes will implement the
 * GetStimulus() function to represent the various type of stimuli to the cardiac
 * cell.
 */
class AbstractStimulusFunction
{
public:
    //Returns stimulus at time 'time'
    virtual double GetStimulus(double time) = 0;
    virtual ~AbstractStimulusFunction()
    {}
};

#endif //_ABSTRACTSTIMULUSFUNCTION_HPP_

