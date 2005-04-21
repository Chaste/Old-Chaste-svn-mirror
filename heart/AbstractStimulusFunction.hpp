#ifndef _ABSTRACTSTIMULUSFUNCTION_HPP_
#define _ABSTRACTSTIMULUSFUNCTION_HPP_


class AbstractStimulusFunction 
{
    
    public:
        //Function that returns stimulus at time r
        virtual double GetStimulus(double time) = 0;
};

#endif //_ABSTRACTSTIMULUSFUNCTION_HPP_

