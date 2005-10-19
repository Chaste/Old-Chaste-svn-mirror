#ifndef _ABSTRACTMONODOMAINPROBLEMSTIMULUS_HPP_
#define _ABSTRACTMONODOMAINPROBLEMSTIMULUS_HPP_

#include "MonodomainPde.hpp"

/**
 * Abstract class that specifies the stimulus protocol for a
 * monodomain problem.
 */

template<int SPACE_DIM>
class AbstractMonodomainProblemStimulus
{
public:
    /**
     * Method to set the stimulus function on a monodomain pde.
     * The method should call SetStimulusAtNode at the relevant
     * nodes with the approriate stimulus function.
     * N.B. If the stimulus function is created within this
     * method, uses the static keyword to ensure the stimulus
     * function persits through the execution of the simulation.
     * 
     * @param pPde Monodomain Pde to which stimulus should be applied.
     */
    virtual void Apply(MonodomainPde<SPACE_DIM> *pPde) = 0;
};

#endif //_ABSTRACTMONODOMAINPROBLEMSTIMULUS_HPP_
