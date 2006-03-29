#ifndef _ABSTRACTMONODOMAINPROBLEMSTIMULUS_HPP_
#define _ABSTRACTMONODOMAINPROBLEMSTIMULUS_HPP_

#include "MonodomainPde.hpp"
#include "ConformingTetrahedralMesh.cpp"

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
     * @param pPde  Monodomain Pde to which stimulus should be applied.
     * @param pMesh The mesh to apply the stimulus on.
     */
    virtual void Apply(MonodomainPde<SPACE_DIM> *pPde,
                       ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> *pMesh) = 0;
    virtual ~AbstractMonodomainProblemStimulus()
    {
    }
};

#endif //_ABSTRACTMONODOMAINPROBLEMSTIMULUS_HPP_
