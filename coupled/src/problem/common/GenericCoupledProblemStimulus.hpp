#ifndef _GENERICCOUPLEDPROBLEMSTIMULUS_HPP_
#define _GENERICCOUPLEDPROBLEMSTIMULUS_HPP_

#include "MonodomainPde.hpp"

template<int SPACE_DIM>
class GenericCoupledProblemStimulus
{
public:
    virtual void Apply(MonodomainPde<SPACE_DIM> *pPde) { }
};

#endif //_GENERICCOUPLEDPROBLEMSTIMULUS_HPP_
