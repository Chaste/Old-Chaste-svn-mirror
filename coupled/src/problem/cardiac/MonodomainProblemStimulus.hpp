#ifndef _MONODOMAINPROBLEMSTIMULUS_HPP_
#define _MONODOMAINPROBLEMSTIMULUS_HPP_

#include "MonodomainPde.hpp"

template<int SPACE_DIM>
class MonodomainProblemStimulus
{
public:
    virtual void Apply(MonodomainPde<SPACE_DIM> *pPde) { }
};

#endif //_MONODOMAINPROBLEMSTIMULUS_HPP_
