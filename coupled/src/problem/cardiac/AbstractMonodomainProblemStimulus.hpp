#ifndef _ABSTRACTMONODOMAINPROBLEMSTIMULUS_HPP_
#define _ABSTRACTMONODOMAINPROBLEMSTIMULUS_HPP_

#include "MonodomainPde.hpp"

template<int SPACE_DIM>
class AbstractMonodomainProblemStimulus
{
public:
    virtual void Apply(MonodomainPde<SPACE_DIM> *pPde) = 0;
};

#endif //_ABSTRACTMONODOMAINPROBLEMSTIMULUS_HPP_
