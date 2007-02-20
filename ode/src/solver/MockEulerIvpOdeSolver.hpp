#ifndef _MOCKEULERIVPODESOLVER_HPP_
#define _MOCKEULERIVPODESOLVER_HPP_

#include "EulerIvpOdeSolver.hpp"

/**
 * This 'mock' class is only used in testing. It is the same
 * as the EulerIvpOdeSolver, but also keeps a count of how many
 * times it has been called. This is useful to check ode solving
 * has been parallelised.
 */


class MockEulerIvpOdeSolver : public EulerIvpOdeSolver
{
private:
    unsigned mCallCount;

protected:
    virtual void InternalSolve(AbstractOdeSystem* pAbstractOdeSystem,
                               std::vector<double>& rCurrentYValues,
                               std::vector<double>& rWorkingMemory,
                               double startTime,
                               double endTime,
                               double timeStep);
    
public:
    MockEulerIvpOdeSolver();
    
    unsigned GetCallCount();
                       
    virtual ~MockEulerIvpOdeSolver()
    {}
    
};

#endif //_MOCKEULERIVPODESOLVER_HPP_
