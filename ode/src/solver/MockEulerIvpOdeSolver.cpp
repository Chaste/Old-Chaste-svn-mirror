

#include "MockEulerIvpOdeSolver.hpp"

MockEulerIvpOdeSolver::MockEulerIvpOdeSolver() : EulerIvpOdeSolver()
{
    mCallCount=0;
}

int MockEulerIvpOdeSolver::GetCallCount()
{
    return mCallCount;
}


void MockEulerIvpOdeSolver::InternalSolve(AbstractOdeSystem* pAbstractOdeSystem,
                                          std::vector<double>& rCurrentYValues,
                                          std::vector<double>& rWorkingMemory,
                                          double startTime,
                                          double endTime,
                                          double timeStep)
{
    mCallCount++;
    EulerIvpOdeSolver::InternalSolve(pAbstractOdeSystem,
                                     rCurrentYValues,
                                     rWorkingMemory,
                                     startTime,
                                     endTime,
                                     timeStep);
}
