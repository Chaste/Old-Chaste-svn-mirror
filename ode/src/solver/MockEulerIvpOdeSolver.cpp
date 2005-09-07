

#include "MockEulerIvpOdeSolver.hpp"

MockEulerIvpOdeSolver::MockEulerIvpOdeSolver() : EulerIvpOdeSolver()
{
    mCallCount=0;
}

OdeSolution MockEulerIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem, 
                          double startTime,
                          double endTime,
                          double timeStep,
                          std::vector<double> initialConditions)
{
    mCallCount++;
    return EulerIvpOdeSolver::Solve(pAbstractOdeSystem, 
                          startTime,
                          endTime,
                          timeStep,
                          initialConditions);
}

int MockEulerIvpOdeSolver::GetCallCount()
{
    return mCallCount;
}
