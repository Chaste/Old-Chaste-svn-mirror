

#include "MockEulerIvpOdeSolver.hpp"

MockEulerIvpOdeSolver::MockEulerIvpOdeSolver() : EulerIvpOdeSolver()
{
    mCallCount=0;
}

OdeSolution MockEulerIvpOdeSolver::Solve(AbstractOdeSystem *pAbstractOdeSystem,
                                         std::vector<double> &rYValues,
                                         double startTime,
                                         double endTime,
                                         double timeStep,
                                         double timeSampling)
{
    return EulerIvpOdeSolver::Solve(pAbstractOdeSystem,
                                    rYValues,
                                    startTime,
                                    endTime,
                                    timeStep,
                                    timeSampling);
}

void MockEulerIvpOdeSolver::Solve(AbstractOdeSystem *pAbstractOdeSystem,
                                  std::vector<double> &rYValues,
                                  double startTime,
                                  double endTime,
                                  double timeStep)
{
    mCallCount++;
    EulerIvpOdeSolver::Solve(pAbstractOdeSystem,
                             rYValues,
                             startTime,
                             endTime,
                             timeStep);
}

int MockEulerIvpOdeSolver::GetCallCount()
{
    return mCallCount;
}
