/**
 * Concrete RungeKutta2IvpOdeSolver class.
 */
#ifndef _RUNGEKUTTA4IVPODESOLVER_HPP_
#define _RUNGEKUTTA4IVPODESOLVER_HPP_

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class RungeKutta4IvpOdeSolver : public AbstractOneStepIvpOdeSolver
{
protected:
    void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                             double timeStep,
                             double time,
                             std::vector<double>& currentYValues,
                             std::vector<double>& nextYValues);

private:
    // Working memory
    std::vector<double> k1;
    std::vector<double> k2;
    std::vector<double> k3;
    std::vector<double> k4;
    std::vector<double> yki;
};

#endif //_RUNGEKUTTA4IVPODESOLVER_HPP_
