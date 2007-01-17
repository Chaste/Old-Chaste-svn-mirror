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
    
public:
    RungeKutta4IvpOdeSolver()
    {};
};

#endif //_RUNGEKUTTA4IVPODESOLVER_HPP_
