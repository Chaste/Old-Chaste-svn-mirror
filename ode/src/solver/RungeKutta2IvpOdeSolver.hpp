/**
 * Concrete RungeKutta2IvpOdeSolver class.
 */
#ifndef _RUNGEKUTTA2IVPODESOLVER_HPP_
#define _RUNGEKUTTA2IVPODESOLVER_HPP_

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class RungeKutta2IvpOdeSolver : public AbstractOneStepIvpOdeSolver
{
protected:
    void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                             double timeStep,
                             double time,
                             std::vector<double>& currentYValues,
                             std::vector<double>& nextYValues);
    
public:
    RungeKutta2IvpOdeSolver()
    {};	//Constructor-does nothing
                                            
};

#endif //_RUNGEKUTTA2IVPODESOLVER_HPP_
