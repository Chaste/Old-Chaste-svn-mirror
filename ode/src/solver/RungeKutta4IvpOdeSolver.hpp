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
public:
    RungeKutta4IvpOdeSolver()
    {};
    
    std::vector<double> CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                            double timeStep,
                                            double time,
                                            std::vector<double> currentYValue);
                                            
};

#endif //_RUNGEKUTTA4IVPODESOLVER_HPP_
