/**
 * Concrete BackwardEulerIvpOdeSolver class.
*/
#ifndef BACKWARDEULERIVPODESOLVER_HPP_
#define BACKWARDEULERIVPODESOLVER_HPP_

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class BackwardEulerIvpOdeSolver : public AbstractOneStepIvpOdeSolver
{
public:
    BackwardEulerIvpOdeSolver()
    {} //Constructor-does nothing
    
    std::vector<double> CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                            double timeStep,
                                            double time,
                                            std::vector<double> currentYValue);
                                            
    virtual ~BackwardEulerIvpOdeSolver()
    {}
    
};


#endif /*BACKWARDEULERIVPODESOLVER_HPP_*/
