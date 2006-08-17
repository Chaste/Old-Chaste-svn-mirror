#ifndef BACKWARDEULERIVPODESOLVEREVENTS_H_
#define BACKWARDEULERIVPODESOLVEREVENTS_H_

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"



class BackwardEulerIvpOdeSolverEvents : public BackwardEulerIvpOdeSolver
{
public:
    BackwardEulerIvpOdeSolverEvents()
            : BackwardEulerIvpOdeSolver()
    {} //Constructor-does nothing
    
    virtual ~BackwardEulerIvpOdeSolverEvents()
    {}
    
    bool CalculateStoppingEvent(AbstractOdeSystem* pAbstractOdeSystem, std::vector<double> currentYValue, double time);
    
};

#endif /*BACKWARDEULERIVPODESOLVEREVENTS_H_*/
