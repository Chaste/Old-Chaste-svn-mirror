/**
 * Abstract One Step Ode Solver class. Sets up variables and functions for all the ODE solvers
 * that only have one timestep.
*/
#ifndef _ABSTRACTONESTEPIVPODESOLVER_HPP_
#define _ABSTRACTONESTEPIVPODESOLVER_HPP_

#include "AbstractIvpOdeSolver.hpp"

class AbstractOneStepIvpOdeSolver : public AbstractIvpOdeSolver
{
    public:
    
    /**
     * This version of solve returns an OdeSolution set and takes in initialConditions. 
     */ 
	virtual OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				              double startTime,
				              double endTime,
				              double timeStep,
				              std::vector<double> initialConditions);
    
    /**
     * This version of solve modifies the StateVariables member of AbstractOdeSystem
     * instead of returning an OdeSolution set. The StateVariables should be initialised 
     * as the initial conditions. 
     */   
    virtual void Solve(AbstractOdeSystem* pAbstractOdeSystem, 
                              double startTime,
                              double endTime,
                              double timeStep);                        
                              
				              
	virtual std::vector<double> CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
								double timeStep,
								double time,
								std::vector<double> currentYValue)=0;
                                
    virtual ~AbstractOneStepIvpOdeSolver()
    {
    }                                
};

#endif //_ABSTRACTONESTEPIVPODESOLVER_HPP_
