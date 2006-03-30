/**
 * Abstract IvpOdeSolver class. Sets up variables and functions for a numerical solution 
 * technique for an initial value ODE problem. 
*/
#ifndef _ABSTRACTIVPODESOLVER_HPP_
#define _ABSTRACTIVPODESOLVER_HPP_

#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class AbstractIvpOdeSolver
{
	public:
	
        
    /**
     * This version of solve returns an OdeSolution set and takes in initialConditions. 
     */ 
    virtual OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				              double startTime,
				              double endTime,
				              double timeStep,
				              std::vector<double> initialConditions  = std::vector<double>()) = 0;
     
    /**
     * This version of solve modifies the StateVariables member of AbstractOdeSystem
     * instead of returning an OdeSolution set. The StateVariables should be initialised 
     * as the initial conditions. 
     */                          
    virtual void Solve2(AbstractOdeSystem* pAbstractOdeSystem, 
                              double startTime,
                              double endTime,
                              double timeStep) = 0;  
                              
    virtual ~AbstractIvpOdeSolver()
    {
    }                                                         
};

#endif //_ABSTRACTIVPODESOLVER_HPP_
