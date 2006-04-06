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
    * Solves a system of ODEs using a specified one-step ODE solver
    * 
    * @param pAbstractOdeSystem points to the concrete ODE system to be solved
    * @param startTime the time at which the initial conditions are specified
    * @param endTime the time to which the system should be solved and the solution 
    * returned
    * @param timeStep the time interval to be used by the solver
    * @param initialConditions a standard vector specifying the intial condition 
    * of each solution variable in the system. The default is the empty std::vector
    * 
    * @return OdeSolution is an object containing an integer of the number of 
    * equations, a stdAbstractOdeSystem::vector of times and a std::vector of std::vectors where 
    * each of those vectors contains the solution for one variable of the ODE 
    * system at those times.
    * 
    * If @param initialConditions is not given, Solve works with the mStateVariable member
    * of the pAbstractOdeSystem. It uses the current value of mStateVariable as the
    * initial condition, and updates mStateVariable to be the solution at time @param endTime.
    * In this case no output is written to @param OdeSolution
    * 
    */
    virtual OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem, 
                              std::vector<double>& rYValues,
                              double startTime,
                              double endTime,
                              double timeStep,
                              double timeSampling)=0;
    
    virtual void Solve(AbstractOdeSystem* pAbstractOdeSystem, 
                              std::vector<double>& rYValues,
                              double startTime,
                              double endTime,
                              double timeStep)=0;
     
    virtual ~AbstractIvpOdeSolver()
    {
    }
};

#endif //_ABSTRACTIVPODESOLVER_HPP_
