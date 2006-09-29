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
protected :
    /**
     * boolean indicating whether the solver quit due to the ODEs
     * stopping event occuring
     */
    bool mStoppingEventOccured;
    
    /** if a stopping event occured the time is stored here */
    double mStoppingTime;
    
    
public :
    /**
     * Solves a system of ODEs using a specified one-step ODE solver
     * 
     * @param pAbstractOdeSystem points to the concrete ODE system to be solved
     * @param startTime the time at which the initial conditions are specified
     * @param endTime the time to which the system should be solved and the solution 
     * returned
     * @param timeStep the time interval to be used by the solver
     * @param initialConditions a standard vector specifying the intial condition 
     * of each solution variable in the system.
     * 
     * @return OdeSolution is an object containing an integer of the number of 
     * equations, a stdAbstractOdeSystem::vector of times and a std::vector of std::vectors where 
     * each of those vectors contains the solution for one variable of the ODE 
     * system at those times.
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
                       
                       
    /**
     * Determine whether the solver quit due to the ODE's stopping event 
     * triggering
     */
    bool StoppingEventOccured()
    {
        return mStoppingEventOccured;
    }
    
    
    AbstractIvpOdeSolver()
            : mStoppingEventOccured(false)
    {}
    
    
    virtual ~AbstractIvpOdeSolver()
    {}
};

#endif //_ABSTRACTIVPODESOLVER_HPP_
