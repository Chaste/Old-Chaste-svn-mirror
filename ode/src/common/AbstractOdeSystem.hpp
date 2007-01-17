/**
 * Abstract OdeSystem class. Sets up variables and functions for a general ODE system.
*/

#ifndef _ABSTRACTODESYSTEM_HPP_
#define _ABSTRACTODESYSTEM_HPP_

#include <vector>
#include <string>
#include "Exception.hpp"

class AbstractOdeSystem
{
protected:
    // this is public so that AbstractIvpOdeSolver has direct access
    //\todo - change this to private/protected and make AbstractIvpOdeSolver a
    // friend, or whatever.
    std::vector<double> mStateVariables;
    std::vector<std::string> mVariableNames;
    std::vector<std::string> mVariableUnits;
    
    unsigned int mNumberOfStateVariables;
    std::vector<double> mInitialConditions;
    bool mUseAnalytic;
    
    
public:
    /**
     * Constructor for an ODE system.
     * 
     * @param numberOfStateVariables  how many ODEs make up the system
     */
    AbstractOdeSystem(unsigned numberOfStateVariables = 0)
    {
        mNumberOfStateVariables = numberOfStateVariables;
        mUseAnalytic = false;
    }
    
    virtual ~AbstractOdeSystem()
    {} /**<  Destructor */
    
    /**
     * Method to evaluate the derivatives of the system.
     * 
     * This version will allocate a fresh std::vector<double> to return the results in.
     * 
     * One of the 2 versions of this method MUST be implemented by concrete subclasses.
     * 
     * @param time  the current time
     * @param rY  the current values of the state variables
     * @return  the derivatives of the system
     */
    virtual std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY)
    {
        std::vector<double> res(rY.size());
        EvaluateYDerivatives(time, rY, res);
        return res;
    }

    /**
     * Method to evaluate the derivatives of the system.
     * 
     * This version returns the results in its third parameter.
     * 
     * One of the 2 versions of this method MUST be implemented by concrete subclasses.
     * 
     * @param time  the current time
     * @param rY  the current values of the state variables
     * @param rDY  storage for the derivatives of the system; will be filled in on return
     */
    virtual void EvaluateYDerivatives(double time, const std::vector<double> &rY,
                                      std::vector<double> &rDY)
    {
        std::vector<double> res = EvaluateYDerivatives(time, rY);
        rDY.assign(res.begin(), res.end());
    }    

    unsigned GetNumberOfStateVariables()
    {
        return mNumberOfStateVariables;
    }
    
    
    
    virtual void SetInitialConditions(std::vector<double> initialConditions)
    {
        if (initialConditions.size() != mNumberOfStateVariables)
        {
            EXCEPTION("The number of initial conditions must be that of the number of state variables");
        }
        mInitialConditions=initialConditions;
    }
    
    virtual void SetInitialConditionsComponent(unsigned index, double initialCondition)
    {
        if ( index >= mNumberOfStateVariables)
        {
            EXCEPTION("Index is greater than the number of state variables");
        }
        mInitialConditions[index]=initialCondition;
    }
    
    
    std::vector<double> GetInitialConditions()
    {
        return mInitialConditions;
    }
    
    void SetStateVariables(std::vector<double> stateVariables)
    {
        if ( mNumberOfStateVariables != stateVariables.size() )
        {
            EXCEPTION("The size of the passed in vector must be that of the number of state variables");
        }
        mStateVariables = stateVariables;
    }
    
    std::vector<double>& rGetStateVariables()
    {
        return mStateVariables;
    }
    
    std::vector<std::string>& rGetVariableNames()
    {
        return mVariableNames;
    }
    
    std::vector<std::string>& rGetVariableUnits()
    {
        return mVariableUnits;
    }
    
    /**
     *  CalculateStoppingEvent() - can be overloaded if the ODE is to be solved
     *  only until a particular event (for example, only until the y value becomes
     *  negative.
     * 
     *  After each timestep the solver will call this method on the ODE to see if 
     *  it should stop there. By default, false is returned here.
     */
    virtual bool CalculateStoppingEvent(double time, const std::vector<double> &rY)
    {
        return false;
    }
    
    virtual bool GetUseAnalytic()
    {
        return mUseAnalytic;
    }
    
};

#endif //_ABSTRACTODESYSTEM_HPP_
