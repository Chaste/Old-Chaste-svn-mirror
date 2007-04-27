/**
 * Abstract OdeSystem class. Sets up variables and functions for a general ODE system.
*/

#ifndef _ABSTRACTODESYSTEM_HPP_
#define _ABSTRACTODESYSTEM_HPP_

#include <vector>
#include <string>
#include "Exception.hpp"
#include <assert.h>

class AbstractOdeSystem
{
protected:
    unsigned mNumberOfStateVariables;
    std::vector<double> mStateVariables;
    std::vector<double> mInitialConditions;
    std::vector<std::string> mVariableNames;
    std::vector<std::string> mVariableUnits;
    
    bool mUseAnalytic;
    
    
public:
    /**
     * Constructor for an ODE system.
     * 
     * @param numberOfStateVariables  how many ODEs make up the system
     */
    AbstractOdeSystem(unsigned numberOfStateVariables = 0);
    
    /**
     * Virtual destructor since we have virtual methods.
     */
    virtual ~AbstractOdeSystem();
    
    /**
     * Method to evaluate the derivatives of the system.
     * 
     * @param time  the current time
     * @param rY  the current values of the state variables
     * @param rDY  storage for the derivatives of the system; will be filled in on return
     */
    virtual void EvaluateYDerivatives(double time, const std::vector<double> &rY,
                                      std::vector<double> &rDY)=0;
                                      
    unsigned GetNumberOfStateVariables() const
    {
        return mNumberOfStateVariables;
    }
    
    
    void SetInitialConditions(const std::vector<double>& rInitialConditions)
    {
        if (rInitialConditions.size() != mNumberOfStateVariables)
        {
            EXCEPTION("The number of initial conditions must be that of the number of state variables");
        }
        mInitialConditions=rInitialConditions;
    }
    
    void SetInitialConditionsComponent(unsigned index, double initialCondition)
    {
        if ( index >= mNumberOfStateVariables)
        {
            EXCEPTION("Index is greater than the number of state variables");
        }
        mInitialConditions[index]=initialCondition;
    }
    
    
    std::vector<double> GetInitialConditions() const
    {
        return mInitialConditions;
    }
    
    void SetStateVariables(const std::vector<double>& rStateVariables)
    {
        if (  mNumberOfStateVariables != rStateVariables.size() )
        {
            EXCEPTION("The size of the passed in vector must be that of the number of state variables");
        }
        mStateVariables = rStateVariables;
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
     *  negative.td::vector<std::string> mVariableUnits;
     * 
     *  After each timestep the solver will call this method on the ODE to see if 
     *  it should stop there. By default, false is returned here.
     */
    virtual bool CalculateStoppingEvent(double time, const std::vector<double> &rY)
    {
        return false;
    }
    
    bool GetUseAnalytic()
    {
        return mUseAnalytic;
    }
    
    unsigned GetStateVariableNumberByName(std::string name)
    {
        unsigned var_number=0;
        while (var_number != mNumberOfStateVariables && mVariableNames[var_number]!=name)
        {
            var_number++;
        }
        if (var_number == mNumberOfStateVariables)
        {
            EXCEPTION("State variable does not exist");
        }
        return var_number;
    }
    
    double GetStateVariableValueByNumber(unsigned varNumber) const
    {
        assert(varNumber < mNumberOfStateVariables);
        return mStateVariables[varNumber];
    }
    
    std::string GetStateVariableUnitsByNumber(unsigned varNumber) const
    {
        assert(varNumber < mNumberOfStateVariables);
        return mVariableUnits[varNumber];
    }
};


#endif //_ABSTRACTODESYSTEM_HPP_
