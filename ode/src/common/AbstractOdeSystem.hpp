#ifndef _ABSTRACTODESYSTEM_HPP_
#define _ABSTRACTODESYSTEM_HPP_

#include <vector>
#include <string>
#include "Exception.hpp"
#include <cassert>

/**
 * Abstract OdeSystem class.
 * Sets up variables and functions for a general ODE system.
 */
class AbstractOdeSystem
{
    friend class TestAbstractOdeSystem;
    
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
    
    
    /**
     * This method is used to establish a state varible's position within
     * the vector of state variables of an ODE system. This number can
     * then be used with the methods GetStateVariableValueByNumber and
     * GetStateVariableUnitsByNumber.
     * 
     * @param name The name of a state variable.
     * @return The state variable's position within
     * the vector of state variables associated with the ODE system.
     */
    unsigned GetStateVariableNumberByName(const std::string name);
    
    /**
     * @param varNumber A state variable's position within
     * the vector of state variables associated with the ODE system.
     * @return The units associated with the state variable.
     */
    double GetStateVariableValueByNumber(unsigned varNumber) const
    {
        assert(varNumber < mNumberOfStateVariables);
        return mStateVariables[varNumber];
    }
    
    /**
     * @param varNumber A state variable's position within
     * the vector of state variables associated with the ODE system.
     * @return The value of the state variable.
     */
    std::string GetStateVariableUnitsByNumber(unsigned varNumber) const
    {
        assert(varNumber < mNumberOfStateVariables);
        return mVariableUnits[varNumber];
    }

protected:
    /**
     * Used to include extra debugging information in exception messages.
     * For example,
     *      EXCEPTION(DumpState("Gating variable out of range"));
     */
    std::string DumpState(const std::string& message);
};


#endif //_ABSTRACTODESYSTEM_HPP_
