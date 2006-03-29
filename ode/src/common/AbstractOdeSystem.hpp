/**
 * Abstract OdeSystem class. Sets up variables and functions for a general ODE system.
*/

#ifndef _ABSTRACTODESYSTEM_HPP_
#define _ABSTRACTODESYSTEM_HPP_

#include <vector>
#include <string>

class AbstractOdeSystem
{ 
    protected:
    
    std::vector<double> mStateVariables;
    std::vector<double> mInitialConditions;
 
    public:
	
    std::vector<std::string> mVariableNames;
    std::vector<std::string> mVariableUnits;


	AbstractOdeSystem() {}; /**< Constructor*/
	
	virtual ~AbstractOdeSystem() {}; /**<  Destructor */  
	
	virtual std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY) = 0;
	
    int GetNumberOfStateVariables()
    { 
        if(mInitialConditions.size()>0)
        {
            return mInitialConditions.size(); 
        }
        else
        {
            return mStateVariables.size(); 
        }
    };
    
    virtual void VerifyVariables(std::vector<double>& odeVars) {}
    
    virtual void SetInitialConditions(std::vector<double>& initialConditions) 
    {
        mInitialConditions=initialConditions;
    }
    
    virtual void SetInitialConditionsComponent(int index, double initialCondition) 
    {
        mInitialConditions[index]=initialCondition;
    }
    
    
    std::vector<double> GetInitialConditions()
    {
        return mInitialConditions;
    }
        
    virtual void SetStateVariables(std::vector<double>& stateVariables) {}
    
    std::vector<double> GetStateVariables()
    {
        return mStateVariables;
    }
    
};

#endif //_ABSTRACTODESYSTEM_HPP_
