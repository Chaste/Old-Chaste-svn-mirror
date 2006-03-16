/**
 * Abstract OdeSystem class. Sets up variables and functions for a general ODE system.
*/

#ifndef _ABSTRACTODESYSTEM_HPP_
#define _ABSTRACTODESYSTEM_HPP_

#include <vector>
#include <string>

class AbstractOdeSystem
{
	public:
	
    std::vector<std::string> mVariableNames;
    std::vector<std::string> mVariableUnits;
    std::vector<double> mInitialConditions;

	AbstractOdeSystem() {}; /**< Constructor*/
	
	virtual ~AbstractOdeSystem() {}; /**<  Destructor */  
	
	virtual std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY) = 0;
	
    int GetNumberOfStateVariables()
    { 
        return mInitialConditions.size(); 
    };
    
    virtual void VerifyVariables(std::vector<double>& odeVars) {}
};

#endif //_ABSTRACTODESYSTEM_HPP_
