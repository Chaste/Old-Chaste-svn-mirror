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
    // this is public so that AbstractIvpOdeSolver has direct access
    //\todo - change this to private/protected and make AbstractIvpOdeSolver a
    // friend, or whatever.
    std::vector<double> mStateVariables;
    std::vector<std::string> mVariableNames;
    std::vector<std::string> mVariableUnits;
    
    unsigned int mNumberOfStateVariables;
    std::vector<double> mInitialConditions;
 
    

	
 
public:
    /**
     * Constructor for an ODE system.
     * 
     * @param numberOfStateVariables  how many ODEs make up the system
     */
	AbstractOdeSystem(unsigned numberOfStateVariables = 0)
    {
        mNumberOfStateVariables = numberOfStateVariables;
    }

	virtual ~AbstractOdeSystem() {} /**<  Destructor */  
	
	virtual std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY) = 0;
	
    unsigned GetNumberOfStateVariables()
    { 
        return mNumberOfStateVariables;
    }
    
    
    
    virtual void SetInitialConditions(std::vector<double> initialConditions) 
    {
        assert(initialConditions.size() == mNumberOfStateVariables);
        mInitialConditions=initialConditions;
    }
    
    virtual void SetInitialConditionsComponent(unsigned index, double initialCondition) 
    {
        assert( index < mNumberOfStateVariables); 
        mInitialConditions[index]=initialCondition;
    }
    
    
    std::vector<double> GetInitialConditions()
    {
        return mInitialConditions;
    }
        
    void SetStateVariables(std::vector<double> stateVariables) 
    {
        assert( mNumberOfStateVariables == stateVariables.size() );
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
    
};

#endif //_ABSTRACTODESYSTEM_HPP_
