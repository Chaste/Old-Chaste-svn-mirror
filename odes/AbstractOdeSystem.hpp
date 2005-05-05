/**
 * Abstract OdeSystem class. Sets up variables and functions for a general ODE system.
*/

#ifndef _ABSTRACTODESYSTEM_HPP_
#define _ABSTRACTODESYSTEM_HPP_

#include <vector>

class AbstractOdeSystem
{
	private:
	
	int mNumberOfEquations; /**< Number of equations in the ODE system */
	
	public:
	
	int GetNumberOfEquations(void ) {return mNumberOfEquations;}
	
	AbstractOdeSystem(int numberOfEquations): mNumberOfEquations(numberOfEquations) {}; /**< Constructor*/
	
	~AbstractOdeSystem() {}; /**<  Destructor */  
	
	virtual std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY) = 0;
	
};

#endif //_ABSTRACTODESYSTEM_HPP_
