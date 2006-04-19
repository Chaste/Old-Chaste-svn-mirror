/**
 * Concrete EulerIvpOdeSolver class. 
*/
#include "EulerIvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <iostream>
#include <vector>
#include <cassert>

/**
 * Solves a system of ODEs using the Forward Euler method
 * 
 * If @currentYValue isn't given, will update the mStateVariables member of pAbstractOdeSystem
 *   
 * To be used in the form:
 * 
 * EulerIvpOdeSolver mySolver;
 * OdeSolution solution=mySolver->Solve(pMyOdeSystem, StartTime, EndTime, TimeStep, yInit);
 *  
*/

std::vector<double> EulerIvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
										double timeStep,
										double time,
										std::vector<double> currentYValue)
{
	/* 
     * for each timestep in AbstractOneStepIvpSolver calculates a vector containing 
     * the next Y value from the current one for each equation in the system.
	 */

	int num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();

    //\todo this isn't necessary for euler, should have
    //some kind of EvaluateYDerivativesAtStateVariableFunction
    bool bUsingStateVariables = false;
    if(currentYValue.empty())
    {
        bUsingStateVariables = true;
        currentYValue.reserve(num_equations);
        currentYValue = pAbstractOdeSystem->rGetStateVariables();
    }

	std::vector<double> dy(num_equations);
    dy = pAbstractOdeSystem->EvaluateYDerivatives(time, currentYValue);
    
    
    //\todo only reserve m1emory if returning this
	std::vector<double> next_y_value(num_equations);
    
	for(int i=0;i<num_equations; i++) 
	{
		if( !bUsingStateVariables )
        {
            next_y_value[i] = currentYValue[i] + timeStep*dy[i];		
        }
        else
        {
            pAbstractOdeSystem->rGetStateVariables()[i] += timeStep*dy[i]; 
        }
	}
	return next_y_value;
}
