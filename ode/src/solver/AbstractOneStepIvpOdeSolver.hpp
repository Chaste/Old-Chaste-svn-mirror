/**
 * Abstract One Step Ode Solver class. Sets up variables and functions for all the ODE solvers
 * that only have one timestep.
*/
#ifndef _ABSTRACTONESTEPIVPODESOLVER_HPP_
#define _ABSTRACTONESTEPIVPODESOLVER_HPP_

#include "AbstractIvpOdeSolver.hpp"

class AbstractOneStepIvpOdeSolver : public AbstractIvpOdeSolver
{
	public: 
	virtual OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				              double startTime,
				              double endTime,
				              double timeStep,
				              std::vector<double> initialConditions);
				              
	virtual std::vector<double> CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
								double timeStep,
								double time,
								std::vector<double> currentYValue)=0;
};

#endif //_ABSTRACTONESTEPIVPODESOLVER_HPP_
