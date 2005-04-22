/**
 * Concrete RungeKutta2IvpOdeSolver class. 
*/
#ifndef _RUNGEKUTTA4IVPODESOLVER_HPP_
#define _RUNGEKUTTA4IVPODESOLVER_HPP_

#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class RungeKutta4IvpOdeSolver : public AbstractIvpOdeSolver
{
	public:
	RungeKutta4IvpOdeSolver() {};
	
	OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				      double startTime,
				      double endTime,
				      double timeStep,
				      std::vector<double> initialConditions);
	
};

#endif //_RUNGEKUTTA4IVPODESOLVER_HPP_
