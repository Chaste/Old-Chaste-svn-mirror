/**
 * Concrete RungeKutta2IvpOdeSolver class. 
*/
#ifndef _RUNGEKUTTA2IVPODESOLVER_HPP_
#define _RUNGEKUTTA2IVPODESOLVER_HPP_

#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class RungeKutta2IvpOdeSolver : public AbstractIvpOdeSolver
{
	public:
	RungeKutta2IvpOdeSolver() {};
	
	OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				      double startTime,
				      double endTime,
				      double timeStep,
				      std::vector<double> initialConditions);
	
};

#endif //_RUNGEKUTTA2IVPODESOLVER_HPP_
