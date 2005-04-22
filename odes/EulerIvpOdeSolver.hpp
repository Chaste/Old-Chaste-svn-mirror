/**
 * Concrete EulerIvpOdeSolver class. 
*/
#ifndef _EULERIVPODESOLVER_HPP_
#define _EULERIVPODESOLVER_HPP_

#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class EulerIvpOdeSolver : public AbstractIvpOdeSolver
{
	public:
	EulerIvpOdeSolver() {};
	
	OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				      double startTime,
				      double endTime,
				      double timeStep,
				      std::vector<double> initialConditions);
	
};

#endif //_EULERIVPODESOLVER_HPP_

