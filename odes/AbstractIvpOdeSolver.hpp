#ifndef _ABSTRACTIVPODESOLVER_HPP_
#define _ABSTRACTIVPODESOLVER_HPP_

#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class AbstractIvpOdeSolver
{
	public:
	virtual OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				              double StartTime,
				              double EndTime,
				              double TimeStep,
				              std::vector<double> InitialConditions) = 0;
};

#endif //_ABSTRACTIVPODESOLVER_HPP_
