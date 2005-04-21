#ifndef _RK2IVPODESOLVER_HPP_
#define _RK2IVPODESOLVER_HPP_

#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class RK2IvpOdeSolver : public AbstractIvpOdeSolver
{
	public:
	RK2IvpOdeSolver() {};
	
	OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				      double startTime,
				      double endTime,
				      double timeStep,
				      std::vector<double> initialConditions);
	
};

#endif //_RK2IVPODESOLVER_HPP_
