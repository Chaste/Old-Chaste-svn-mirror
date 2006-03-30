/**
 * Concrete AdamsBashforthIvpOdeSolver class. 
*/
#ifndef _ADAMSBASHFORTHIVPODESOLVER_HPP_
#define _ADAMSBASHFORTHIVPODESOLVER_HPP_

#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class AdamsBashforthIvpOdeSolver : public AbstractIvpOdeSolver
{
	public:
	AdamsBashforthIvpOdeSolver() {};
	
	OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				      double startTime,
				      double endTime,
				      double timeStep,
				      std::vector<double> initialConditions = std::vector<double>());
	
    void Solve2(AbstractOdeSystem* pAbstractOdeSystem, 
                double startTime,
                double endTime,
                double timeStep);
    
};

#endif //_ADAMSBASHFORTHIVPODESOLVER_HPP_

