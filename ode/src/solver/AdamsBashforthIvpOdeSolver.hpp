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
    AdamsBashforthIvpOdeSolver()
    {}
    ~AdamsBashforthIvpOdeSolver()
    {}
    
    OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem,
                      std::vector<double>& rYValues,
                      double startTime,
                      double endTime,
                      double timeStep,
                      double timeSampling);
                      
    void Solve(AbstractOdeSystem* pAbstractOdeSystem,
               std::vector<double>& rYValues,
               double startTime,
               double endTime,
               double timeStep);
               
};

#endif //_ADAMSBASHFORTHIVPODESOLVER_HPP_

