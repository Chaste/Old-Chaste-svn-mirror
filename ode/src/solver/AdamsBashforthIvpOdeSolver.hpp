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
               
    // prevent this method being called with Adams-Bashforth because it leads
    // to a seg fault
    // TODO: fix. (possibly lots to fix with adams-bashforth)
    void SolveAndUpdateStateVariable(AbstractOdeSystem* pAbstractOdeSystem,
                                     double startTime,
                                     double endTime,
                                     double timeStep)
    {
        EXCEPTION("This method is not implemented with Adams-Bashforth (causes crashes)");
    }
    
};

#endif //_ADAMSBASHFORTHIVPODESOLVER_HPP_

