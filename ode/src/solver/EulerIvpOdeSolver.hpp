/**
 * Concrete EulerIvpOdeSolver class.
 */
#ifndef _EULERIVPODESOLVER_HPP_
#define _EULERIVPODESOLVER_HPP_

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class EulerIvpOdeSolver : public AbstractOneStepIvpOdeSolver
{
protected:
    void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                             double timeStep,
                             double time,
                             std::vector<double>& currentYValues,
                             std::vector<double>& nextYValues);
                             
public:
    EulerIvpOdeSolver()
    {}; //Constructor-does nothing
    
    virtual ~EulerIvpOdeSolver()
    {}
};

#endif //_EULERIVPODESOLVER_HPP_

