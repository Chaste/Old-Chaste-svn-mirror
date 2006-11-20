#ifndef _ABSTRACTANALYTICJACOBIAN_HPP_
#define _ABSTRACTANALYTICJACOBIAN_HPP_

#include "petscvec.h"
#include "petscmat.h"

#include "AbstractOdeSystem.hpp"

/**
 * Abstract Analytic Jacobian
 * 
 * Represents an ODE system with an analytic Jacobian available,
 * which can be computed using the method AnalyticJacobian.
 */
class AbstractOdeSystemWithAnalyticJacobian: public AbstractOdeSystem
{
protected:

    
public:
    AbstractOdeSystemWithAnalyticJacobian(unsigned numberOfStateVariables = 0)
        :AbstractOdeSystem(numberOfStateVariables)
        {
        }
    
    virtual PetscErrorCode AnalyticJacobian(Vec solutionGuess, Mat *pJacobian, double time, double timeStep) = 0;
    
};

#endif //_ABSTRACTANALYTICJACOBIAN_HPP_
