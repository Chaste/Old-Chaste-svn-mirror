#ifndef _ABSTRACTANALYTICJACOBIAN_HPP_
#define _ABSTRACTANALYTICJACOBIAN_HPP_

//#include "PetscSetupAndFinalize.hpp" this breaks it.
#include <vector>
#include <string>
#include "Exception.hpp"
#include "petscvec.h"
#include "petscmat.h"

#include "AbstractOdeSystem.hpp"
#include "PetscException.hpp"

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
        mUseAnalytic = true;
    }
    
    virtual void AnalyticJacobian(std::vector<double> &solutionGuess, double** jacobian, double time, double timeStep) = 0;
    
};

#endif //_ABSTRACTANALYTICJACOBIAN_HPP_
