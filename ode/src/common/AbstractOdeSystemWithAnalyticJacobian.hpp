/**
 * Abstract Analytic Jacobian
*/

#ifndef _ABSTRACTANALYTICJACOBIAN_HPP_
#define _ABSTRACTANALYTICJACOBIAN_HPP_

#include <vector>
#include <string>
#include "Exception.hpp"
#include "AbstractOdeSystem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"

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
