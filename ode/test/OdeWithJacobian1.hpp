/*
 * Concrete Jacobian class
 */
#ifndef _ODEWITHJACOBIAN1_HPP_
#define _ODEWITHJACOBIAN1_HPP_
#include "AbstractOdeSystemWithAnalyticJacobian.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "PetscException.hpp"

class OdeWithJacobian1 : public AbstractOdeSystemWithAnalyticJacobian
{
public :

    OdeWithJacobian1()
            : AbstractOdeSystemWithAnalyticJacobian(1) // 1 here is the number of variables
    {
        mInitialConditions.push_back(0.0);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0]=rY[0]*rY[0];
    }
    
    void AnalyticJacobian(const std::vector<double> &solutionGuess, double** jacobian, double time, double timeStep)
    {
        jacobian[0][0] = 1 - 2.0*timeStep*solutionGuess[0];
    }
    
};


#endif //_ODEWITHJACOBIAN1_HPP_
