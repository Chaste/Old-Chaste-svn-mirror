#ifndef ODE5JACOBIAN_HPP_
#define ODE5JACONIAN_HPP_

#include "AbstractOdeSystemWithAnalyticJacobian.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "PetscException.hpp"


class Ode5Jacobian : public AbstractOdeSystemWithAnalyticJacobian
{
private : 
    double mAlpha;
    
public :
    Ode5Jacobian() : AbstractOdeSystemWithAnalyticJacobian(1)  // 1 here is the number of unknowns
    {
        mInitialConditions.push_back(0.2);
        mAlpha = 100;
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0]=mAlpha*rY[0]*(1-rY[0]);
    }
    
    void AnalyticJacobian(std::vector<double> &solutionGuess, double** jacobian, double time, double timeStep)
    {
        jacobian[0][0] = 1.0 - timeStep*mAlpha + 2.0*timeStep*mAlpha*solutionGuess[0];
    }
    
};

#endif /*ODE5JACOBIAN_HPP_*/
