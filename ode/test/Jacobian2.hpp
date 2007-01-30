/*
 * Concrete Jacobian class
 */
#ifndef _JACOBIAN2_HPP_
#define _JACOBIAN2_HPP_
#include "AbstractOdeSystemWithAnalyticJacobian.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "PetscException.hpp"

class Jacobian2 : public AbstractOdeSystemWithAnalyticJacobian
{
public :

    Jacobian2()
            : AbstractOdeSystemWithAnalyticJacobian(2) // 1 here is the number of variables
    {
        mInitialConditions.push_back(0.0);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0]=rY[0]*rY[0] + rY[1]*rY[1];
        rDY[1]=rY[0]*rY[0] + 2*rY[1]*rY[1];
    }
    
    void AnalyticJacobian(std::vector<double> &solutionGuess, double** jacobian, double time, double timeStep)
    {
        jacobian[0][0] = 1 - 2.0*timeStep*solutionGuess[0];
        jacobian[0][1] =   - 2.0*timeStep*solutionGuess[1];
        jacobian[1][0] =   - 2.0*timeStep*solutionGuess[0];
        jacobian[1][1] = 1 - 4.0*timeStep*solutionGuess[1];
    } 
        
};


#endif //_JACOBIAN1_HPP_
