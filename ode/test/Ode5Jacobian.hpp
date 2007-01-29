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
    
    PetscErrorCode AnalyticJacobian(Vec solutionGuess, Mat *pJacobian, double time, double timeStep)
    {
        int num_equations = mNumberOfStateVariables;
        // Copies info from Petsc vector to a vecor we can use!
        ReplicatableVector solution_guess_replicated;
        solution_guess_replicated.ReplicatePetscVector(solutionGuess);
        std::vector<double> current_guess;
        current_guess.reserve(num_equations);
        for (int i=0; i<num_equations; i++)
        {
            current_guess.push_back(solution_guess_replicated[i]);
        }
        
        // Put dx1/dt = x1^2 in...
        double value = 1.0 - timeStep*mAlpha + 2.0*timeStep*mAlpha*current_guess[0];
        int row = 0;
        int col = 0;
        MatSetValue(*pJacobian, row, col, value, INSERT_VALUES);
        
        MatAssemblyBegin(*pJacobian,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(*pJacobian,MAT_FINAL_ASSEMBLY);
        return 0;   
    }
    
    void BetterAnalyticJacobian(std::vector<double> &solutionGuess, double** jacobian, double time, double timeStep)
    {
        jacobian[0][0] = 1.0 - timeStep*mAlpha + 2.0*timeStep*mAlpha*solutionGuess[0];
    }
    
};

#endif /*ODE5JACOBIAN_HPP_*/
