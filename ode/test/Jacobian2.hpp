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
        
        // Put dx1/dt = x1^2 + x2^2 in...
        // dR1/dx1
        double value = 1 - 2.0*timeStep*current_guess[0];
        int row = 0;
        int col = 0;
        MatSetValue(*pJacobian, row, col, value, INSERT_VALUES);
        // dR1/dx2
        value = - 2.0*timeStep*current_guess[1];
        row = 0;
        col = 1;
        MatSetValue(*pJacobian, row, col, value, INSERT_VALUES);
        // Put dx2/dt  = x1^2 + 2*x2^2 in...
        // dR2/dx1
        value = - 2.0*timeStep*current_guess[0];
        row = 1;
        col = 0;
        MatSetValue(*pJacobian, row, col, value, INSERT_VALUES);
        // dR2/dx2
        value = 1 - 4.0*timeStep*current_guess[1];
        row = 1;
        col = 1;
        MatSetValue(*pJacobian, row, col, value, INSERT_VALUES);
        
        MatAssemblyBegin(*pJacobian,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(*pJacobian,MAT_FINAL_ASSEMBLY);
        return 0;   
    }
    
        
};


#endif //_JACOBIAN1_HPP_
