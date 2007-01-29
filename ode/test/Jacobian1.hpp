/*
 * Concrete Jacobian class
 */
#ifndef _JACOBIAN1_HPP_
#define _JACOBIAN1_HPP_
#include "AbstractOdeSystemWithAnalyticJacobian.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "PetscException.hpp"

class Jacobian1 : public AbstractOdeSystemWithAnalyticJacobian
{
public :

    Jacobian1()
            : AbstractOdeSystemWithAnalyticJacobian(1) // 1 here is the number of variables
    {
        mInitialConditions.push_back(0.0);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0]=rY[0]*rY[0];
    }
    
    void BetterAnalyticJacobian(std::vector<double> &solutionGuess, double** jacobian, double time, double timeStep)
    {
        jacobian[0][0] = 1 - 2.0*timeStep*solutionGuess[0];
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
        double value = 1 - 2.0*timeStep*current_guess[0];
        int row = 0;
        int col = 0;
        MatSetValue(*pJacobian, row, col, value, INSERT_VALUES);
        
        MatAssemblyBegin(*pJacobian,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(*pJacobian,MAT_FINAL_ASSEMBLY);
        return 0;   
    }
    
        
};


#endif //_JACOBIAN1_HPP_
