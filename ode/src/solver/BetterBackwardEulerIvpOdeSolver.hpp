#ifndef BETTERBACKWARDEULERIVPODESOLVER_HPP_
#define BETTERBACKWARDEULERIVPODESOLVER_HPP_

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"
#include <cassert>
#include <vector>

template<int SIZE>
class BetterBackwardEulerIvpOdeSolver  : public AbstractOneStepIvpOdeSolver
{
private:
    /** the epsilon to use in calculating numerical jacobians */
    double mEpsilon;
    bool mForceUseOfNumericalJacobian;
    
    /** Working memory : residual vector */
    double mResidual[SIZE];
    /** Working memory : Jacobian matrix */
    double** mJacobian;
    /** Working memory : update vector */
    double mUpdate[SIZE];
    


protected:
    void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                             double timeStep,
                             double time,
                             std::vector<double>& currentYValues,
                             std::vector<double>& nextYValues)
    {
        // check the size of the ode system matches the solvers expected
        assert(SIZE == pAbstractOdeSystem->GetNumberOfStateVariables());
        
        
        unsigned counter = 0;
//        const double eps = 1e-6 * rCurrentGuess[0]; // Our tolerance (should use min(guess) perhaps?)
        const double eps = 1e-6; // JonW tolerance
        double norm = 2*eps;
        
        std::vector<double> current_guess(SIZE);
        current_guess.assign(currentYValues.begin(), currentYValues.end());
        
        while (norm > eps)
        {
            // Calculate Jacobian and mResidual for current guess
            ComputeResidual(pAbstractOdeSystem, timeStep, time, currentYValues, current_guess);
            ComputeJacobian(pAbstractOdeSystem, timeStep, time, currentYValues, current_guess);
//            // Update norm (our style)
//            norm = ComputeNorm(mResidual);
            
            // Solve Newton linear system
            SolveLinearSystem();
    
            // Update norm (JonW style)
            norm = ComputeNorm(mUpdate);
            
            // Update current guess
            for (unsigned i=0; i<SIZE; i++)
            {
                current_guess[i] -= mUpdate[i];
            }
                    
            assert(counter++ < 15); // avoid infinite loops
        }
        nextYValues.assign(current_guess.begin(), current_guess.end());
        
    }    

public:


    
    BetterBackwardEulerIvpOdeSolver()
    {
        // default epsilon
        mEpsilon = 1e-6;
        mForceUseOfNumericalJacobian = false;
        mJacobian = new double*[SIZE];
        for(unsigned i=0 ; i<SIZE ; i++)
        {
            mJacobian[i] = new double[SIZE];
        }
    }
    
    ~BetterBackwardEulerIvpOdeSolver()
    {
        for(unsigned i=0 ; i<SIZE ; i++)
        {
            delete mJacobian[i];
        }
        delete mJacobian;
    }
    
    void SetEpsilonForNumericalJacobian(double epsilon)
    {
        assert(epsilon > 0);
        mEpsilon = epsilon;
    }
     
    /** Force the solver to use the numerical Jacobian even if the 
     *  ode system is one with an analytical jacobian provided
     */
    void ForceUseOfNumericalJacobian()
    {
        mForceUseOfNumericalJacobian = true;
    }
    
  /* virtual ~BetterBackwardEulerIvpOdeSolver()
    {
    }   
    
  */
    void ComputeResidual(AbstractOdeSystem* pAbstractOdeSystem,
                         double timeStep,
                         double time,
                         std::vector<double>& currentYValues,
                         std::vector<double>& currentGuess)
    {
        std::vector<double> dy(SIZE);//For JC to optimize
        pAbstractOdeSystem->EvaluateYDerivatives(time, currentGuess, dy);
        for(unsigned i=0; i<SIZE; i++)
        {
            
            mResidual[i] = currentGuess[i] - timeStep * dy[i] - currentYValues[i];
        // mResidual = mCurrentGuess - dt* f(mCurrentGuess) - currentYValues
        }                     
    }                    
    

    void ComputeJacobian(AbstractOdeSystem* pAbstractOdeSystem,
                         double timeStep,
                         double time,
                         std::vector<double>& currentYValues,
                         std::vector<double>& currentGuess)
    {
//        Vec current_guess;
//        VecCreate(PETSC_COMM_WORLD, &current_guess);
//        VecSetSizes(current_guess, PETSC_DECIDE,SIZE);
//        //VecSetType(initial_guess, VECSEQ);
//        VecSetFromOptions(current_guess);
//        for(unsigned i=0; i<SIZE; i++)
//        {
//            VecSetValue(current_guess, i, currentGuess[i], INSERT_VALUES);
//        }
//        VecAssemblyBegin(current_guess);
//        VecAssemblyEnd(current_guess);
//        
//        
//        AbstractOdeSystemWithAnalyticJacobian *p_ode_system 
//         = static_cast<AbstractOdeSystemWithAnalyticJacobian*>(pAbstractOdeSystem);
//
//        Mat jacobian;
//        //MatStructure mat_structure;
//        jacobian=MatCreateSeqAIJ(SIZE,SIZE);        
//        MatSetFromOptions(jacobian);
//        
//        p_ode_system->AnalyticJacobian(current_guess, &jacobian, time, timeStep);
//        
//        
//        for(unsigned i=0; i<SIZE; i++)
//        {
//            for(unsigned j=0; j<SIZE; j++)
//            {
//                int row_as_array[1];
//                row_as_array[0] = i;
//                int col_as_array[1];
//                col_as_array[0] = j;
//                
//                double ret_array[1];
//                
//                MatGetValues(jacobian, 1, row_as_array, 1, col_as_array, ret_array);
//                
//                mJacobian[i][j] = ret_array[0];
//            }
//        }
        
        
        for(unsigned i = 0 ; i<SIZE ; i++)
        {
            for(unsigned j = 0 ; j <SIZE ; j++)
            {
                mJacobian[i][j]=0.0;
            }   
        }

        AbstractOdeSystemWithAnalyticJacobian *p_ode_system 
                 = static_cast<AbstractOdeSystemWithAnalyticJacobian*>(pAbstractOdeSystem);

        p_ode_system->BetterAnalyticJacobian(currentGuess, mJacobian, time, timeStep);
        
    }
    
    void SolveLinearSystem()
    {
        double fact;
        for (unsigned i=0; i<SIZE; i++)
        {
            for (unsigned ii=i+1; ii<SIZE; ii++)
            {
                fact = mJacobian[ii][i]/mJacobian[i][i];
                for (unsigned j=i; j<SIZE; j++)
                {
                    mJacobian[ii][j] -= fact*mJacobian[i][j];
                }
                mResidual[ii] -= fact*mResidual[i];
            }
        }
        for (int i=SIZE-1; i>-1; i--)
        {
            mUpdate[i] = mResidual[i];
            for (unsigned j=i+1; j<SIZE; j++)
            {
                mUpdate[i] -= mJacobian[i][j]*mUpdate[j];
            }
            mUpdate[i] /= mJacobian[i][i];
        }
    }
    
    double ComputeNorm(double vector[SIZE])
    {
        double norm = 0.0;
        for (unsigned i=0; i<SIZE; i++)
        {
            if (fabs(vector[i]) > norm)
            {
                norm = fabs(vector[i]);
            }
        }
        return norm;
    }

};


#endif /*BETTERBACKWARDEULERIVPODESOLVER_HPP_*/
