#ifndef CARDIACNEWTONSOLVER_HPP_
#define CARDIACNEWTONSOLVER_HPP_

#include <cmath>

#include "AbstractBackwardEulerCardiacCell.hpp"

/**
 * Specialised Newton solver for solving the nonlinear systems arising when
 * simulating a cardiac cell using Backward Euler.
 *
 * The class is templated by the size of the nonlinear system, and uses the
 * singleton pattern to ensure only 1 solver for any given system size is created.
 * This allows us to be both computationally and memory efficient.
 *
 * It would be nice to have a test of this class directly, but you need a cardiac
 * cell in order to test it.  So all tests occur when testing particular cardiac
 * cells, e.g. the BackwardEulerLuoRudyIModel1991.
 */
template<unsigned SIZE>
class CardiacNewtonSolver
{
public:
    /**
     * Call this method to obtain a solver instance.
     * 
     * @return a single instance of the class
     */
    static CardiacNewtonSolver<SIZE>* Instance()
    {
        static CardiacNewtonSolver<SIZE> inst;
        return &inst;
    }
    
    /**
     * Use Newton's method to solve the given cell for the next timestep.
     * 
     * @param rCell  the cell to solve
     * @param rCurrentGuess  the current guess at a solution.  Will be mUpdated on exit.
     */
    void Solve(AbstractBackwardEulerCardiacCell<SIZE> &rCell,
               double rCurrentGuess[SIZE])
    {
        unsigned counter = 0;
//        const double eps = 1e-6 * rCurrentGuess[0]; // Our tolerance (should use min(guess) perhaps?)
        const double eps = 1e-6; // JonW tolerance
        double norm = 2*eps;
        
        while (norm > eps)
        {
            // Calculate Jacobian and mResidual for current guess
            rCell.ComputeResidual(rCurrentGuess, mResidual);
            rCell.ComputeJacobian(rCurrentGuess, mJacobian);
            
//            // Update norm (our style)
//            norm = ComputeNorm(mResidual);

            // Solve Newton linear system
            SolveLinearSystem();
            
            // Update norm (JonW style)
            norm = ComputeNorm(mUpdate);
            
            // Update current guess
            for (unsigned i=0; i<SIZE; i++)
            {
                rCurrentGuess[i] -= mUpdate[i];
            }
            
            counter++;
            assert(counter < 15); // avoid infinite loops
        }
    }
    
protected:
    CardiacNewtonSolver()
    {}
    CardiacNewtonSolver(const CardiacNewtonSolver<SIZE>&);
    CardiacNewtonSolver<SIZE>& operator= (const CardiacNewtonSolver<SIZE>&);
    
    /**
     * Compute a norm of a vector.
     */
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
    
    /**
     * Solve a linear system to calculate the Newton update step
     */
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
        /*This must be int, since an unsigned down-loop wouldn't terminate*/
        for (int i=SIZE-1; i>=0; i--)
        {
            mUpdate[i] = mResidual[i];
            for (unsigned j=i+1; j<SIZE; j++)
            {
                mUpdate[i] -= mJacobian[i][j]*mUpdate[j];
            }
            mUpdate[i] /= mJacobian[i][i];
        }
    }
    
private:
    /** Working memory : residual vector */
    double mResidual[SIZE];
    /** Working memory : Jacobian matrix */
    double mJacobian[SIZE][SIZE];
    /** Working memory : update vector */
    double mUpdate[SIZE];
};

#endif /*CARDIACNEWTONSOLVER_HPP_*/
