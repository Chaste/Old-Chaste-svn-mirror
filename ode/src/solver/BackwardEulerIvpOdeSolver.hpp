/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef BACKWARDEULERIVPODESOLVER_HPP_
#define BACKWARDEULERIVPODESOLVER_HPP_

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "AbstractOdeSystemWithAnalyticJacobian.hpp"
#include "OdeSolution.hpp"
#include <cassert>
#include <vector>

/**
 * A concrete one step ODE solver class that employs the backward Euler 
 * method. This numerical method is implicit and hence unconditionally stable.
 */
class BackwardEulerIvpOdeSolver  : public AbstractOneStepIvpOdeSolver
{
private:

    /** The number of state variables in the ODE system. */
    unsigned mSizeOfOdeSystem;

    /** The epsilon to use in calculating the numerical Jacobian of the ODE system. */
    double mNumericalJacobianEpsilon;

    /**
     * Whether to force the solver to use the numerical Jacobian even if 
     * the ODE system object provides an analytical Jacobian.
     */
    bool mForceUseOfNumericalJacobian;

    /*
     * NOTE: we use (unsafe) double pointers here rather than
     * std::vectors because using std::vectors would lead to a
     * slow down by a factor of about 4.
     */

    /** Working memory : residual vector */
    double* mResidual;

    /** Working memory : Jacobian matrix */
    double** mJacobian;

    /** Working memory : update vector */
    double* mUpdate;

    /**
     * Compute the current residual.
     * 
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param timeStep  dt
     * @param time  the current time
     * @param rCurrentYValues  the current (initial) state
     * @param rCurrentGuess  current guess for the state at the next timestep
     */
    void ComputeResidual(AbstractOdeSystem* pAbstractOdeSystem,
                         double timeStep,
                         double time,
                         std::vector<double>& rCurrentYValues,
                         std::vector<double>& rCurrentGuess)
    {
        std::vector<double> dy(mSizeOfOdeSystem);//For JC to optimize
        pAbstractOdeSystem->EvaluateYDerivatives(time, rCurrentGuess, dy);
        for (unsigned i=0; i<mSizeOfOdeSystem; i++)
        {

            mResidual[i] = rCurrentGuess[i] - timeStep * dy[i] - rCurrentYValues[i];
        }
    }

    /**
     * Compute the Jacobian of the ODE system.
     * 
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param timeStep  dt
     * @param time  the current time
     * @param rCurrentYValues  the current (initial) state
     * @param rCurrentGuess  current guess for the state at the next timestep
     */
    void ComputeJacobian(AbstractOdeSystem* pAbstractOdeSystem,
                         double timeStep,
                         double time,
                         std::vector<double>& rCurrentYValues,
                         std::vector<double>& rCurrentGuess)
    {
        for (unsigned i=0; i<mSizeOfOdeSystem; i++)
        {
            for (unsigned j=0; j<mSizeOfOdeSystem; j++)
            {
                mJacobian[i][j] = 0.0;
            }
        }

        if (pAbstractOdeSystem->GetUseAnalyticJacobian() && !mForceUseOfNumericalJacobian)
        {
            // The ODE system has an analytic jacobian, so use that
            AbstractOdeSystemWithAnalyticJacobian *p_ode_system
            = static_cast<AbstractOdeSystemWithAnalyticJacobian*>(pAbstractOdeSystem);
            p_ode_system->AnalyticJacobian(rCurrentGuess, mJacobian, time, timeStep);
        }
        else
        {
            ComputeNumericalJacobian(pAbstractOdeSystem,
                                     timeStep,
                                     time,
                                     rCurrentYValues,
                                     rCurrentGuess);
        }
    }

    /**
     * Solve a linear system of equations to update the 
     * current guess for the solution to the ODE system at
     * the next timestep.
     * Used by the method CalculateNextYValue.
     */
    void SolveLinearSystem()
    {
        double fact;
        for (unsigned i=0; i<mSizeOfOdeSystem; i++)
        {
            for (unsigned ii=i+1; ii<mSizeOfOdeSystem; ii++)
            {
                fact = mJacobian[ii][i]/mJacobian[i][i];
                for (unsigned j=i; j<mSizeOfOdeSystem; j++)
                {
                    mJacobian[ii][j] -= fact*mJacobian[i][j];
                }
                mResidual[ii] -= fact*mResidual[i];
            }
        }
        /* This needs to int, since a downloop in unsigned won't
         * terminate properly*/
        for (int i=mSizeOfOdeSystem-1; i>=0; i--)
        {
            mUpdate[i] = mResidual[i];
            for (unsigned j=i+1; j<mSizeOfOdeSystem; j++)
            {
                mUpdate[i] -= mJacobian[i][j]*mUpdate[j];
            }
            mUpdate[i] /= mJacobian[i][i];
        }
    }

    /**
     * Compute the infinity/maximum norm of a vector.
     * Used by the method CalculateNextYValue.
     * 
     * @param vector  a pointer to a vector
     * @return the vector's norm.
     */
    double ComputeNorm(double* vector)
    {
        double norm = 0.0;
        for (unsigned i=0; i<mSizeOfOdeSystem; i++)
        {
            if (fabs(vector[i]) > norm)
            {
                norm = fabs(vector[i]);
            }
        }
        return norm;
    }

    /**
     * Compute the Jacobian of the ODE system numerically.
     * 
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param timeStep  dt
     * @param time  the current time
     * @param rCurrentYValues  the current (initial) state
     * @param rCurrentGuess  current guess for the state at the next timestep
     */
    void ComputeNumericalJacobian(AbstractOdeSystem* pAbstractOdeSystem,
                                  double timeStep,
                                  double time,
                                  std::vector<double>& rCurrentYValues,
                                  std::vector<double>& rCurrentGuess)
    {
        std::vector<double> residual(mSizeOfOdeSystem);
        std::vector<double> residual_perturbed(mSizeOfOdeSystem);
        std::vector<double> guess_perturbed(mSizeOfOdeSystem);

        double epsilon = mNumericalJacobianEpsilon;

        ComputeResidual(pAbstractOdeSystem, timeStep, time, rCurrentYValues, rCurrentGuess);
        for (unsigned i=0; i<mSizeOfOdeSystem; i++)
        {
            residual[i] = mResidual[i];
        }

        for (unsigned global_column=0; global_column<mSizeOfOdeSystem; global_column++)
        {
            for (unsigned i=0; i<mSizeOfOdeSystem; i++)
            {
                guess_perturbed[i] = rCurrentGuess[i];
            }

            guess_perturbed[global_column] += epsilon;

            ComputeResidual(pAbstractOdeSystem, timeStep, time, rCurrentYValues, guess_perturbed);
            for (unsigned i=0; i<mSizeOfOdeSystem; i++)
            {
                residual_perturbed[i] = mResidual[i];
            }

            // Compute residual_perturbed - residual
            double one_over_eps = 1.0/epsilon;
            for (unsigned i=0; i<mSizeOfOdeSystem; i++)
            {
                mJacobian[i][global_column] = one_over_eps*(residual_perturbed[i] - residual[i]);
            }
        }
    }

protected:

    /**
     * Calculate the solution to the ODE system at the next timestep.
     * 
     * A usage example:
     *     BackwardEulerIvpOdeSolver mySolver;
     *     OdeSolution solution = mySolver.Solve(pMyOdeSystem, yInit, StartTime, EndTime, TimeStep, SamplingTime);
     * 
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param timeStep  dt
     * @param time  the current time
     * @param rCurrentYValues  the current (initial) state
     * @param rNextYValues  the state at the next timestep
     */
    void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                             double timeStep,
                             double time,
                             std::vector<double>& rCurrentYValues,
                             std::vector<double>& rNextYValues)
    {
        // check the size of the ode system matches the solvers expected
        assert(mSizeOfOdeSystem == pAbstractOdeSystem->GetNumberOfStateVariables());

        unsigned counter = 0;
//        const double eps = 1e-6 * rCurrentGuess[0]; // Our tolerance (should use min(guess) perhaps?)
        const double eps = 1e-6; // JonW tolerance
        double norm = 2*eps;

        std::vector<double> current_guess(mSizeOfOdeSystem);
        current_guess.assign(rCurrentYValues.begin(), rCurrentYValues.end());

        while (norm > eps)
        {
            // Calculate Jacobian and mResidual for current guess
            ComputeResidual(pAbstractOdeSystem, timeStep, time, rCurrentYValues, current_guess);
            ComputeJacobian(pAbstractOdeSystem, timeStep, time, rCurrentYValues, current_guess);
//            // Update norm (our style)
//            norm = ComputeNorm(mResidual);

            // Solve Newton linear system
            SolveLinearSystem();

            // Update norm (JonW style)
            norm = ComputeNorm(mUpdate);

            // Update current guess
            for (unsigned i=0; i<mSizeOfOdeSystem; i++)
            {
                current_guess[i] -= mUpdate[i];
            }

            counter++;
            assert(counter < 20); // avoid infinite loops
        }
        rNextYValues.assign(current_guess.begin(), current_guess.end());
    }

public:

    /**
     * Constructor.
     * 
     * @param sizeOfOdeSystem  the number of state variables in the ODE system
     */
    BackwardEulerIvpOdeSolver(unsigned sizeOfOdeSystem)
    {
        mSizeOfOdeSystem = sizeOfOdeSystem;

        // default epsilon
        mNumericalJacobianEpsilon = 1e-6;
        mForceUseOfNumericalJacobian = false;

        // allocate memory
        mResidual = new double[mSizeOfOdeSystem];
        mUpdate = new double[mSizeOfOdeSystem];

        mJacobian = new double*[mSizeOfOdeSystem];
        for (unsigned i=0; i<mSizeOfOdeSystem; i++)
        {
            mJacobian[i] = new double[mSizeOfOdeSystem];
        }
    }

    /**
     * Destructor.
     */
    ~BackwardEulerIvpOdeSolver()
    {
        // Delete pointers
        delete[] mResidual;
        delete[] mUpdate;

        for (unsigned i=0; i<mSizeOfOdeSystem; i++)
        {
            delete[] mJacobian[i];
        }
        delete[] mJacobian;
    }

    /**
     * Set the epsilon to use in calculating the 
     * numerical Jacobian of the ODE system.
     * 
     * @param epsilon
     */
    void SetEpsilonForNumericalJacobian(double epsilon)
    {
        assert(epsilon > 0);
        mNumericalJacobianEpsilon = epsilon;
    }

    /**
     * Force the solver to use the numerical Jacobian even if 
     * the ODE system object provides an analytical Jacobian.
     */
    void ForceUseOfNumericalJacobian()
    {
        mForceUseOfNumericalJacobian = true;
    }
};


#endif /*BACKWARDEULERIVPODESOLVER_HPP_*/
