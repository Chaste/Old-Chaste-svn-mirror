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
#include "NhsModelWithImplicitSolver.hpp"
#include "TimeStepper.hpp"
#include <cmath>
#include "LogFile.hpp"
#include <iostream>
/*
 *=========================== PRIVATE METHODS ==============================
 */

void NhsModelWithImplicitSolver::ImplicitSolveForActiveTension()
{
    // solve a 1d nonlinear problem for the active tension

    // current active tension
    double current_active_tension = mActiveTensionInitialGuess;

    // see what the residual is
    double residual = CalcActiveTensionResidual(current_active_tension);

    // solve using Newton's method, no damping. Stop if num iterations
    // reaches 15 (very conservative)
    unsigned counter = 0;
    std::vector<double> old_residuals(15);
    while ((fabs(residual)>mTolerance) && (counter++<15))
    {
        // numerically approximate the jacobian
        double h = std::max(fabs(current_active_tension/100),1e-8);

        double jac = (CalcActiveTensionResidual(current_active_tension+h) - CalcActiveTensionResidual(current_active_tension))/h;

        current_active_tension -= residual/jac;
        residual = CalcActiveTensionResidual(current_active_tension);
        old_residuals.push_back(residual);
    }
    //assert(counter<15);
    if(counter >= 15) // not sure what to do here, after having got a case where resid stagnated on +/- 5.500502e-10
    {
        #define COVERAGE_IGNORE
        LOG(1, "\nWARNING in NhsModelWithImplicitSolver::ImplicitSolveForActiveTension(), counter="<<counter<<",resids=");
        for (unsigned i=0; i<old_residuals.size(); i++)
        {
            LOG(1,old_residuals[i]);
        }
        LOG(1,"Final residual = " << residual);
        if(residual > 100*mTolerance)
        {
            LOG(1,"Residual > 100*mTol, throwing exception\n");
            EXCEPTION("NhsModelWithImplicitSolver::ImplicitSolveForActiveTension() failed to converge");
        }
        #undef COVERAGE_IGNORE
    }


    // save the active tension initial guess for next time round
    mActiveTensionInitialGuess = current_active_tension;

    //// could have something like this, an an overloaded GetActiveTension(),
    //// but it will be the same as it parent::GetActiveTension is called after
    //// the state vars have been updated
    mActiveTensionSolution = current_active_tension;
}

double NhsModelWithImplicitSolver::CalcActiveTensionResidual(double activeTensionGuess)
{
    // to calculate the active tension residual, we solve use the current active tension,
    // solve for Ca_trop implicitly, then z implicitly, then Q implicitly, then see
    // what the resulting active tension was

    // solve for new Ca_trop implicitly
    double new_Ca_Trop = ImplicitSolveForCaTrop(activeTensionGuess);

    // store the current Ca_trop solution, in case it turns out to be the current one
    mTemporaryStateVariables[0] = new_Ca_Trop;

    // solve for new z implicitly (or with a mixed implicit-explicit scheme, if asked for)
    double new_z;
    if(!mUseImplicitExplicitSolveForZ)
    {
        new_z = ImplicitSolveForZ(new_Ca_Trop);
    }
    else
    {
        new_z = ImplicitExplicitSolveForZ(new_Ca_Trop);
    }

    // store the current z solution, in case it turns out to be the current one
    mTemporaryStateVariables[1] = new_z;

    // get the new T0
    double new_T0 = CalculateT0(new_z);

    // solve for the new Qi's (and therefore Q) implicitly (note Q1, Q2, Q3 are stored in
    // this method
    double new_Q = ImplicitSolveForQ();                        // <-- SHOULD BE MOVED UP!!
    //double new_Q = QuasiStaticSolveForQ();


    // compute the new active tension and return the difference
    double new_active_tension = 0;
    if(new_Q > 0)
    {
        new_active_tension = new_T0*(1+(2+mA)*new_Q)/(1+new_Q);
    }
    else
    {
        #define COVERAGE_IGNORE // not quite sure how to cover this, of if it is worth it
        new_active_tension = new_T0*(1+mA*new_Q)/(1-new_Q);
        #undef COVERAGE_IGNORE
    }

    return new_active_tension - activeTensionGuess;
}


// Assume the active tension is known and solve for the Ca_trop at the next time
// implicitly using backward euler. This can be done directly as the rhs is linear
// in Ca_trop
double NhsModelWithImplicitSolver::ImplicitSolveForCaTrop(double newActiveTension)
{
    double old_Ca_trop = mCurrentStateVars[0];

    double numer = old_Ca_trop + mDt * mKon * mCalciumI  * mCalciumTroponinMax;
    double denom = 1 + mDt*mKon*mCalciumI + mDt*mKrefoff*(1-newActiveTension/(mGamma*mTref));

    return numer/denom;
}

// Assume the Ca_trop is known and solve for the z at the next time
// implicitly using backward euler. Uses Newton's method
double NhsModelWithImplicitSolver::ImplicitSolveForZ(double newCaTrop)
{
    double current_z = mCurrentStateVars[1];
    double residual = CalcZResidual(current_z, newCaTrop);
    unsigned counter = 0;

    while ((fabs(residual)>mTolerance) && (counter++<15))
    {
        double h = std::max(fabs(current_z/100),1e-8);
        double jac = (CalcZResidual(current_z + h, newCaTrop) - CalcZResidual(current_z, newCaTrop));
        jac/= h;

        current_z -= residual/jac;
        residual = CalcZResidual(current_z, newCaTrop);
    }
    assert(counter<15);

    return current_z;
}


double NhsModelWithImplicitSolver::CalcZResidual(double z, double newCaTrop)
{
    double ca_trop_to_ca_trop50_ratio_to_n = pow(newCaTrop/mCalciumTrop50, mN);
    double dzdt =  mAlpha0 * ca_trop_to_ca_trop50_ratio_to_n * (1-z)
                 - mAlphaR1 * z
                 - mAlphaR2 * pow(z,mNr) / (pow(z,mNr) + pow(mKZ,mNr));

    return z - mCurrentStateVars[1] - mDt*dzdt;
}


// Solve for z semi-implicitly instead of fully implicitly. If we assume we know
// Ca_trop solving for z is a 1d nonlinear problem. Call this to treat the problem
// implicitly in the linear terms on the rhs of dzdt (the (1-z) and (z) terms), and
// explicitly in the nonlinear term (the z^nr/(z^nr + K^nr) term. This means the
// problem can be solved directly and no Newton iterations are needed.
double NhsModelWithImplicitSolver::ImplicitExplicitSolveForZ(double newCaTrop)
{
    double ca_trop_to_ca_trop50_ratio_to_n = pow(newCaTrop/mCalciumTrop50, mN);
    double old_z = mCurrentStateVars[1];

    double numer =   old_z + mDt * mAlpha0 * ca_trop_to_ca_trop50_ratio_to_n
                   - mDt * mAlphaR2 * pow(old_z,mNr) / (pow(old_z,mNr) + pow(mKZ,mNr));
    double denom =   1 + mDt * mAlpha0 * ca_trop_to_ca_trop50_ratio_to_n + mAlphaR1 *mDt;

    return numer/denom;
}


// Solve for Q1, Q2, Q3 implicitly and directly
double NhsModelWithImplicitSolver::ImplicitSolveForQ()
{
    double old_Q1 = mCurrentStateVars[2];
    double old_Q2 = mCurrentStateVars[3];
    double old_Q3 = mCurrentStateVars[4];

    double new_Q1 = (old_Q1 + mDt*mA1*mDLambdaDt)/(1 + mAlpha1*mDt);
    double new_Q2 = (old_Q2 + mDt*mA2*mDLambdaDt)/(1 + mAlpha2*mDt);
    double new_Q3 = (old_Q3 + mDt*mA3*mDLambdaDt)/(1 + mAlpha3*mDt);

    mTemporaryStateVariables[2] = new_Q1;
    mTemporaryStateVariables[3] = new_Q2;
    mTemporaryStateVariables[4] = new_Q3;

    return new_Q1 + new_Q2 + new_Q3;
}


//double NhsModelWithImplicitSolver::QuasiStaticSolveForQ()
//{
//    double new_Q1 = mA1*mDLambdaDt/mAlpha1;
//    double new_Q2 = mA2*mDLambdaDt/mAlpha2;
//    double new_Q3 = mA3*mDLambdaDt/mAlpha3;
//
//    mTemporaryStateVariables[2] = new_Q1;
//    mTemporaryStateVariables[3] = new_Q2;
//    mTemporaryStateVariables[4] = new_Q3;
//
//    return new_Q1 + new_Q2 + new_Q3;
//}



/*
 *=========================== PUBLIC METHODS ==============================
 */

NhsModelWithImplicitSolver::NhsModelWithImplicitSolver()
    : NhsContractionModel()
{
    mActiveTensionInitialGuess = GetActiveTension();
    mUseImplicitExplicitSolveForZ = false;

    mActiveTensionSolution = NhsContractionModel::GetActiveTension();
}

void NhsModelWithImplicitSolver::SetActiveTensionInitialGuess(double activeTensionInitialGuess)
{
    mActiveTensionInitialGuess = activeTensionInitialGuess;
}


void NhsModelWithImplicitSolver::RunDoNotUpdate(double startTime,
                                                 double endTime,
                                                 double timestep)
{
    assert(startTime < endTime);

    mDt = timestep;

    // the implicit routines look in mCurrentStateVars for the current values
    mCurrentStateVars = mStateVariables;

    // loop in time
    TimeStepper stepper(startTime, endTime, timestep);
    while ( !stepper.IsTimeAtEnd() )
    {
        ImplicitSolveForActiveTension();
        mCurrentStateVars = mTemporaryStateVariables;

        stepper.AdvanceOneTimeStep();
    }

    // note: state variables NOT updated.
}


void NhsModelWithImplicitSolver::UseImplicitExplicitSolveForZ(bool useImplicitExplicitSolveForZ)
{
    mUseImplicitExplicitSolveForZ = useImplicitExplicitSolveForZ;
}

double NhsModelWithImplicitSolver::GetNextActiveTension()
{
    return mActiveTensionSolution;
}

