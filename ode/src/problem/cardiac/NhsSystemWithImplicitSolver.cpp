
#include "NhsSystemWithImplicitSolver.hpp"
#include "TimeStepper.hpp"
#include <cmath>

/*
 *=========================== PRIVATE METHODS ==============================
 */

void NhsSystemWithImplicitSolver::ImplicitSolveForActiveTension()
{
    double current_active_tension = mActiveTensionInitialGuess;
    double residual = CalcActiveTensionResidual(current_active_tension);

    unsigned counter = 0;
    while ((fabs(residual)>mTolerance) && (counter++<15))
    {
        double h = std::max(fabs(current_active_tension/100),1e-8);

        double jac = (CalcActiveTensionResidual(current_active_tension+h) - CalcActiveTensionResidual(current_active_tension))/h;
        
        current_active_tension -= residual/jac;
        residual = CalcActiveTensionResidual(current_active_tension);
    } 
    assert(counter<15);
    
    // save the active tension initial guess for next time round
    mActiveTensionInitialGuess = current_active_tension;

    //// could have something like this, an an overloaded GetActiveTension(),
    //// but it will be the same as it parent::GetActiveTension is called after
    //// the state vars have been updated
    mActiveTensionSolution = current_active_tension; 
}

double NhsSystemWithImplicitSolver::CalcActiveTensionResidual(double activeTensionGuess)
{
    // solve for new Ca_trop implicitly
    double new_Ca_Trop = ImplicitSolveForCaTrop(activeTensionGuess);
    
    mTempStoredStateVariables[0] = new_Ca_Trop;
    
    // solve for new Ca_trop implicitly (or with a mixed implicit-explicit scheme, if asked for)
    double new_z;
    if(!mUseImplicitExplicitSolveForZ)
    {
        new_z = ImplicitSolveForZ(new_Ca_Trop);
    }
    else
    {
        new_z = ImplicitExplicitSolveForZ(new_Ca_Trop);
    }
    
    mTempStoredStateVariables[1] = new_z;

    // get the new T0 
    double new_T0 = CalculateT0(new_z); 

    // solve for the new Qi's (and therefore Q) implicitly
    double new_Q = ImplicitSolveForQ();  // should be moved up
                
    double new_active_tension = 0;
    if(new_Q > 0)
    {
        new_active_tension = new_T0*(1+(2+mA)*new_Q)/(1+new_Q);
    }
    else
    {
        new_active_tension = new_T0*(1+mA*new_Q)/(1-new_Q);
    }
    
    return new_active_tension - activeTensionGuess;
}

        

double NhsSystemWithImplicitSolver::ImplicitSolveForCaTrop(double newActiveTension)
{
    double old_Ca_trop = mCurrentStateVars[0];
    
    double numer = old_Ca_trop + mDt * mKon * mCalciumI  * mCalciumTroponinMax;
    double denom = 1 + mDt*mKon*mCalciumI + mDt*mKrefoff*(1-newActiveTension/(mGamma*mTref));

    return numer/denom;
}


double NhsSystemWithImplicitSolver::ImplicitSolveForZ(double newCaTrop)
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

double NhsSystemWithImplicitSolver::CalcZResidual(double z, double newCaTrop)
{
    double ca_trop_to_ca_trop50_ratio_to_n = pow(newCaTrop/mCalciumTrop50, mN);
    double dzdt =  mAlpha0 * ca_trop_to_ca_trop50_ratio_to_n * (1-z) 
                 - mAlphaR1 * z
                 - mAlphaR2 * pow(z,mNr) / (pow(z,mNr) + pow(mKZ,mNr));

    return z - mCurrentStateVars[1] - mDt*dzdt;                     
}


double NhsSystemWithImplicitSolver::ImplicitExplicitSolveForZ(double newCaTrop)
{
    double ca_trop_to_ca_trop50_ratio_to_n = pow(newCaTrop/mCalciumTrop50, mN);
    double old_z = mCurrentStateVars[1];

    double numer =   old_z + mDt * mAlpha0 * ca_trop_to_ca_trop50_ratio_to_n 
                   - mDt * mAlphaR2 * pow(old_z,mNr) / (pow(old_z,mNr) + pow(mKZ,mNr));
    double denom =   1 + mDt * mAlpha0 * ca_trop_to_ca_trop50_ratio_to_n + mAlphaR1 *mDt;
    
    return numer/denom;
}



double NhsSystemWithImplicitSolver::ImplicitSolveForQ()
{
    double old_Q1 = mCurrentStateVars[2];
    double old_Q2 = mCurrentStateVars[3];
    double old_Q3 = mCurrentStateVars[4];
    
    double new_Q1 = (old_Q1 + mDt*mA1*mDLambda1Dt)/(1 + mAlpha1*mDt);
    double new_Q2 = (old_Q2 + mDt*mA2*mDLambda1Dt)/(1 + mAlpha2*mDt);
    double new_Q3 = (old_Q3 + mDt*mA3*mDLambda1Dt)/(1 + mAlpha3*mDt);
    
    mTempStoredStateVariables[2] = new_Q1;
    mTempStoredStateVariables[3] = new_Q2;
    mTempStoredStateVariables[4] = new_Q3;
    
    return new_Q1 + new_Q2 + new_Q3; 
}

            


/*
 *=========================== PRIVATE METHODS ==============================
 */

NhsSystemWithImplicitSolver::NhsSystemWithImplicitSolver()
    : NhsCellularMechanicsOdeSystem()
{
    mTempStoredStateVariables.resize(GetNumberOfStateVariables());
    mActiveTensionInitialGuess = GetActiveTension();
    mUseImplicitExplicitSolveForZ = false;
    
    mActiveTensionSolution = NhsCellularMechanicsOdeSystem::GetActiveTension();
}

void NhsSystemWithImplicitSolver::SetActiveTensionInitialGuess(double activeTensionInitialGuess)
{
    mActiveTensionInitialGuess = activeTensionInitialGuess;
}


void NhsSystemWithImplicitSolver::SolveDoNotUpdate(double startTime, 
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
        mCurrentStateVars = mTempStoredStateVariables;
    
        stepper.AdvanceOneTimeStep();
    }
}

void NhsSystemWithImplicitSolver::UpdateStateVariables()
{
    for(unsigned i=0; i<mStateVariables.size(); i++)
    {
        mStateVariables[i] = mCurrentStateVars[i];
    }      
}

void NhsSystemWithImplicitSolver::UseImplicitExplicitSolveForZ(bool useImplicitExplicitSolveForZ)
{
    mUseImplicitExplicitSolveForZ = useImplicitExplicitSolveForZ;
}
