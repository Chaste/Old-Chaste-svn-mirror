#ifndef NHSSYSTEMWITHIMPLICITSOLVER_HPP_
#define NHSSYSTEMWITHIMPLICITSOLVER_HPP_

#include "NHSCellularMechanicsOdeSystem.hpp"
#include "TimeStepper.hpp"

// todo: doxygen

/** 
 *  NHS system with build in implicit solver. Jon Whiteley's method, which breaks down
 *  the multivariable implicit solve into a sequence of 1d implicit solves
 */
class NhsSystemWithImplicitSolver : public NHSCellularMechanicsOdeSystem
{
private:
    const static double mTolerance = 1e-6; 
    double mDt;
    bool mUseImplicitExplicitSolveForZ;

    std::vector<double> mTempStoredStateVariables;
    std::vector<double> mCurrentStateVars;
    
    double mActiveTensionInitialGuess;
    double mActiveTensionSolution;

    void ImplicitSolveForActiveTension()
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
        
        mActiveTensionSolution = current_active_tension;
        mActiveTensionInitialGuess = mActiveTensionSolution;
    }

    double CalcActiveTensionResidual(double activeTensionGuess)
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

            
    
    double ImplicitSolveForCaTrop(double newActiveTension)
    {
        double old_Ca_trop = mCurrentStateVars[0];
        
        double numer = old_Ca_trop + mDt * mKon * mCalciumI  * mCalciumTroponinMax;
        double denom = 1 + mDt*mKon*mCalciumI + mDt*mKrefoff*(1-newActiveTension/(mGamma*mTref));

        return numer/denom;
    }


    double ImplicitSolveForZ(double newCaTrop)
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
    
    double CalcZResidual(double z, double newCaTrop)
    {
        double ca_trop_to_ca_trop50_ratio_to_n = pow(newCaTrop/mCalciumTrop50, mN);
        double dzdt =  mAlpha0 * ca_trop_to_ca_trop50_ratio_to_n * (1-z) 
                     - mAlphaR1 * z
                     - mAlphaR2 * pow(z,mNr) / (pow(z,mNr) + pow(mKZ,mNr));

        return z - mCurrentStateVars[1] - mDt*dzdt;                     
    }


    double ImplicitExplicitSolveForZ(double newCaTrop)
    {
        double ca_trop_to_ca_trop50_ratio_to_n = pow(newCaTrop/mCalciumTrop50, mN);
        double old_z = mCurrentStateVars[1];

        double numer =   old_z + mDt * mAlpha0 * ca_trop_to_ca_trop50_ratio_to_n 
                       - mDt * mAlphaR2 * pow(old_z,mNr) / (pow(old_z,mNr) + pow(mKZ,mNr));
        double denom =   1 + mDt * mAlpha0 * ca_trop_to_ca_trop50_ratio_to_n + mAlphaR1 *mDt;
        
        return numer/denom;
    }



    double ImplicitSolveForQ()
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

            

public :  
    NhsSystemWithImplicitSolver()
        : NHSCellularMechanicsOdeSystem()
    {
        mTempStoredStateVariables.resize(GetNumberOfStateVariables());
        mActiveTensionInitialGuess = 0.0;
        mUseImplicitExplicitSolveForZ = false;
    }
    
    void SetActiveTensionInitialGuess(double activeTensionInitialGuess)
    {
        mActiveTensionInitialGuess = activeTensionInitialGuess;
    }
    
    void SolveDoNotUpdate(double startTime, 
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
    
    void UpdateStateVariables()
    {
        for(unsigned i=0; i<mStateVariables.size(); i++)
        {
            mStateVariables[i] = mCurrentStateVars[i];
        }      
    }
    
    void UseImplicitExplicitSolveForZ(bool useImplicitExplicitSolveForZ = true)
    {
        mUseImplicitExplicitSolveForZ = useImplicitExplicitSolveForZ;
    }
};

#endif /*NHSSYSTEMWITHIMPLICITSOLVER_HPP_*/
