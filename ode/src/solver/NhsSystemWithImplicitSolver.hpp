#ifndef NHSSYSTEMWITHIMPLICITSOLVER_HPP_
#define NHSSYSTEMWITHIMPLICITSOLVER_HPP_

#include "NHSCellularMechanicsOdeSystem.hpp"
#include "EulerIvpOdeSolver.hpp"

class NhsSystemWithImplicitSolver : public NHSCellularMechanicsOdeSystem
{
private:
    const static double mTolerance = 1e-6; 
    double mDt;

    std::vector<double> mStoredStateVariables;
    double mActiveTensionInitialGuess;
    double mActiveTensionSolution;

    void ImplicitSolveForActiveTension()
    {
        double current_active_tension = mActiveTensionInitialGuess;
        
        double residual = CalcActiveTensionResidual(current_active_tension);
        std::cout << "ImplicitSolveForActiveTension: " << residual << "\n";

        unsigned counter = 0;
        while ((fabs(residual)>mTolerance) && (counter++<15))
        {
            double h = current_active_tension/100;
            double jac = (CalcActiveTensionResidual(current_active_tension+h) - CalcActiveTensionResidual(current_active_tension))/h;
            
            current_active_tension -= jac/residual;
            residual = CalcActiveTensionResidual(current_active_tension);
            std::cout << "ImplicitSolveForActiveTension: " << residual << "\n";
        } 
        assert(counter<15);
        
        mActiveTensionSolution = current_active_tension;
        mActiveTensionInitialGuess = mActiveTensionSolution;
    }

    double CalcActiveTensionResidual(double activeTensionGuess)
    {
        // solve for new Ca_trop implicitly
        double new_Ca_Trop = ImplicitSolveForCaTrop(activeTensionGuess);
        
        mStoredStateVariables[0] = new_Ca_Trop;
        
        // solve for new Ca_trop implicitly
        double new_z = ImplicitSolveForZ(new_Ca_Trop);

        mStoredStateVariables[1] = new_z;

        // get the new T0 
        double new_T0 = CalculateT0(new_z); std::cout << "new T0 = " << new_T0 <<"\n";
    
        // solve for the new Qi's (and therefore Q) implicitly
        double new_Q = ImplicitSolveForQ(); std::cout << "new Q = " << new_Q <<"\n";  // can be moved up
                
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
        double old_Ca_trop = mStateVariables[0];
        
        double numer = old_Ca_trop + mDt * mKon * mCalciumI  * mCalciumTroponinMax;
        double denom = 1 + mDt*mKon*mCalciumI + mDt*mKrefoff*(1-newActiveTension/(mGamma*mTref));

        std::cout << "ImplicitSolveForCaTrop: ca=" << numer/denom  << "\n";
        
        return numer/denom;
    }


    double ImplicitSolveForZ(double newCaTrop)
    {
        double current_z = mStateVariables[1];
        double residual = CalcZResidual(current_z, newCaTrop);
            std::cout << "ImplicitSolveForZ: z,resid=" <<current_z << ", " << residual << "\n";
        unsigned counter = 0;
        while ((fabs(residual)>mTolerance) && (counter++<15))
        {
            
            double h = current_z/100;
            double jac = (CalcZResidual(current_z + h, newCaTrop) - CalcZResidual(current_z, newCaTrop))/h;
            
            current_z -= jac/residual;
            residual = CalcZResidual(current_z, newCaTrop);
            std::cout << "ImplicitSolveForZ: resid=" <<residual << "\n";
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

        return z - mStateVariables[1] - mDt*dzdt;                     
    }

//    double ImplicitExplicitSolveForZ(double newCaTrop)
//    {
//        // fill in -    
//    }



    double ImplicitSolveForQ()
    {
        double old_Q1 = mStateVariables[2];
        double old_Q2 = mStateVariables[3];
        double old_Q3 = mStateVariables[4];
        
        double new_Q1 = (old_Q1 + mDt*mA1*mDLambda1Dt)/(1 + mAlpha1*mDt);
        double new_Q2 = (old_Q2 + mDt*mA2*mDLambda1Dt)/(1 + mAlpha2*mDt);
        double new_Q3 = (old_Q3 + mDt*mA3*mDLambda1Dt)/(1 + mAlpha3*mDt);
        
        mStoredStateVariables[2] = new_Q1;
        mStoredStateVariables[3] = new_Q2;
        mStoredStateVariables[4] = new_Q3;
        
        std::cout << "ImplicitSolveForQ: Q=" << new_Q1 + new_Q2 + new_Q3  << "\n";
        
        return new_Q1 + new_Q2 + new_Q3; 
    }

            

public :  
    NhsSystemWithImplicitSolver()
        : NHSCellularMechanicsOdeSystem()
    {
        mStoredStateVariables.resize(GetNumberOfStateVariables());
        mActiveTensionInitialGuess = 0.0;
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
        
        // store the initial state vars
        std::vector<double> old_state_vars = rGetStateVariables();
        
//todo: proper implicit solve
        // solve
//        EulerIvpOdeSolver solver;
  //      solver.SolveAndUpdateStateVariable(this, startTime, endTime, timestep);

        double time = startTime;
        while(time < endTime)
        {
            std::cout << " ** TIME = " << time << "\n";
            ImplicitSolveForActiveTension();
            time += timestep;
        }


        
        // rewrite with the old state vars
        for(unsigned i=0; i<rGetStateVariables().size(); i++)
        {
            rGetStateVariables()[i] = old_state_vars[i];
        }        
    }
    
    void UpdateSystem()
    {
        for(unsigned i=0; i<mStateVariables.size(); i++)
        {
            mStateVariables[i] = mStoredStateVariables[i];
        }      
    }
};

#endif /*NHSSYSTEMWITHIMPLICITSOLVER_HPP_*/
