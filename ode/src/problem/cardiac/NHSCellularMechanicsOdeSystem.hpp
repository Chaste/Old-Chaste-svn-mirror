#ifndef NHSCELLULARMECHANICSODESYSTEM_HPP_
#define NHSCELLULARMECHANICSODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"

class NHSCellularMechanicsOdeSystem  : public AbstractOdeSystem
{
friend class TestCellularMechanicsOdeSystems;

private :
    double mLambda1;
    double mDLambda1Dt;
    double mCalciumI;
    
    // parameters
    /** what are these? */

    /** FILL IN. (microMols)^-1 sec^-1 */
    double mKon;

    /** FILL IN. sec^-1 */
    double mKrefoff;

    /** FILL IN. Dimensionless */
    double mGamma;

    /** FILL IN. microMols */
    double mCalciumTroponinMax;

    /** FILL IN. sec^-1 */
    double mAlphaR1;

    /** FILL IN. sec^-1 */
    double mAlphaR2;

    /** FILL IN. Dimensionless */
    double mKZ;

    /** FILL IN. Dimensionless */
    double mNr;

    /** FILL IN. Dimensionless */
    double mBeta1;

    /** FILL IN. sec^-1 */
    double mAlpha0;

    /** FILL IN. Dimensionless */
    double mN;

    /** FILL IN. Dimensionless */
    double mZp;

    /** FILL IN. microMols */
    double mCalcium50ref;

    /** FILL IN. kPa */
    double mTref;

    /** FILL IN. Dimensionless */
    double mBeta0;

    /** FILL IN. Dimensionless */
    double mA;  

    /** FILL IN. Dimensionless */
    double mA1;

    /** FILL IN. Dimensionless */
    double mA2;

    /** FILL IN. Dimensionless */
    double mA3;

    /** FILL IN. sec^-1 */
    double mAlpha1;

    /** FILL IN. sec^-1 */
    double mAlpha2;

    /** FILL IN. sec^-1 */
    double mAlpha3;
    
    double CalculateCalciumTop50()
    {
        double one_plus_beta1_times_lam_minus_one = 1 + mBeta1*(mLambda1-1);

        double calcium_trop_50 = mCalciumTroponinMax * mCalcium50ref * one_plus_beta1_times_lam_minus_one;
        calcium_trop_50 /= mCalcium50ref*one_plus_beta1_times_lam_minus_one + (1-one_plus_beta1_times_lam_minus_one/(2*mGamma))*mKrefoff/mKon;
        
        return calcium_trop_50;
    }


    double CalculateT0(double z)
    {
        double calcium_trop_50 = CalculateCalciumTop50();

        double zp_to_n_plus_K_to_n = pow(mZp,mNr) + pow(mKZ,mNr);
        
        double k1 = mAlphaR2 * pow(mZp,mNr-1) * mNr * pow(mKZ,mNr);
        k1 /= zp_to_n_plus_K_to_n * zp_to_n_plus_K_to_n;
        
        double k2 = mAlphaR2 * pow(mZp,mNr)/zp_to_n_plus_K_to_n;
        k2 *= 1 - mNr*pow(mKZ,mNr)/zp_to_n_plus_K_to_n;
        
        double calcium_ratio_to_n = pow(calcium_trop_50/mCalciumTroponinMax, mN);
        
        double z_max = mAlpha0 - k2*calcium_ratio_to_n;
        z_max /= mAlpha0 + (mAlphaR1 + k1)*calcium_ratio_to_n;

        double T0 = z * mTref * (1+mBeta0*(mLambda1-1)) / z_max;
        
        return T0;
    }
     
public : 
    NHSCellularMechanicsOdeSystem()
        :   AbstractOdeSystem(5) // five state variables
    {
        mVariableNames.push_back("CalciumTroponin");
        mVariableUnits.push_back("microMols");
        mStateVariables.push_back(0);
    
        mVariableNames.push_back("z");
        mVariableUnits.push_back("");
        mStateVariables.push_back(0);
    
        mVariableNames.push_back("Q1");
        mVariableUnits.push_back("");
        mStateVariables.push_back(0);

        mVariableNames.push_back("Q2");
        mVariableUnits.push_back("");
        mStateVariables.push_back(0);

        mVariableNames.push_back("Q3");
        mVariableUnits.push_back("");
        mStateVariables.push_back(0);
        
        mKon = 100;
        mKrefoff = 200;
        mGamma = 2;
        mCalciumTroponinMax = 70;
        mAlphaR1 = 2;
        mAlphaR2 = 1.75;
        mKZ = 0.15;
        mNr = 3;
        mBeta1 = -4;
        mAlpha0 = 8;
        mN = 3;
        mZp = 0.85;
        mCalcium50ref = 1.05;
        mTref = 56.2;
        mBeta0 = 4.9;
        mA = 0.35;  
        mA1 = -29;
        mA2 = 138;
        mA3 = 129;
        mAlpha1 = 30;
        mAlpha2 = 130;
        mAlpha3 = 625;
        
        mLambda1 = 1.0;
        mDLambda1Dt = 0.0;
        mCalciumI = 0.0;            
    }
    
    void SetLambda1DerivativeAndCalciumI(double lambda1, double dlambda1Dt, double calciumI)
    {
        mLambda1 = lambda1;
        mDLambda1Dt = dlambda1Dt;
        mCalciumI = calciumI;
    }
    

#define COVERAGE_IGNORE

    double GetCalciumTroponinValue()
    {
        return mStateVariables[0];
    }
    
#undef COVERAGE_IGNORE

    void EvaluateYDerivatives(double time,
                              const std::vector<double> &rY,
                              std::vector<double> &rDY)
    {
        const double& calcium_troponin = rY[0];
        const double& z = rY[1];
        const double& Q1 = rY[2];
        const double& Q2 = rY[3];
        const double& Q3 = rY[4];
        
        double Q = Q1 + Q2 + Q3;
        double T0 = CalculateT0(z);
        
        double Ta;
        if(Q>0)
        {
            Ta = T0*(1+(2+mA)*Q)/(1+Q);
        }
        else
        {
            Ta = T0*(1+mA*Q)/(1-Q);
        }

        rDY[0] =   mKon * mCalciumI * ( mCalciumTroponinMax - calcium_troponin)
                 - mKrefoff * (1-Ta/(mGamma*mTref)) * calcium_troponin;

        double calcium_trop_50 = CalculateCalciumTop50();
        double ca_trop_to_ca_trop50_ratio_to_n = pow(calcium_troponin/calcium_trop_50, mN);
                
        rDY[1] =   mAlpha0 * ca_trop_to_ca_trop50_ratio_to_n * (1-z) 
                 - mAlphaR1 * z
                  -mAlphaR2 * pow(z,mNr) / (pow(z,mNr) + pow(mKZ,mNr));
        
               
        rDY[2] = mAlpha1 * mDLambda1Dt - mAlpha1 * Q1;
        rDY[3] = mAlpha2 * mDLambda1Dt - mAlpha2 * Q2;
        rDY[4] = mAlpha3 * mDLambda1Dt - mAlpha3 * Q3;
    }
    
    
    double GetActiveTension()
    {
        double T0 = CalculateT0(mStateVariables[1]);
        double Q = mStateVariables[2]+mStateVariables[3]+mStateVariables[4];
        
        if(Q>0)
        {
            return T0*(1+(2+mA)*Q)/(1+Q);
        }
        else
        {
            return T0*(1+mA*Q)/(1-Q);
        }
    }
};
#endif /*NHSCELLULARMECHANICSODESYSTEM_HPP_*/
