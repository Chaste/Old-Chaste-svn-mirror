
#include "NHSCellularMechanicsOdeSystem.hpp"
#include <cmath>

/*
 * ============================== PRIVATE FUNCTIONS =====================================
 */
void NHSCellularMechanicsOdeSystem::CalculateCalciumTrop50()
{
    double one_plus_beta1_times_lam_minus_one = 1 + mBeta1*(mLambda1-1);

    mCalciumTrop50 = mCalciumTroponinMax * mCalcium50ref * one_plus_beta1_times_lam_minus_one;
    mCalciumTrop50 /= mCalcium50ref*one_plus_beta1_times_lam_minus_one + (1-one_plus_beta1_times_lam_minus_one/(2*mGamma))*mKrefoff/mKon;
}


double NHSCellularMechanicsOdeSystem::CalculateT0(double z)
{
    double calcium_ratio_to_n = pow(mCalciumTrop50/mCalciumTroponinMax, mN);
    
    double z_max = mAlpha0 - mK2*calcium_ratio_to_n;
    z_max /= mAlpha0 + (mAlphaR1 + mK1)*calcium_ratio_to_n;

    return z * mTref * (1+mBeta0*(mLambda1-1)) / z_max;
}



/*
 * ============================== PRIVATE FUNCTIONS =====================================
 */

NHSCellularMechanicsOdeSystem::NHSCellularMechanicsOdeSystem()
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
            
    mLambda1 = 1.0;
    mDLambda1Dt = 0.0;
    mCalciumI = 0.0;            
    
    // Initialise mCalciumTrop50!!
    CalculateCalciumTrop50();
    
    double zp_to_n_plus_K_to_n = pow(mZp,mNr) + pow(mKZ,mNr);
    
    mK1 = mAlphaR2 * pow(mZp,mNr-1) * mNr * pow(mKZ,mNr);
    mK1 /= zp_to_n_plus_K_to_n * zp_to_n_plus_K_to_n;
    
    mK2 = mAlphaR2 * pow(mZp,mNr)/zp_to_n_plus_K_to_n;
    mK2 *= 1 - mNr*pow(mKZ,mNr)/zp_to_n_plus_K_to_n;
}

void NHSCellularMechanicsOdeSystem::SetLambda1AndDerivative(double lambda1, double dlambda1Dt)
{
    assert(lambda1>0.0);
    mLambda1 = lambda1;
    mDLambda1Dt = dlambda1Dt;
    // lambda changed so update mCalciumTrop50!!
    CalculateCalciumTrop50();
}

void NHSCellularMechanicsOdeSystem::SetIntracellularCalciumConcentration(double calciumI)
{
    assert(calciumI > 0.0);
    mCalciumI = calciumI;
}


double NHSCellularMechanicsOdeSystem::GetCalciumTroponinValue()
{
    return mStateVariables[0];
}

void NHSCellularMechanicsOdeSystem::EvaluateYDerivatives(double time,
                                                         const std::vector<double> &rY,
                                                         std::vector<double> &rDY)
{
    const double& calcium_troponin = rY[0];
    const double& z = rY[1];
    const double& Q1 = rY[2];
    const double& Q2 = rY[3];
    const double& Q3 = rY[4];
    
    // check the state vars are in the expected range
    #define COVERAGE_IGNORE
    if(calcium_troponin < 0)
    {
        EXCEPTION("CalciumTrop concentration went negative");
    }
    if(z<0)
    {
        EXCEPTION("z went negative");
    }
    if(z>1)
    {
        EXCEPTION("z became greater than 1");
    }
    #undef COVERAGE_IGNORE

            
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

    double ca_trop_to_ca_trop50_ratio_to_n = pow(calcium_troponin/mCalciumTrop50, mN);
            
    rDY[1] =   mAlpha0 * ca_trop_to_ca_trop50_ratio_to_n * (1-z) 
             - mAlphaR1 * z
             - mAlphaR2 * pow(z,mNr) / (pow(z,mNr) + pow(mKZ,mNr));
    
           
    rDY[2] = mA1 * mDLambda1Dt - mAlpha1 * Q1;
    rDY[3] = mA2 * mDLambda1Dt - mAlpha2 * Q2;
    rDY[4] = mA3 * mDLambda1Dt - mAlpha3 * Q3;
}


double NHSCellularMechanicsOdeSystem::GetActiveTension()
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
