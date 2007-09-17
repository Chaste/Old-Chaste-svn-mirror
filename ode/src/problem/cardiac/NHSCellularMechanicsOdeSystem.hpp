#ifndef NHSCELLULARMECHANICSODESYSTEM_HPP_
#define NHSCELLULARMECHANICSODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"

class NHSCellularMechanicsOdeSystem  : public AbstractOdeSystem
{
friend class TestCellularMechanicsOdeSystems;

protected :
    double mLambda1;
    double mDLambda1Dt;
    double mCalciumI;
    double mCalciumTrop50; // only dependent on constants and lambda, so updated whenever lambda is updated

    double mK1;
    double mK2;

    /** FILL IN. (mMols)^-1 (ms)^-1 */
    static const double mKon = 100;

    /** FILL IN. (ms)^-1 */
    static const double mKrefoff = 0.2;

    /** FILL IN. Dimensionless */
    static const double mGamma = 2;

    /** FILL IN. mMols */
    static const double mCalciumTroponinMax = 0.07;

    /** FILL IN. (ms)^-1 */
    static const double mAlphaR1 = 0.002;

    /** FILL IN. (ms)^-1 */
    static const double mAlphaR2 = 0.00175;

    /** FILL IN. Dimensionless */
    static const double mKZ = 0.15;

    /** FILL IN. Dimensionless */
    static const double mNr = 3;

    /** FILL IN. Dimensionless */
    static const double mBeta1 = -4;

    /** FILL IN. (ms)^-1 */
    static const double mAlpha0 = 0.008;

    /** FILL IN. Dimensionless */
    static const double mN = 3;

    /** FILL IN. Dimensionless */
    static const double mZp = 0.85;

    /** FILL IN. mMols */
    static const double mCalcium50ref = 0.00105;

    /** FILL IN. kPa */
    static const double mTref = 56.2;

    /** FILL IN. Dimensionless */
    static const double mBeta0 = 4.9;

    /** FILL IN. Dimensionless */
    static const double mA = 0.35;  

    /** FILL IN. Dimensionless */
    static const double mA1 = -29;

    /** FILL IN. Dimensionless */
    static const double mA2 = 138;

    /** FILL IN. Dimensionless */
    static const double mA3 = 129;

    /** FILL IN. (ms)^-1 */
    static const double mAlpha1 = 0.03;

    /** FILL IN. (ms)^-1 */
    static const double mAlpha2 = 0.130;

    /** FILL IN. (ms)^-1 */
    static const double mAlpha3 = 0.625;
    

    void CalculateCalciumTrop50()
    {
        double one_plus_beta1_times_lam_minus_one = 1 + mBeta1*(mLambda1-1);

        mCalciumTrop50 = mCalciumTroponinMax * mCalcium50ref * one_plus_beta1_times_lam_minus_one;
        mCalciumTrop50 /= mCalcium50ref*one_plus_beta1_times_lam_minus_one + (1-one_plus_beta1_times_lam_minus_one/(2*mGamma))*mKrefoff/mKon;
    }


    double CalculateT0(double z)
    {
        double calcium_ratio_to_n = pow(mCalciumTrop50/mCalciumTroponinMax, mN);
        
        double z_max = mAlpha0 - mK2*calcium_ratio_to_n;
        z_max /= mAlpha0 + (mAlphaR1 + mK1)*calcium_ratio_to_n;

        return z * mTref * (1+mBeta0*(mLambda1-1)) / z_max;
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
    
    void SetLambda1AndDerivative(double lambda1, double dlambda1Dt)
    {
        assert(lambda1>0.0);
        mLambda1 = lambda1;
        mDLambda1Dt = dlambda1Dt;
        // lambda changed so update mCalciumTrop50!!
        CalculateCalciumTrop50();
    }
    
    void SetIntracellularCalciumConcentration(double calciumI)
    {
        assert(calciumI > 0.0);
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
