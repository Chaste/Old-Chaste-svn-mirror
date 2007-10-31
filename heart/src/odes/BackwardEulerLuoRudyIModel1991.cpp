#include "BackwardEulerLuoRudyIModel1991.hpp"
#include "Exception.hpp"
#include "CardiacNewtonSolver.hpp"

#include <cmath>
#include <cassert>
//#include <iostream>


/**
 * Constructor
 */
BackwardEulerLuoRudyIModel1991::BackwardEulerLuoRudyIModel1991(
    double dt,
    AbstractStimulusFunction *pIntracellularStimulus,
    AbstractStimulusFunction *pExtracellularStimulus)
        : AbstractBackwardEulerCardiacCell<1>(8, 4, dt, pIntracellularStimulus,
                                              pExtracellularStimulus)
{
    Init();
}


/**
 * Constructor with the same signature as the forward cell models
 */
BackwardEulerLuoRudyIModel1991::BackwardEulerLuoRudyIModel1991(
   AbstractIvpOdeSolver *pSolver,
   double dt,
   AbstractStimulusFunction *pIntracellularStimulus,
   AbstractStimulusFunction *pExtracellularStimulus)
        : AbstractBackwardEulerCardiacCell<1>(8, 4, dt, pIntracellularStimulus,
                                              pExtracellularStimulus)
{
    Init();
}


/**
 * Destructor
 */
BackwardEulerLuoRudyIModel1991::~BackwardEulerLuoRudyIModel1991(void)
{
}


void BackwardEulerLuoRudyIModel1991::Init()
{
    // set final parameter member variable
    fast_sodium_current_E_Na = ((membrane_R * membrane_T) / membrane_F) *
                                log(ionic_concentrations_Nao / ionic_concentrations_Nai);

    // state variables
    mVariableNames.push_back("h");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(0.9804713);
    
    mVariableNames.push_back("j");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(0.98767124);
    
    mVariableNames.push_back("m");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(0.00187018);
    
    mVariableNames.push_back("CaI");
    mVariableUnits.push_back("mMol");
    mInitialConditions.push_back(0.0002);
    
    mVariableNames.push_back("V");
    mVariableUnits.push_back("mV");
    mInitialConditions.push_back(-83.853);
    
    mVariableNames.push_back("d");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(0.00316354);
    
    mVariableNames.push_back("f");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(0.99427859);
    
    mVariableNames.push_back("x");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(0.16647703);

    AbstractCardiacCell::Init();
}

double BackwardEulerLuoRudyIModel1991::GetIntracellularCalciumConcentration()
{
    return mStateVariables[3];
}


void BackwardEulerLuoRudyIModel1991::UpdateTransmembranePotential(double time)
{
    // Compute next value of V
    std::vector<double> &rY = rGetStateVariables();
    
    double fast_sodium_current_h_gate_h = rY[0];
    double fast_sodium_current_j_gate_j = rY[1];
    double fast_sodium_current_m_gate_m = rY[2];
    double intracellular_calcium_concentration_Cai = rY[3];
    double membrane_V = rY[4];
    double slow_inward_current_d_gate_d = rY[5];
    double slow_inward_current_f_gate_f = rY[6];
    double time_dependent_potassium_current_X_gate_X = rY[7];
    
    double background_current_i_b = background_current_g_b*(membrane_V-background_current_E_b);
    double fast_sodium_current_i_Na = fast_sodium_current_g_Na*pow(fast_sodium_current_m_gate_m, 3.0)*fast_sodium_current_h_gate_h*fast_sodium_current_j_gate_j*(membrane_V-fast_sodium_current_E_Na);
    double slow_inward_current_E_si = 7.7-13.0287*log(intracellular_calcium_concentration_Cai);
    double slow_inward_current_i_si = 0.09*slow_inward_current_d_gate_d*slow_inward_current_f_gate_f*(membrane_V-slow_inward_current_E_si);
    
    double time_dependent_potassium_current_g_K = 0.282*sqrt(ionic_concentrations_Ko/5.4);
    double time_dependent_potassium_current_Xi_gate_Xi;
    
    if (membrane_V > -100.0)
    {
        time_dependent_potassium_current_Xi_gate_Xi = 2.837*(exp(0.04*(membrane_V+77.0))-1.0)/((membrane_V+77.0)*exp(0.04*(membrane_V+35.0)));
    }
    else
    {
        #define COVERAGE_IGNORE
        time_dependent_potassium_current_Xi_gate_Xi = 1.0;
        #undef COVERAGE_IGNORE
    }
    
    double time_dependent_potassium_current_E_K = ((membrane_R*membrane_T)/membrane_F)*log((ionic_concentrations_Ko+time_dependent_potassium_current_PR_NaK*ionic_concentrations_Nao)/(ionic_concentrations_Ki+time_dependent_potassium_current_PR_NaK*ionic_concentrations_Nai));
    double time_dependent_potassium_current_i_K = time_dependent_potassium_current_g_K*time_dependent_potassium_current_X_gate_X*time_dependent_potassium_current_Xi_gate_Xi*(membrane_V-time_dependent_potassium_current_E_K);
    double time_independent_potassium_current_g_K1 = 0.6047*sqrt(ionic_concentrations_Ko/5.4);
    double time_independent_potassium_current_E_K1 =((membrane_R*membrane_T)/membrane_F)*log(ionic_concentrations_Ko/ionic_concentrations_Ki);
    double time_independent_potassium_current_K1_gate_alpha_K1 = 1.02/(1.0+exp(0.2385*(membrane_V-time_independent_potassium_current_E_K1-59.215)));
    double time_independent_potassium_current_K1_gate_beta_K1 = (0.49124*exp(0.08032*(membrane_V+5.476-time_independent_potassium_current_E_K1))+exp(0.06175*(membrane_V-(time_independent_potassium_current_E_K1+594.31))))/(1.0+exp(-0.5143*(membrane_V-time_independent_potassium_current_E_K1+4.753)));
    double time_independent_potassium_current_K1_gate_K1_infinity = time_independent_potassium_current_K1_gate_alpha_K1/(time_independent_potassium_current_K1_gate_alpha_K1+time_independent_potassium_current_K1_gate_beta_K1);
    double time_independent_potassium_current_i_K1 = time_independent_potassium_current_g_K1*time_independent_potassium_current_K1_gate_K1_infinity*(membrane_V-time_independent_potassium_current_E_K1);
    double plateau_potassium_current_Kp = 1.0/(1.0+exp((7.488-membrane_V)/5.98));
    double plateau_potassium_current_E_Kp = time_independent_potassium_current_E_K1;
    double plateau_potassium_current_i_Kp = plateau_potassium_current_g_Kp*plateau_potassium_current_Kp*(membrane_V-plateau_potassium_current_E_Kp);
    double i_stim = GetStimulus(time);
    
    //calculate dV
    double membrane_V_prime = (-1.0/membrane_C)*(fast_sodium_current_i_Na+slow_inward_current_i_si+time_dependent_potassium_current_i_K+time_independent_potassium_current_i_K1+plateau_potassium_current_i_Kp+background_current_i_b + i_stim);
    
    rY[4] += mDt * membrane_V_prime;
}

void BackwardEulerLuoRudyIModel1991::ComputeOneStepExceptVoltage(double tStart)
{
    // This method computes one timestep for all state variables except V, using
    // backward Euler.
    std::vector<double>& rY = rGetStateVariables();
    double fast_sodium_current_h_gate_h = rY[0];
    double fast_sodium_current_j_gate_j = rY[1];
    double fast_sodium_current_m_gate_m = rY[2];
    
    double membrane_V = rY[4];
    double slow_inward_current_d_gate_d = rY[5];
    double slow_inward_current_f_gate_f = rY[6];
    double time_dependent_potassium_current_X_gate_X = rY[7];
    
    double fast_sodium_current_h_gate_alpha_h;
    
    if (membrane_V < -40.0)
    {
        fast_sodium_current_h_gate_alpha_h = 0.135*exp((80.0+membrane_V)/-6.8);
    }
    else
    {
        fast_sodium_current_h_gate_alpha_h = 0.0;
    }
    
    double fast_sodium_current_h_gate_beta_h;
    
    if (membrane_V < -40.0)
    {
        fast_sodium_current_h_gate_beta_h = 3.56*exp(0.079*membrane_V)+3.1e5*exp(0.35*membrane_V);
    }
    else
    {
        fast_sodium_current_h_gate_beta_h = 1.0/(0.13*(1.0+exp((membrane_V+10.66)/-11.1)));
    }
    
    double fast_sodium_current_j_gate_alpha_j;
    
    if (membrane_V < -40.0)
    {
        fast_sodium_current_j_gate_alpha_j = (-1.2714e5*exp(0.2444*membrane_V)-3.474e-5*exp(-0.04391*membrane_V))*(membrane_V+37.78)/(1.0+exp(0.311*(membrane_V+79.23)));
    }
    else
    {
        fast_sodium_current_j_gate_alpha_j = 0.0;
    }
    
    double fast_sodium_current_j_gate_beta_j;
    
    if (membrane_V < -40.0)
    {
        fast_sodium_current_j_gate_beta_j = 0.1212*exp(-0.01052*membrane_V)/(1.0+exp(-0.1378*(membrane_V+40.14)));
    }
    else
    {
        fast_sodium_current_j_gate_beta_j = 0.3*exp(-2.535e-7*membrane_V)/(1.0+exp(-0.1*(membrane_V+32.0)));
    }
    
    double fast_sodium_current_m_gate_alpha_m = 0.32*(membrane_V+47.13)/(1.0-exp(-0.1*(membrane_V+47.13)));
    double fast_sodium_current_m_gate_beta_m = 0.08*exp(-membrane_V/11.0);
    
    double slow_inward_current_d_gate_alpha_d = 0.095*exp(-0.01*(membrane_V-5.0))/(1.0+exp(-0.072*(membrane_V-5.0)));
    double slow_inward_current_d_gate_beta_d = 0.07*exp(-0.017*(membrane_V+44.0))/(1.0+exp(0.05*(membrane_V+44.0)));
    
    double slow_inward_current_f_gate_alpha_f = 0.012*exp(-0.008*(membrane_V+28.0))/(1.0+exp(0.15*(membrane_V+28.0)));
    double slow_inward_current_f_gate_beta_f = 0.0065*exp(-0.02*(membrane_V+30.0))/(1.0+exp(-0.2*(membrane_V+30.0)));
    
    double time_dependent_potassium_current_X_gate_alpha_X = 0.0005*exp(0.083*(membrane_V+50.0))/(1.0+exp(0.057*(membrane_V+50.0)));
    double time_dependent_potassium_current_X_gate_beta_X = 0.0013*exp(-0.06*(membrane_V+20.0))/(1.0+exp(-0.04*(membrane_V+20.0)));
    
    // update h
    fast_sodium_current_h_gate_h =   (fast_sodium_current_h_gate_h + fast_sodium_current_h_gate_alpha_h*mDt)
                                     / (1 + (fast_sodium_current_h_gate_alpha_h + fast_sodium_current_h_gate_beta_h)*mDt);
                                     
    // update m
    fast_sodium_current_m_gate_m =   (fast_sodium_current_m_gate_m + fast_sodium_current_m_gate_alpha_m*mDt)
                                     / (1 + (fast_sodium_current_m_gate_alpha_m + fast_sodium_current_m_gate_beta_m)*mDt);
                                     
    // update j
    fast_sodium_current_j_gate_j =   (fast_sodium_current_j_gate_j + fast_sodium_current_j_gate_alpha_j*mDt)
                                     / (1 + (fast_sodium_current_j_gate_alpha_j + fast_sodium_current_j_gate_beta_j)*mDt);
                                     
    // update d
    slow_inward_current_d_gate_d =   (slow_inward_current_d_gate_d + slow_inward_current_d_gate_alpha_d*mDt)
                                     / (1 + (slow_inward_current_d_gate_alpha_d + slow_inward_current_d_gate_beta_d)*mDt);
                                     
    // update f
    slow_inward_current_f_gate_f =   (slow_inward_current_f_gate_f + slow_inward_current_f_gate_alpha_f*mDt)
                                     / (1 + (slow_inward_current_f_gate_alpha_f + slow_inward_current_f_gate_beta_f)*mDt);
                                     
    // update X
    time_dependent_potassium_current_X_gate_X =   ( time_dependent_potassium_current_X_gate_X + time_dependent_potassium_current_X_gate_alpha_X*mDt)
                                                  / ( 1 + (time_dependent_potassium_current_X_gate_alpha_X+time_dependent_potassium_current_X_gate_beta_X)*mDt);
                                                  
    rY[0] =  fast_sodium_current_h_gate_h;
    rY[1] =  fast_sodium_current_j_gate_j;
    rY[2] =  fast_sodium_current_m_gate_m;
    
    rY[5] =  slow_inward_current_d_gate_d;
    rY[6] =  slow_inward_current_f_gate_f;
    rY[7] =  time_dependent_potassium_current_X_gate_X;
    
    // finally, use the newton method to calculate the Calcium concentration at the next time
    double guess[1] = {rY[3]};
    
    CardiacNewtonSolver<1> *solver = CardiacNewtonSolver<1>::Instance();
    solver->Solve(*this, guess);
    
    rY[3] = guess[0];
    
    // Timings:
    //  Hardcoded one-variable newton:
    //    Forward: 1.66
    //    Backward: 1.29
    //    Backward (long dt): 0.03
    //  Original DoNewton:
    //    Forward: 1.69
    //    Backward: 1.28
    //    Backward (long dt): 0.02
    //  Generic DoNewton:
    //    Forward: 1.66, 1.67, 1.67
    //    Backward: 1.36, 1.31, 1.32
    //    Backward (long dt): 0.02, 0.03, 0.03
    //  Generic DoNewton with JonW norm:
    //    Forward: 1.67, 1.67, 1.67
    //    Backward: 1.30, 1.29, 1.30
    //    Backward (long dt): 0.03, 0.03, 0.02
    //  CardiacNewtonSolver with JonW norm:
    //    Forward: 1.70, 1.69, 1.70, 1.71
    //    Backward: 1.29, 1.29, 1.30, 1.32
    //    Backward (long dt): 0.03, 0.03, 0.02, 0.03
    // So the method call itself has no significant impact, nor the use of generic solver.
    // The choice of norm also seems immaterial for this problem.
}

void BackwardEulerLuoRudyIModel1991::ComputeResidual(const double rCurrentGuess[1], double rResidual[1])
{
    std::vector<double>& rY = rGetStateVariables();
    double membrane_V = rY[4];
    double slow_inward_current_d_gate_d = rY[5];
    double slow_inward_current_f_gate_f = rY[6];
    
    double slow_inward_current_E_si = 7.7-13.0287*log(rCurrentGuess[0]);
    double slow_inward_current_i_si = 0.09*slow_inward_current_d_gate_d*slow_inward_current_f_gate_f*(membrane_V-slow_inward_current_E_si);
    
    double dCdt_using_current_guess = -1e-4*slow_inward_current_i_si+0.07*(1e-4-rCurrentGuess[0]);
    
    rResidual[0] = rCurrentGuess[0] - rY[3] - mDt*dCdt_using_current_guess;
}

void BackwardEulerLuoRudyIModel1991::ComputeJacobian(const double rCurrentGuess[1], double rJacobian[1][1])
{
    std::vector<double>& rY = rGetStateVariables();
    double slow_inward_current_d_gate_d = rY[5];
    double slow_inward_current_f_gate_f = rY[6];
    
    rJacobian[0][0] = 1 - mDt*(-0.1172583e-3*slow_inward_current_d_gate_d*slow_inward_current_f_gate_f/rCurrentGuess[0]
                               - 0.07);
}

double BackwardEulerLuoRudyIModel1991::GetIIonic()
{
    double fast_sodium_current_h_gate_h = mStateVariables[0];
    double fast_sodium_current_j_gate_j = mStateVariables[1];
    double fast_sodium_current_m_gate_m = mStateVariables[2];
    double intracellular_calcium_concentration_Cai = mStateVariables[3];
    double membrane_V = mStateVariables[4];
    double slow_inward_current_d_gate_d = mStateVariables[5];
    double slow_inward_current_f_gate_f = mStateVariables[6];
    double time_dependent_potassium_current_X_gate_X = mStateVariables[7];
    
    /*
     * Compute the LuoRudyIModel1991OdeSystem model
     */
    double background_current_i_b = background_current_g_b*(membrane_V-background_current_E_b);
    
    double fast_sodium_current_i_Na = fast_sodium_current_g_Na*pow(fast_sodium_current_m_gate_m, 3.0)*fast_sodium_current_h_gate_h*fast_sodium_current_j_gate_j*(membrane_V-fast_sodium_current_E_Na);
    double slow_inward_current_E_si = 7.7-13.0287*log(intracellular_calcium_concentration_Cai);

    double slow_inward_current_i_si = 0.09*slow_inward_current_d_gate_d*slow_inward_current_f_gate_f*(membrane_V-slow_inward_current_E_si);
    
    double time_dependent_potassium_current_g_K = 0.282*sqrt(ionic_concentrations_Ko/5.4);
    double time_dependent_potassium_current_Xi_gate_Xi;
    
    // Although the equation below looks strange (particularly the arguments of the
    // exponentials, it is in fact correct.
    if (membrane_V > -100.0)
    {
        time_dependent_potassium_current_Xi_gate_Xi = 2.837*(exp(0.04*(membrane_V+77.0))-1.0)/((membrane_V+77.0)*exp(0.04*(membrane_V+35.0)));
    }
    else
    {
        #define COVERAGE_IGNORE
        time_dependent_potassium_current_Xi_gate_Xi = 1.0;
        #undef COVERAGE_IGNORE
    }

    double time_dependent_potassium_current_E_K = ((membrane_R*membrane_T)/membrane_F)*log((ionic_concentrations_Ko+time_dependent_potassium_current_PR_NaK*ionic_concentrations_Nao)/(ionic_concentrations_Ki+time_dependent_potassium_current_PR_NaK*ionic_concentrations_Nai));
    double time_dependent_potassium_current_i_K = time_dependent_potassium_current_g_K*time_dependent_potassium_current_X_gate_X*time_dependent_potassium_current_Xi_gate_Xi*(membrane_V-time_dependent_potassium_current_E_K);
    
    double time_independent_potassium_current_g_K1 = 0.6047*sqrt(ionic_concentrations_Ko/5.4);
    double time_independent_potassium_current_E_K1 =((membrane_R*membrane_T)/membrane_F)*log(ionic_concentrations_Ko/ionic_concentrations_Ki);
    double time_independent_potassium_current_K1_gate_alpha_K1 = 1.02/(1.0+exp(0.2385*(membrane_V-time_independent_potassium_current_E_K1-59.215)));
    double time_independent_potassium_current_K1_gate_beta_K1 = (0.49124*exp(0.08032*(membrane_V+5.476-time_independent_potassium_current_E_K1))+exp(0.06175*(membrane_V-(time_independent_potassium_current_E_K1+594.31))))/(1.0+exp(-0.5143*(membrane_V-time_independent_potassium_current_E_K1+4.753)));
    double time_independent_potassium_current_K1_gate_K1_infinity = time_independent_potassium_current_K1_gate_alpha_K1/(time_independent_potassium_current_K1_gate_alpha_K1+time_independent_potassium_current_K1_gate_beta_K1);
    double time_independent_potassium_current_i_K1 = time_independent_potassium_current_g_K1*time_independent_potassium_current_K1_gate_K1_infinity*(membrane_V-time_independent_potassium_current_E_K1);
    
    double plateau_potassium_current_Kp = 1.0/(1.0+exp((7.488-membrane_V)/5.98));
    double plateau_potassium_current_E_Kp = time_independent_potassium_current_E_K1;
    double plateau_potassium_current_i_Kp = plateau_potassium_current_g_Kp*plateau_potassium_current_Kp*(membrane_V-plateau_potassium_current_E_Kp);
    
    double i_ionic = fast_sodium_current_i_Na+slow_inward_current_i_si+time_dependent_potassium_current_i_K+time_independent_potassium_current_i_K1+plateau_potassium_current_i_Kp+background_current_i_b;

    assert( !isnan(i_ionic));
    return i_ionic;
}


void BackwardEulerLuoRudyIModel1991::VerifyStateVariables()
{
//#ifndef NDEBUG
    const std::vector<double>& rY = rGetStateVariables();
 
    const double fast_sodium_current_h_gate_h = rY[0];            // gating
    const double fast_sodium_current_j_gate_j = rY[1];            // gating
    const double fast_sodium_current_m_gate_m = rY[2];            // gating
    const double intracellular_calcium_concentration_Cai = rY[3]; // concentration
    const double slow_inward_current_d_gate_d = rY[5];            // gating
    const double slow_inward_current_f_gate_f = rY[6];            // gating
    const double time_dependent_potassium_current_X_gate_X = rY[7]; // gating
    
    #define COVERAGE_IGNORE
    if (!(0.0<=fast_sodium_current_h_gate_h && fast_sodium_current_h_gate_h<=1.0))
    {
        EXCEPTION("h gate for fast sodium current has gone out of range. Check model parameters, for example spatial stepsize");
    }
    
    if (!(0.0<=fast_sodium_current_j_gate_j && fast_sodium_current_j_gate_j<=1.0))
    {
        EXCEPTION("j gate for fast sodium current has gone out of range. Check model parameters, for example spatial stepsize");
    }
    
    if (!(0.0<=fast_sodium_current_m_gate_m && fast_sodium_current_m_gate_m<=1.0))
    {
        EXCEPTION("m gate for fast sodium current has gone out of range. Check model parameters, for example spatial stepsize");
    }

    if (!(0.0<intracellular_calcium_concentration_Cai))
    {
        EXCEPTION("intracellular_calcium_concentration_Cai has become non-positive, ie gone out of range. Check model parameters, for example spatial stepsize");
    }

    if (!(0.0<=slow_inward_current_d_gate_d && slow_inward_current_d_gate_d<=1.0))
    {
        EXCEPTION("d gate for slow inward current has gone out of range. Check model parameters, for example spatial stepsize");
    }
    
    if (!(0.0<=slow_inward_current_f_gate_f && slow_inward_current_f_gate_f<=1.0))
    {
        EXCEPTION("f gate for slow inward current has gone out of range. Check model parameters, for example spatial stepsize");
    }
    
    if (!(0.0<=time_dependent_potassium_current_X_gate_X && time_dependent_potassium_current_X_gate_X<=1.0))
    {
        EXCEPTION("X gate for time dependent potassium current has gone out of range. Check model parameters, for example spatial stepsize");
    }
    #undef COVERAGE_IGNORE
//#endif
}

