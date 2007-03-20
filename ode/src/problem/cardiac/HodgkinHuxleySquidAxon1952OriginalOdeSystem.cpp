#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"
#include "AbstractOdeSystem.hpp"
#include <cmath>


/**
 * Constructor
 */
HodgkinHuxleySquidAxon1952OriginalOdeSystem::HodgkinHuxleySquidAxon1952OriginalOdeSystem(
    AbstractIvpOdeSolver *pOdeSolver,
    double dt,
    AbstractStimulusFunction *pIntracellularStimulus,
    AbstractStimulusFunction *pExtracellularStimulus)
        : AbstractCardiacCell(pOdeSolver,4,0,dt,pIntracellularStimulus,pExtracellularStimulus)
{
    /*
     * State variable
     */
    mVariableNames.push_back("V");
    mVariableUnits.push_back("mV");
    mInitialConditions.push_back(-75.0);
    
    mVariableNames.push_back("n");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(0.325);
    
    mVariableNames.push_back("h");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(0.6);
    
    mVariableNames.push_back("m");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(0.05);
    
    Init();
}

/**
 * Destructor
 */
HodgkinHuxleySquidAxon1952OriginalOdeSystem::~HodgkinHuxleySquidAxon1952OriginalOdeSystem(void)
{
    // Do nothing
}

/**
 * Function returns a vector representing the RHS of the HodgkinHuxleySquidAxon1952OriginalOdeSystem system of Odes at each time step, y' = [y1' ... yn'].
 * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
 *
 * @param rDY filled in derivatives using HodgkinHuxleySquidAxon1952OriginalOdeSystem system of equations
 */
void HodgkinHuxleySquidAxon1952OriginalOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
{
    /*
     * Typical initial conditions for the HodgkinHuxleySquidAxon1952OriginalOdeSystem model
     *
     * membrane_V = 0.0
     * potassium_channel_n_gate_n = 0.325
     * sodium_channel_h_gate_h = 0.6
     * sodium_channel_m_gate_m = 0.05
     */
    
    double membrane_V = rY[0];
    double potassium_channel_n_gate_n = rY[1];
    double sodium_channel_h_gate_h = rY[2];
    double sodium_channel_m_gate_m = rY[3];
    
    /*
     * Compute the HodgkinHuxleySquidAxon1952OriginalOdeSystem model
     */
    
    double leakage_current_E_L = membrane_E_R+10.613;
    double leakage_current_i_L = leakage_current_g_L*(membrane_V-leakage_current_E_L);
    
    double membrane_i_Stim = GetStimulus(time);
    
    double sodium_channel_E_Na = membrane_E_R+115.0;
    double sodium_channel_i_Na = sodium_channel_g_Na*sodium_channel_m_gate_m*sodium_channel_m_gate_m*sodium_channel_m_gate_m*sodium_channel_h_gate_h*(membrane_V-sodium_channel_E_Na);
    double potassium_channel_E_K = membrane_E_R-12.0;
    double potassium_channel_i_K = potassium_channel_g_K*potassium_channel_n_gate_n*potassium_channel_n_gate_n*potassium_channel_n_gate_n*potassium_channel_n_gate_n*(membrane_V-potassium_channel_E_K);
    
    double membrane_V_prime = -(-membrane_i_Stim+sodium_channel_i_Na+potassium_channel_i_K+leakage_current_i_L)/membrane_Cm;
    // do not update voltage if the mSetVoltageDerivativeToZero flag has been set
    if (mSetVoltageDerivativeToZero)
    {
        membrane_V_prime = 0;
    }
    
    double potassium_channel_n_gate_alpha_n;
    if (-65.0001<membrane_V && membrane_V<-64.9999)
    {
        potassium_channel_n_gate_alpha_n = 0.1;
    }
    else
    {
        potassium_channel_n_gate_alpha_n = -0.01*(membrane_V+65.0)/(exp(-(membrane_V+65.0)/10.0)-1.0);
    }
    
    double potassium_channel_n_gate_beta_n = 0.125*exp((membrane_V+75.0)/80.0);
    double potassium_channel_n_gate_n_prime = potassium_channel_n_gate_alpha_n*(1.0-potassium_channel_n_gate_n)-potassium_channel_n_gate_beta_n*potassium_channel_n_gate_n;
    double sodium_channel_h_gate_alpha_h = 0.07*exp(-(membrane_V+75.0)/20.0);
    double sodium_channel_h_gate_beta_h = 1.0/(exp(-(membrane_V+45.0)/10.0)+1.0);
    double sodium_channel_h_gate_h_prime = sodium_channel_h_gate_alpha_h*(1.0-sodium_channel_h_gate_h)-sodium_channel_h_gate_beta_h*sodium_channel_h_gate_h;
    
    double sodium_channel_m_gate_alpha_m;
    if (-50.0001<membrane_V && membrane_V<-49.9999)
    {
        sodium_channel_m_gate_alpha_m = 1;
    }
    else
    {
        sodium_channel_m_gate_alpha_m = -0.1*(membrane_V+50.0)/(exp(-(membrane_V+50.0)/10.0)-1.0);
    }
    double sodium_channel_m_gate_beta_m = 4.0*exp(-(membrane_V+75.0)/18.0);
    double sodium_channel_m_gate_m_prime = sodium_channel_m_gate_alpha_m*(1.0-sodium_channel_m_gate_m)-sodium_channel_m_gate_beta_m*sodium_channel_m_gate_m;
    
    rDY[0] = membrane_V_prime;
    rDY[1] = potassium_channel_n_gate_n_prime;
    rDY[2] = sodium_channel_h_gate_h_prime;
    rDY[3] = sodium_channel_m_gate_m_prime;
}

void HodgkinHuxleySquidAxon1952OriginalOdeSystem::Init()
{
    AbstractCardiacCell::Init();
    leakage_current_g_L = 0.3;   // mS/cm2
    membrane_Cm = 1.0;   // 1 uF/cm2
    membrane_E_R = -75.0;   // mV
    potassium_channel_g_K = 36.0;   // mS/cm2
    sodium_channel_g_Na = 120.0;   // mS/cm2
}

double HodgkinHuxleySquidAxon1952OriginalOdeSystem::GetIIonic()
{
    double membrane_V = mStateVariables[mVoltageIndex];
    double potassium_channel_n_gate_n = mStateVariables[1];
    double sodium_channel_h_gate_h = mStateVariables[2];
    double sodium_channel_m_gate_m = mStateVariables[3];
    
    /*
     * Compute the HodgkinHuxleySquidAxon1952OriginalOdeSystem model
     */
    
    double leakage_current_E_L = membrane_E_R+10.613;
    double leakage_current_i_L = leakage_current_g_L*(membrane_V-leakage_current_E_L);
    
    double sodium_channel_E_Na = membrane_E_R+115.0;
    double sodium_channel_i_Na = sodium_channel_g_Na*sodium_channel_m_gate_m*sodium_channel_m_gate_m*sodium_channel_m_gate_m*sodium_channel_h_gate_h*(membrane_V-sodium_channel_E_Na);
    double potassium_channel_E_K = membrane_E_R-12.0;
    double potassium_channel_i_K = potassium_channel_g_K*potassium_channel_n_gate_n*potassium_channel_n_gate_n*potassium_channel_n_gate_n*potassium_channel_n_gate_n*(membrane_V-potassium_channel_E_K);
    double i_ionic = sodium_channel_i_Na+potassium_channel_i_K+leakage_current_i_L;
    return i_ionic;
    
}
