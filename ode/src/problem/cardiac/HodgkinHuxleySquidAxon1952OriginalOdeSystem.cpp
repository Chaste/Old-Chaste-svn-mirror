#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"
#include "AbstractOdeSystem.hpp"
#include <cmath>
#include <cassert>

/**
 * Constructor
 */
HodgkinHuxleySquidAxon1952OriginalOdeSystem::HodgkinHuxleySquidAxon1952OriginalOdeSystem(AbstractStimulusFunction *stimulus):
AbstractOdeSystem()
{
   mpStimulus= stimulus;

   /*
    * Constants for the HodgkinHuxleySquidAxon1952OriginalOdeSystem model
    */

   leakage_current_g_L = 0.3;
   membrane_Cm = 1.0;
   membrane_E_R = 0.0;
   potassium_channel_g_K = 36.0;
   sodium_channel_g_Na = 120.0;

   /*
    * State variable
    */
   
    mVariableNames.push_back("V");
    mVariableUnits.push_back("mV");
    mInitialConditions.push_back(0.0);
    
    mVariableNames.push_back("n");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(0.31768);
    
    mVariableNames.push_back("h");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(0.59612);
    
    mVariableNames.push_back("m");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(0.05293);
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
 * @return std::vector<double> RHS of HodgkinHuxleySquidAxon1952OriginalOdeSystem system of equations
 */
std::vector<double> HodgkinHuxleySquidAxon1952OriginalOdeSystem::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
   /*
    * Typical initial conditions for the HodgkinHuxleySquidAxon1952OriginalOdeSystem model
    *
    * membrane_V = 0.0
    * potassium_channel_n_gate_n = 0.325
    * sodium_channel_h_gate_h = 0.6
    * sodium_channel_m_gate_m = 0.05
    */

   /*
    * Throw an exception if the initial vector is larger than the number of equations
    */

   assert(rY.size() == 4);

   double membrane_V = rY[0];
   double potassium_channel_n_gate_n = rY[1];
   double sodium_channel_h_gate_h = rY[2];
   double sodium_channel_m_gate_m = rY[3];

   /*
    * Compute the HodgkinHuxleySquidAxon1952OriginalOdeSystem model
    */

   double leakage_current_E_L = membrane_E_R-10.613;
   double leakage_current_i_L = leakage_current_g_L*(membrane_V-leakage_current_E_L);

   double membrane_i_Stim = mpStimulus->GetStimulus(time);
   
   double sodium_channel_E_Na = membrane_E_R-115.0;
   double sodium_channel_i_Na = sodium_channel_g_Na*pow(sodium_channel_m_gate_m, 3.0)*sodium_channel_h_gate_h*(membrane_V-sodium_channel_E_Na);
   double potassium_channel_E_K = membrane_E_R+12.0;
   double potassium_channel_i_K = potassium_channel_g_K*pow(potassium_channel_n_gate_n, 4.0)*(membrane_V-potassium_channel_E_K);
   double membrane_V_prime = -(-membrane_i_Stim+sodium_channel_i_Na+potassium_channel_i_K+leakage_current_i_L)/membrane_Cm;
   double potassium_channel_n_gate_alpha_n = 0.01*(membrane_V+10.0)/(exp((membrane_V+10.0)/10.0)-1.0);
   double potassium_channel_n_gate_beta_n = 0.125*exp(membrane_V/80.0);
   double potassium_channel_n_gate_n_prime = potassium_channel_n_gate_alpha_n*(1.0-potassium_channel_n_gate_n)-potassium_channel_n_gate_beta_n*potassium_channel_n_gate_n;
   double sodium_channel_h_gate_alpha_h = 0.07*exp(membrane_V/20.0);
   double sodium_channel_h_gate_beta_h = 1.0/(exp((membrane_V+30.0)/10.0)+1.0);
   double sodium_channel_h_gate_h_prime = sodium_channel_h_gate_alpha_h*(1.0-sodium_channel_h_gate_h)-sodium_channel_h_gate_beta_h*sodium_channel_h_gate_h;
   double sodium_channel_m_gate_alpha_m = 0.1*(membrane_V+25.0)/(exp((membrane_V+25.0)/10.0)-1.0);
   double sodium_channel_m_gate_beta_m = 4.0*exp(membrane_V/18.0);
   double sodium_channel_m_gate_m_prime = sodium_channel_m_gate_alpha_m*(1.0-sodium_channel_m_gate_m)-sodium_channel_m_gate_beta_m*sodium_channel_m_gate_m;

   std::vector<double> returnRHS;

   returnRHS.push_back(membrane_V_prime);
   returnRHS.push_back(potassium_channel_n_gate_n_prime);
   returnRHS.push_back(sodium_channel_h_gate_h_prime);
   returnRHS.push_back(sodium_channel_m_gate_m_prime);

   return returnRHS;
}
