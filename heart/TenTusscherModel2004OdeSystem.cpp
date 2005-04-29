#include "TenTusscherModel2004OdeSystem.hpp"
#include "AbstractOdeSystem.hpp"
#include <cmath>
#include <cassert>

/**
 * Constructor
 */
TenTusscherModel2004OdeSystem::TenTusscherModel2004OdeSystem(AbstractStimulusFunction *stimulus):
AbstractOdeSystem(17)
{
   mpStimulus= stimulus;

   /*
    * Constants for the TenTusscherModel2004OdeSystem model
    */

   calcium_background_current_g_bCa = 0.000592;
   calcium_dynamics_a_rel = 16.464;
   calcium_dynamics_b_rel = 0.25;
   calcium_dynamics_Bufc = 0.15;
   calcium_dynamics_Bufsr = 10.0;
   calcium_dynamics_c_rel = 8.232;
   calcium_dynamics_K_up = 0.00025;
   calcium_dynamics_Kbufc = 0.001;
   calcium_dynamics_Kbufsr = 0.3;
   calcium_dynamics_tau_g = 2.0;    // Couldn't find in paper
   calcium_dynamics_V_leak = 8E-5;
   calcium_dynamics_Vmax_up = 0.000425;
   calcium_dynamics_Vsr = 1094.0;
   calcium_pump_current_g_pCa = 0.025;
   calcium_pump_current_K_pCa = 0.0005;
   fast_sodium_current_g_Na = 14.838;
   inward_rectifier_potassium_current_g_K1 = 5.405;
   L_type_calcium_current_f_Ca_gate_tau_f_Ca = 2.0; // Couldn't find in paper
   L_type_calcium_current_g_CaL = 1.75E-4;
   membrane_Cm = 2.0;
   model_parameters_Ca_o = 2.0;
   model_parameters_F = 96.4867; //coulombs per mill. Mol.
   model_parameters_K_o = 5.4;
   model_parameters_Na_o = 140.0;
   model_parameters_R = 8.3143;
   model_parameters_T = 310.0;
   model_parameters_Vc = 16404.0;
   Na_Ca_exchanger_current_alpha = 2.5;
   Na_Ca_exchanger_current_gamma = 0.35;
   Na_Ca_exchanger_current_K_mCa = 1.38;
   Na_Ca_exchanger_current_K_mNa_i = 87.5;
   Na_Ca_exchanger_current_K_NaCa = 1000.0;
   Na_Ca_exchanger_current_K_sat = 0.1;
   potassium_pump_current_g_pK = 0.0146;
   rapid_delayed_rectifier_current_g_Kr = 0.096;
   reversal_potentials_p_KNa = 0.03;
   slow_delayed_rectifier_current_g_Ks = 0.245;
   sodium_background_current_g_bNa = 0.00029;
   sodium_potassium_pump_current_K_mK = 1.0;
   sodium_potassium_pump_current_K_mNa = 40.0;
   sodium_potassium_pump_current_P_NaK = 1.362;
   transient_outward_current_g_to = 0.294;
}

/**
 * Destructor
 */
TenTusscherModel2004OdeSystem::~TenTusscherModel2004OdeSystem(void)
{
   // Do nothing
}

/**
 * Function returns a vector representing the RHS of the TenTusscherModel2004OdeSystem system of Odes at each time step, y' = [y1' ... yn'].
 * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
 *
 * @return std::vector<double> RHS of TenTusscherModel2004OdeSystem system of equations
 */
std::vector<double> TenTusscherModel2004OdeSystem::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
   /*
    * Typical initial conditions for the TenTusscherModel2004OdeSystem model
    *
    * calcium_dynamics_Ca_i = 0.0002
    * calcium_dynamics_Ca_SR = 0.2
    * calcium_dynamics_g = 1.0
    * fast_sodium_current_h_gate_h = 0.75
    * fast_sodium_current_j_gate_j = 0.75
    * fast_sodium_current_m_gate_m = 0.0
    * L_type_calcium_current_d_gate_d = 0.0
    * L_type_calcium_current_f_Ca_gate_f_Ca = 1.0
    * L_type_calcium_current_f_gate_f = 1.0
    * membrane_V = -86.2
    * potassium_dynamics_K_i = 138.3
    * rapid_delayed_rectifier_current_X_r1_gate_X_r1 = 0.0
    * rapid_delayed_rectifier_current_X_r2_gate_X_r2 = 1.0
    * slow_delayed_rectifier_current_X_s_gate_X_s = 0.0
    * sodium_dynamics_Na_i = 11.6
    * transient_outward_current_r_gate_r = 0.0
    * transient_outward_current_s_gate_s = 1.0
    */

   /*
    * Throw an exception if the initial vector is larger than the number of equations
    */

   assert(rY.size() == 17);

   double calcium_dynamics_Ca_i = rY[0];
   double calcium_dynamics_Ca_SR = rY[1];
   double calcium_dynamics_g = rY[2];
   double fast_sodium_current_h_gate_h = rY[3];
   double fast_sodium_current_j_gate_j = rY[4];
   double fast_sodium_current_m_gate_m = rY[5];
   double L_type_calcium_current_d_gate_d = rY[6];
   double L_type_calcium_current_f_Ca_gate_f_Ca = rY[7];
   double L_type_calcium_current_f_gate_f = rY[8];
   double membrane_V = rY[9];
   double potassium_dynamics_K_i = rY[10];
   double rapid_delayed_rectifier_current_X_r1_gate_X_r1 = rY[11];
   double rapid_delayed_rectifier_current_X_r2_gate_X_r2 = rY[12];
   double slow_delayed_rectifier_current_X_s_gate_X_s = rY[13];
   double sodium_dynamics_Na_i = rY[14];
   double transient_outward_current_r_gate_r = rY[15];
   double transient_outward_current_s_gate_s = rY[16];

   /*
    * Compute the TenTusscherModel2004OdeSystem model
    */

   double reversal_potentials_E_Ca = (model_parameters_R*model_parameters_T/(2.0*model_parameters_F))*log(model_parameters_Ca_o/calcium_dynamics_Ca_i);
   double calcium_background_current_i_bCa = calcium_background_current_g_bCa*(membrane_V-reversal_potentials_E_Ca);
   double calcium_dynamics_i_rel = (calcium_dynamics_a_rel*pow(calcium_dynamics_Ca_SR,2.0)/(pow(calcium_dynamics_b_rel,2.0)+pow(calcium_dynamics_Ca_SR,2.0))+calcium_dynamics_c_rel)*L_type_calcium_current_d_gate_d*calcium_dynamics_g;
   double calcium_dynamics_i_up = calcium_dynamics_Vmax_up/(1.0+pow(calcium_dynamics_K_up,2.0)/pow(calcium_dynamics_Ca_i,2.0));
   double calcium_dynamics_i_leak = calcium_dynamics_V_leak*(calcium_dynamics_Ca_SR-calcium_dynamics_Ca_i);

   double calcium_dynamics_g_infinity;

   if (calcium_dynamics_Ca_i <= 0.00035)
   {
      calcium_dynamics_g_infinity = 1.0/(1.0+pow(calcium_dynamics_Ca_i, 6.0)/pow(0.00035, 6.0));
   }
   else
   {
      calcium_dynamics_g_infinity = 1.0/(1.0+pow(calcium_dynamics_Ca_i, 16.0)/pow(0.00035, 16.0));
   }

   double calcium_dynamics_k;

   if ((calcium_dynamics_g_infinity > calcium_dynamics_g) && (membrane_V > -60.0))
   {
      calcium_dynamics_k = 0.0;
   }
   else
   {
      calcium_dynamics_k = 1.0;
   }

   double calcium_dynamics_g_prime = calcium_dynamics_k*(calcium_dynamics_g_infinity-calcium_dynamics_g)/calcium_dynamics_tau_g;
   double calcium_dynamics_Ca_i_bufc = calcium_dynamics_Ca_i*calcium_dynamics_Bufc/(calcium_dynamics_Ca_i+calcium_dynamics_Kbufc);
   double calcium_dynamics_Ca_SR_bufsr = calcium_dynamics_Ca_SR*calcium_dynamics_Bufsr/(calcium_dynamics_Ca_SR+calcium_dynamics_Kbufsr);
   double L_type_calcium_current_i_CaL = L_type_calcium_current_g_CaL*L_type_calcium_current_d_gate_d*L_type_calcium_current_f_gate_f*L_type_calcium_current_f_Ca_gate_f_Ca*4.0*membrane_V*pow(model_parameters_F,2.0)/(model_parameters_R*model_parameters_T)*(calcium_dynamics_Ca_i*exp(2.0*membrane_V*model_parameters_F/(model_parameters_R*model_parameters_T))-0.341*model_parameters_Ca_o)/(exp(2.0*membrane_V*model_parameters_F/(model_parameters_R*model_parameters_T))-1.0);
   double calcium_pump_current_i_pCa = calcium_pump_current_g_pCa*calcium_dynamics_Ca_i/(calcium_pump_current_K_pCa+calcium_dynamics_Ca_i);
   double Na_Ca_exchanger_current_i_NaCa = Na_Ca_exchanger_current_K_NaCa*(exp(Na_Ca_exchanger_current_gamma*membrane_V*model_parameters_F/(model_parameters_R*model_parameters_T))*pow(sodium_dynamics_Na_i, 3.0)*model_parameters_Ca_o-exp((Na_Ca_exchanger_current_gamma-1.0)*membrane_V*model_parameters_F/(model_parameters_R*model_parameters_T))*pow(model_parameters_Na_o, 3.0)*calcium_dynamics_Ca_i*Na_Ca_exchanger_current_alpha)/((pow(Na_Ca_exchanger_current_K_mNa_i, 3.0)+pow(model_parameters_Na_o, 3.0))*(Na_Ca_exchanger_current_K_mCa+model_parameters_Ca_o)*(1.0+Na_Ca_exchanger_current_K_sat*exp((Na_Ca_exchanger_current_gamma-1.0)*membrane_V*model_parameters_F/(model_parameters_R*model_parameters_T))));
   double calcium_dynamics_Ca_i_prime = calcium_dynamics_i_leak+calcium_dynamics_i_rel-((L_type_calcium_current_i_CaL+calcium_background_current_i_bCa+calcium_pump_current_i_pCa-2.0*Na_Ca_exchanger_current_i_NaCa)/(2.0*model_parameters_Vc*model_parameters_F)*membrane_Cm+calcium_dynamics_i_up);
   double calcium_dynamics_Ca_SR_prime = (model_parameters_Vc/calcium_dynamics_Vsr)*(calcium_dynamics_i_up-(calcium_dynamics_i_rel+calcium_dynamics_i_leak));
   double reversal_potentials_E_Na = (model_parameters_R*model_parameters_T/model_parameters_F)*log(model_parameters_Na_o/sodium_dynamics_Na_i);
   double fast_sodium_current_i_Na = fast_sodium_current_g_Na*pow(fast_sodium_current_m_gate_m, 3.0)*fast_sodium_current_h_gate_h*fast_sodium_current_j_gate_j*(membrane_V-reversal_potentials_E_Na);
   double fast_sodium_current_h_gate_h_infinity = 1.0/pow((1.0+exp((71.55+membrane_V)/7.43)),2.0);

   double fast_sodium_current_h_gate_alpha_h;

   if (membrane_V < -40.0)
   {
      fast_sodium_current_h_gate_alpha_h = 0.057*exp(-(80.0+membrane_V)/6.8);
   }
   else
   {
      fast_sodium_current_h_gate_alpha_h = 0.0;
   }

   double fast_sodium_current_h_gate_beta_h;

   if (membrane_V < -40.0)
   {
      fast_sodium_current_h_gate_beta_h = 2.7*exp(0.079*membrane_V)+3.1e5 * exp(0.3485*membrane_V);
   }
   else
   {
      fast_sodium_current_h_gate_beta_h = 0.77/(0.13*(1.0+exp(-(membrane_V+10.66)/11.1)));
   }

   double fast_sodium_current_h_gate_tau_h = 1.0/(fast_sodium_current_h_gate_alpha_h+fast_sodium_current_h_gate_beta_h);
   double fast_sodium_current_h_gate_h_prime = fast_sodium_current_h_gate_alpha_h*(1.0-fast_sodium_current_h_gate_h)-fast_sodium_current_h_gate_beta_h*fast_sodium_current_h_gate_h;
   double fast_sodium_current_j_gate_j_infinity = 1.0/pow((1.0+exp((71.55+membrane_V)/7.43)),2.0);

   double fast_sodium_current_j_gate_alpha_j;

   if (membrane_V < -40.0)
   {
      fast_sodium_current_j_gate_alpha_j = (-2.5428e4 * exp(0.2444*membrane_V)-6.948E-6*exp(-0.04391*membrane_V))*(membrane_V+37.78)/(1.0+exp(0.311*(membrane_V+79.23)));
   }
   else
   {
      fast_sodium_current_j_gate_alpha_j = 0.0;
   }

   double fast_sodium_current_j_gate_beta_j;

   if (membrane_V < -40.0)
   {
      fast_sodium_current_j_gate_beta_j = 0.02424*exp(-0.01052*membrane_V)/(1.0+exp(-0.1378*(membrane_V+40.14)));
   }
   else
   {
      fast_sodium_current_j_gate_beta_j = 0.6*exp(0.057*membrane_V)/(1.0+exp(-0.1*(membrane_V+32.0)));
   }

   double fast_sodium_current_j_gate_tau_j = 1.0/(fast_sodium_current_j_gate_alpha_j+fast_sodium_current_j_gate_beta_j);
   double fast_sodium_current_j_gate_j_prime = fast_sodium_current_j_gate_alpha_j*(1.0-fast_sodium_current_j_gate_j)-fast_sodium_current_j_gate_beta_j*fast_sodium_current_j_gate_j;
   
   double fast_sodium_current_m_gate_m_infinity = 1.0/pow((1.0+exp((-56.86-membrane_V)/9.03)),2.0);
   double fast_sodium_current_m_gate_alpha_m = 1.0/(1.0+exp((-60.0-membrane_V)/5.0));
   double fast_sodium_current_m_gate_beta_m = 0.1/(1.0+exp((35.0+membrane_V)/5.0))+0.1/(1.0+exp((membrane_V-50.0)/200.0));
   double fast_sodium_current_m_gate_tau_m = fast_sodium_current_m_gate_alpha_m*fast_sodium_current_m_gate_beta_m;
   double fast_sodium_current_m_gate_m_prime = fast_sodium_current_m_gate_alpha_m*(1.0-fast_sodium_current_m_gate_m)-fast_sodium_current_m_gate_beta_m*fast_sodium_current_m_gate_m;
   double reversal_potentials_E_K = (model_parameters_R*model_parameters_T/model_parameters_F)*log(model_parameters_K_o/potassium_dynamics_K_i);
   double inward_rectifier_potassium_current_K1_gate_alpha_K1 = 0.1/(1.0+exp(0.06*(membrane_V-(reversal_potentials_E_K+200.0))));
   double inward_rectifier_potassium_current_K1_gate_beta_K1 = (3.0*exp(0.0002*(membrane_V-reversal_potentials_E_K+100.0))+exp(0.1*(membrane_V-(reversal_potentials_E_K+10.0))))/(1.0+exp(-0.5*(membrane_V-reversal_potentials_E_K)));
   double inward_rectifier_potassium_current_K1_gate_K1_infinity = inward_rectifier_potassium_current_K1_gate_alpha_K1/(inward_rectifier_potassium_current_K1_gate_alpha_K1+inward_rectifier_potassium_current_K1_gate_beta_K1);
   double inward_rectifier_potassium_current_i_K1 = inward_rectifier_potassium_current_g_K1*sqrt(model_parameters_K_o/5.4)*inward_rectifier_potassium_current_K1_gate_K1_infinity*(membrane_V-reversal_potentials_E_K);
   double L_type_calcium_current_d_gate_d_infinity = 1.0/(1.0+exp((-5.0-membrane_V)/7.5));
   // Change below from /7.5 to /13.0
   double L_type_calcium_current_d_gate_alpha_d = 1.4/(1.0+exp((-35.0-membrane_V)/13.0))+0.25;
   double L_type_calcium_current_d_gate_beta_d = 1.4/(1.0+exp((5.0+membrane_V)/5.0));
   // Change below from 1.4 to 1.0
   double L_type_calcium_current_d_gate_gamma_d = 1.0/(1.0+exp((50.0-membrane_V)/20.0));
   double L_type_calcium_current_d_gate_tau_d = L_type_calcium_current_d_gate_alpha_d*L_type_calcium_current_d_gate_beta_d+L_type_calcium_current_d_gate_gamma_d;
   double L_type_calcium_current_d_gate_d_prime = L_type_calcium_current_d_gate_alpha_d*(1.0-L_type_calcium_current_d_gate_d)-L_type_calcium_current_d_gate_beta_d*L_type_calcium_current_d_gate_d;
   // Line below has an exp removed.............
   double L_type_calcium_current_f_Ca_gate_alpha_f_Ca = 1.0/(1.0+pow(calcium_dynamics_Ca_i/0.000325, 8.0));
   double L_type_calcium_current_f_Ca_gate_beta_f_Ca = 0.1/(1.0+exp((calcium_dynamics_Ca_i-0.0005)/0.0001));
   double L_type_calcium_current_f_Ca_gate_gamma_f_Ca = 0.2/(1.0+exp((calcium_dynamics_Ca_i-0.00075)/0.0008));
   double L_type_calcium_current_f_Ca_gate_f_Ca_infinity = (L_type_calcium_current_f_Ca_gate_alpha_f_Ca+L_type_calcium_current_f_Ca_gate_beta_f_Ca+L_type_calcium_current_f_Ca_gate_gamma_f_Ca+0.23)/1.46;

   double L_type_calcium_current_f_Ca_gate_k;

   if ((L_type_calcium_current_f_Ca_gate_f_Ca_infinity > L_type_calcium_current_f_Ca_gate_f_Ca) && (membrane_V > -60.0))
   {
      L_type_calcium_current_f_Ca_gate_k = 0.0;
   }
   else
   {
      L_type_calcium_current_f_Ca_gate_k = 1.0;
   }

   double L_type_calcium_current_f_Ca_gate_f_Ca_prime = L_type_calcium_current_f_Ca_gate_k*(L_type_calcium_current_f_Ca_gate_f_Ca_infinity-L_type_calcium_current_f_Ca_gate_f_Ca)/L_type_calcium_current_f_Ca_gate_tau_f_Ca;
   double L_type_calcium_current_f_gate_f_infinity = 1.0/(1.0+exp((20.0+membrane_V)/7.0));
   double L_type_calcium_current_f_gate_tau_f = 1125.0*exp(-pow(membrane_V+27.0,2.0)/240.0)+165.0/(1.0+exp((25.0-membrane_V)/10.0))+80.0;
   double L_type_calcium_current_f_gate_f_prime = (L_type_calcium_current_f_gate_f_infinity-L_type_calcium_current_f_gate_f)/L_type_calcium_current_f_gate_tau_f;

   double membrane_i_stim = mpStimulus->GetStimulus(time);

   double transient_outward_current_i_to = transient_outward_current_g_to*transient_outward_current_r_gate_r*transient_outward_current_s_gate_s*(membrane_V-reversal_potentials_E_K);
   double rapid_delayed_rectifier_current_i_Kr = rapid_delayed_rectifier_current_g_Kr*sqrt(model_parameters_K_o/5.4)*rapid_delayed_rectifier_current_X_r1_gate_X_r1*rapid_delayed_rectifier_current_X_r2_gate_X_r2*(membrane_V-reversal_potentials_E_K);
   double reversal_potentials_E_Ks = (model_parameters_R*model_parameters_T/model_parameters_F)*log((model_parameters_K_o+reversal_potentials_p_KNa*model_parameters_Na_o)/(potassium_dynamics_K_i+reversal_potentials_p_KNa*sodium_dynamics_Na_i));
   double slow_delayed_rectifier_current_i_Ks = slow_delayed_rectifier_current_g_Ks*pow(slow_delayed_rectifier_current_X_s_gate_X_s,2.0)*(membrane_V-reversal_potentials_E_Ks);
   double sodium_potassium_pump_current_i_NaK = sodium_potassium_pump_current_P_NaK*model_parameters_K_o*sodium_dynamics_Na_i/((model_parameters_K_o+sodium_potassium_pump_current_K_mK)*(sodium_dynamics_Na_i+sodium_potassium_pump_current_K_mNa)*(1.0+0.1245*exp(-0.1*membrane_V*model_parameters_F/(model_parameters_R*model_parameters_T))+0.0353*exp(-membrane_V*model_parameters_F/(model_parameters_R*model_parameters_T))));
   double sodium_background_current_i_bNa = sodium_background_current_g_bNa*(membrane_V-reversal_potentials_E_Na);
   double potassium_pump_current_i_pK = potassium_pump_current_g_pK*(membrane_V-reversal_potentials_E_K)/(1.0+exp((25.0-membrane_V)/5.98));
   
   double membrane_V_prime = -(fast_sodium_current_i_Na+inward_rectifier_potassium_current_i_K1
        +transient_outward_current_i_to+rapid_delayed_rectifier_current_i_Kr+
        slow_delayed_rectifier_current_i_Ks+L_type_calcium_current_i_CaL+
        Na_Ca_exchanger_current_i_NaCa+sodium_potassium_pump_current_i_NaK+
        sodium_background_current_i_bNa+calcium_background_current_i_bCa+
        calcium_pump_current_i_pCa+potassium_pump_current_i_pK+membrane_i_stim)/membrane_Cm;
   // Took a *membrane_Cm off the end of next one, and it is missing an I_ax , but can't find one of them...
   double potassium_dynamics_K_i_prime = -(inward_rectifier_potassium_current_i_K1+transient_outward_current_i_to+rapid_delayed_rectifier_current_i_Kr+slow_delayed_rectifier_current_i_Ks+potassium_pump_current_i_pK+membrane_i_stim-2.0*sodium_potassium_pump_current_i_NaK)/(model_parameters_Vc*model_parameters_F);
   double rapid_delayed_rectifier_current_X_r1_gate_alpha_X_r1 = 450.0/(1.0+exp((-45.0-membrane_V)/10.0));
   double rapid_delayed_rectifier_current_X_r1_gate_beta_X_r1 = 6.0/(1.0+exp((membrane_V+30.0)/11.5));
   double rapid_delayed_rectifier_current_X_r1_gate_X_r1_prime = rapid_delayed_rectifier_current_X_r1_gate_alpha_X_r1*(1.0-rapid_delayed_rectifier_current_X_r1_gate_X_r1)-rapid_delayed_rectifier_current_X_r1_gate_beta_X_r1*rapid_delayed_rectifier_current_X_r1_gate_X_r1;
   double rapid_delayed_rectifier_current_X_r1_gate_X_r1_infinity = 1.0/(1.0+exp((-26.0-membrane_V)/7.0));
   double rapid_delayed_rectifier_current_X_r1_gate_tau_X_r1 = rapid_delayed_rectifier_current_X_r1_gate_alpha_X_r1*rapid_delayed_rectifier_current_X_r1_gate_beta_X_r1;
   double rapid_delayed_rectifier_current_X_r2_gate_alpha_X_r2 = 3.0/(1.0+exp((-60.0-membrane_V)/20.0));
   double rapid_delayed_rectifier_current_X_r2_gate_beta_X_r2 = 1.12/(1.0+exp((membrane_V-60.0)/20.0));
   double rapid_delayed_rectifier_current_X_r2_gate_X_r2_prime = rapid_delayed_rectifier_current_X_r2_gate_alpha_X_r2*(1.0-rapid_delayed_rectifier_current_X_r2_gate_X_r2)-rapid_delayed_rectifier_current_X_r2_gate_beta_X_r2*rapid_delayed_rectifier_current_X_r2_gate_X_r2;
   double rapid_delayed_rectifier_current_X_r2_gate_X_r2_infinity = 1.0/(1.0+exp((88.0+membrane_V)/24.0));
   double rapid_delayed_rectifier_current_X_r2_gate_tau_X_r2 = rapid_delayed_rectifier_current_X_r2_gate_alpha_X_r2*rapid_delayed_rectifier_current_X_r2_gate_beta_X_r2;
   
   double slow_delayed_rectifier_current_X_s_gate_alpha_X_s = 1100.0/sqrt(1.0+exp((-10.0-membrane_V)/6.0));
   double slow_delayed_rectifier_current_X_s_gate_beta_X_s = 1.0/(1.0+exp((membrane_V-60.0)/20.0));
   double slow_delayed_rectifier_current_X_s_gate_X_s_prime = slow_delayed_rectifier_current_X_s_gate_alpha_X_s*(1.0-slow_delayed_rectifier_current_X_s_gate_X_s)-slow_delayed_rectifier_current_X_s_gate_beta_X_s*slow_delayed_rectifier_current_X_s_gate_X_s;
   double slow_delayed_rectifier_current_X_s_gate_X_s_infinity = 1.0/(1.0+exp((-5.0-membrane_V)/14.0));
   double slow_delayed_rectifier_current_X_s_gate_tau_X_s = slow_delayed_rectifier_current_X_s_gate_alpha_X_s*slow_delayed_rectifier_current_X_s_gate_beta_X_s;
   // took a * Cm off the end of this next one
   double sodium_dynamics_Na_i_prime = -(fast_sodium_current_i_Na+sodium_background_current_i_bNa+3.0*sodium_potassium_pump_current_i_NaK+3.0*Na_Ca_exchanger_current_i_NaCa)/(model_parameters_Vc*model_parameters_F);
   
   double transient_outward_current_r_gate_r_infinity = 1.0/(1.0+exp((20.0-membrane_V)/6.0));   // good
   double transient_outward_current_r_gate_tau_r = 9.5*exp(-(pow(membrane_V+40.0,2.0)/1800.0))+0.8;
   
   double transient_outward_current_r_gate_r_prime = (transient_outward_current_r_gate_r_infinity-transient_outward_current_r_gate_r)/transient_outward_current_r_gate_tau_r;
   
   double transient_outward_current_s_gate_s_infinity = 1.0/(1.0+exp((20.0+membrane_V)/5.0));
   double transient_outward_current_s_gate_tau_s = 85.0*exp(-(pow(membrane_V+45.0,2.0)/320.0))+5.0/(1.0+exp((membrane_V-20.0)/5.0))+3.0;
   
   double transient_outward_current_s_gate_s_prime = (transient_outward_current_s_gate_s_infinity-transient_outward_current_s_gate_s)/transient_outward_current_s_gate_tau_s;
   

   std::vector<double> returnRHS;

   returnRHS.push_back(calcium_dynamics_Ca_i_prime);
   returnRHS.push_back(calcium_dynamics_Ca_SR_prime);
   returnRHS.push_back(calcium_dynamics_g_prime);
   returnRHS.push_back(fast_sodium_current_h_gate_h_prime);
   returnRHS.push_back(fast_sodium_current_j_gate_j_prime);
   returnRHS.push_back(fast_sodium_current_m_gate_m_prime);
   returnRHS.push_back(L_type_calcium_current_d_gate_d_prime);
   returnRHS.push_back(L_type_calcium_current_f_Ca_gate_f_Ca_prime);
   returnRHS.push_back(L_type_calcium_current_f_gate_f_prime);
   returnRHS.push_back(membrane_V_prime);
   returnRHS.push_back(potassium_dynamics_K_i_prime);
   returnRHS.push_back(rapid_delayed_rectifier_current_X_r1_gate_X_r1_prime);
   returnRHS.push_back(rapid_delayed_rectifier_current_X_r2_gate_X_r2_prime);
   returnRHS.push_back(slow_delayed_rectifier_current_X_s_gate_X_s_prime);
   returnRHS.push_back(sodium_dynamics_Na_i_prime);
   returnRHS.push_back(transient_outward_current_r_gate_r_prime);
   returnRHS.push_back(transient_outward_current_s_gate_s_prime);

   return returnRHS;
}
