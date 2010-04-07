/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#include "TenTusscher2006OdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include <cmath>
#include <cstdio>
#include "Exception.hpp"


//
// Model-scope constant parameters
//
const double TenTusscher2006OdeSystem::L_type_Ca_current_g_CaL = 0.0000398;   // nanoS_per_picoF
const double TenTusscher2006OdeSystem::calcium_background_current_g_bca = 0.000592;   // nanoS_per_picoF
const double TenTusscher2006OdeSystem::calcium_dynamics_Buf_c = 0.2;   // millimolar
const double TenTusscher2006OdeSystem::calcium_dynamics_Buf_sr = 10.0;   // millimolar
const double TenTusscher2006OdeSystem::calcium_dynamics_Buf_ss = 0.4;   // millimolar
const double TenTusscher2006OdeSystem::calcium_dynamics_Ca_o = 2.0;   // millimolar
const double TenTusscher2006OdeSystem::calcium_dynamics_EC = 1.5;   // millimolar
const double TenTusscher2006OdeSystem::calcium_dynamics_K_buf_c = 0.001;   // millimolar
const double TenTusscher2006OdeSystem::calcium_dynamics_K_buf_sr = 0.3;   // millimolar
const double TenTusscher2006OdeSystem::calcium_dynamics_K_buf_ss = 0.00025;   // millimolar
const double TenTusscher2006OdeSystem::calcium_dynamics_K_up = 0.00025;   // millimolar
const double TenTusscher2006OdeSystem::calcium_dynamics_V_leak = 0.00036;   // millimolar_per_millisecond
const double TenTusscher2006OdeSystem::calcium_dynamics_V_rel = 0.102;   // millimolar_per_millisecond
const double TenTusscher2006OdeSystem::calcium_dynamics_V_sr = 0.001094;   // micrometre3
const double TenTusscher2006OdeSystem::calcium_dynamics_V_ss = 0.00005468;   // micrometre3
const double TenTusscher2006OdeSystem::calcium_dynamics_V_xfer = 0.0038;   // millimolar_per_millisecond
const double TenTusscher2006OdeSystem::calcium_dynamics_Vmax_up = 0.006375;   // millimolar_per_millisecond
const double TenTusscher2006OdeSystem::calcium_dynamics_k1_prime = 0.15;   // per_millimolar2_per_millisecond
const double TenTusscher2006OdeSystem::calcium_dynamics_k2_prime = 0.045;   // per_millimolar_per_millisecond
const double TenTusscher2006OdeSystem::calcium_dynamics_k3 = 0.06;   // per_millisecond
const double TenTusscher2006OdeSystem::calcium_dynamics_k4 = 0.005;   // per_millisecond
const double TenTusscher2006OdeSystem::calcium_dynamics_max_sr = 2.5;   // dimensionless
const double TenTusscher2006OdeSystem::calcium_dynamics_min_sr = 1.0;   // dimensionless
const double TenTusscher2006OdeSystem::calcium_pump_current_K_pCa = 0.0005;   // millimolar
const double TenTusscher2006OdeSystem::calcium_pump_current_g_pCa = 0.1238;   // nanoS_per_picoF
const double TenTusscher2006OdeSystem::fast_sodium_current_g_Na = 14.838;   // nanoS_per_picoF
const double TenTusscher2006OdeSystem::inward_rectifier_potassium_current_g_K1 = 5.405;   // nanoS_per_picoF
const double TenTusscher2006OdeSystem::membrane_Cm = 0.185;   // microF_per_cm2
const double TenTusscher2006OdeSystem::membrane_F = 96485.3415;   // coulomb_per_millimole
const double TenTusscher2006OdeSystem::membrane_R = 8314.472;   // joule_per_mole_kelvin
const double TenTusscher2006OdeSystem::membrane_T = 310.0;   // kelvin
const double TenTusscher2006OdeSystem::membrane_V_c = 0.016404;   // micrometre3
const double TenTusscher2006OdeSystem::potassium_dynamics_K_o = 5.4;   // millimolar
const double TenTusscher2006OdeSystem::potassium_pump_current_g_pK = 0.0146;   // nanoS_per_picoF
const double TenTusscher2006OdeSystem::rapid_time_dependent_potassium_current_g_Kr = 0.153;   // nanoS_per_picoF
const double TenTusscher2006OdeSystem::reversal_potentials_P_kna = 0.03;   // nanoA_per_millimolar
const double TenTusscher2006OdeSystem::slow_time_dependent_potassium_current_g_Ks = 0.392;   // nanoS_per_picoF
const double TenTusscher2006OdeSystem::sodium_background_current_g_bna = 0.00029;   // nanoS_per_picoF
const double TenTusscher2006OdeSystem::sodium_calcium_exchanger_current_K_NaCa = 1000.0;   // picoA_per_picoF
const double TenTusscher2006OdeSystem::sodium_calcium_exchanger_current_K_sat = 0.1;   // dimensionless
const double TenTusscher2006OdeSystem::sodium_calcium_exchanger_current_Km_Ca = 1.38;   // millimolar
const double TenTusscher2006OdeSystem::sodium_calcium_exchanger_current_Km_Nai = 87.5;   // millimolar
const double TenTusscher2006OdeSystem::sodium_calcium_exchanger_current_alpha = 2.5;   // dimensionless
const double TenTusscher2006OdeSystem::sodium_calcium_exchanger_current_gamma = 0.35;   // dimensionless
const double TenTusscher2006OdeSystem::sodium_dynamics_Na_o = 140.0;   // millimolar
const double TenTusscher2006OdeSystem::sodium_potassium_pump_current_K_mNa = 40.0;   // millimolar
const double TenTusscher2006OdeSystem::sodium_potassium_pump_current_K_mk = 1.0;   // millimolar
const double TenTusscher2006OdeSystem::sodium_potassium_pump_current_P_NaK = 2.724;   // picoA_per_picoF
const double TenTusscher2006OdeSystem::transient_outward_current_g_to = 0.294;   // nanoS_per_picoF



/*Constructor*/
TenTusscher2006OdeSystem::TenTusscher2006OdeSystem(
        boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
        boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : AbstractCardiacCell(pSolver, 19, 11, pIntracellularStimulus)
{

    mpSystemInfo = OdeSystemInformation<TenTusscher2006OdeSystem>::Instance();
    //Initialise the scale factors
    mScaleFactorGks=1.0;
    mScaleFactorIto=1.0;
    mScaleFactorGkr=1.0;

    Init();
}

//Destructor
TenTusscher2006OdeSystem::~TenTusscher2006OdeSystem(void)
{
}

void TenTusscher2006OdeSystem::SetScaleFactorGks(double sfgks)
{
    mScaleFactorGks=sfgks;
}
void TenTusscher2006OdeSystem::SetScaleFactorIto(double sfito)
{
    mScaleFactorIto=sfito;
}
void TenTusscher2006OdeSystem::SetScaleFactorGkr(double sfgkr)
{
    mScaleFactorGkr=sfgkr;
}

void TenTusscher2006OdeSystem::EvaluateYDerivatives(double time,
                                                    const std::vector<double> &rY,
                                                    std::vector<double> &rDY)
{
   //---------------------------------------------------------------------------
   // State variables, initial value and names in comments beside
   //---------------------------------------------------------------------------
      //Vector for state variables//
      double Y[19];
      double dY[19];

    Y[0] = rY[0];// 3.373e-5;   // L_t Ype_Ca_current_d_gate_d (dimensionless)
    Y[1] = rY[1];// 0.9755;   // L_t Ype_Ca_current_f2_gate_f2 (dimensionless)
    Y[2] = rY[2];// 0.9953;   // L_t Ype_Ca_current_fCass_gate_fCass (dimensionless)
    Y[3] = rY[3];// 0.7888;   // L_t Ype_Ca_current_f_gate_f (dimensionless)
    Y[4] = rY[4];// 3.64;   // calcium_d Ynamics_Ca_SR (millimolar)
    Y[5] = rY[5];// 0.000126;   // calcium_d Ynamics_Ca_i (millimolar)
    Y[6] = rY[6];// 0.00036;   // calcium_d Ynamics_Ca_ss (millimolar)
    Y[7] = rY[7];// 0.9073;   // calcium_d Ynamics_R_prime (dimensionless)
    Y[8] = rY[8];// 0.7444;   // fast_sodium_current_h_gate_h (dimensionless)
    Y[9] = rY[9];// 0.7045;   // fast_sodium_current_j_gate_j (dimensionless)
    Y[10] = rY[10];// 0.00172;   // fast_sodium_current_m_gate_m (dimensionless)
    Y[11] = rY[11];// -85.23;   // membrane_V (millivolt)
    Y[12] = rY[12];// 136.89;   // potassium_d Ynamics_K_i (millimolar)
    Y[13] = rY[13];// 0.00621;   // rapid_time_dependent_potassium_current_Xr1_gate_Xr1 (dimensionless)
    Y[14] = rY[14];// 0.4712;   // rapid_time_dependent_potassium_current_Xr2_gate_Xr2 (dimensionless)
    Y[15] = rY[15];// 0.0095;   // slow_time_dependent_potassium_current_Xs_gate_Xs (dimensionless)
    Y[16] = rY[16];// 8.604;   // sodium_d Ynamics_Na_i (millimolar)
    Y[17] = rY[17];// 2.42e-8;   // transient_outward_current_r_gate_r (dimensionless)
    Y[18] = rY[18];// 0.999998;   // transient_outward_current_s_gate_s (dimensionless)

    VerifyStateVariables();

   double L_type_Ca_current_i_CaL = L_type_Ca_current_g_CaL*Y[0]*Y[3]*Y[1]*Y[2]*4.0*(Y[11]-15.0)*pow(membrane_F, 2.0)/(membrane_R*membrane_T)*(0.25*Y[6]*exp(2.0*(Y[11]-15.0)*membrane_F/(membrane_R*membrane_T))-calcium_dynamics_Ca_o)/(exp(2.0*(Y[11]-15.0)*membrane_F/(membrane_R*membrane_T))-1.0);
   double L_type_Ca_current_d_gate_d_inf = 1.0/(1.0+exp((-8.0-Y[11])/7.5));
   double L_type_Ca_current_d_gate_alpha_d = 1.4/(1.0+exp((-35.0-Y[11])/13.0))+0.25;
   double L_type_Ca_current_d_gate_beta_d = 1.4/(1.0+exp((Y[11]+5.0)/5.0));
   double L_type_Ca_current_d_gate_gamma_d = 1.0/(1.0+exp((50.0-Y[11])/20.0));
   double L_type_Ca_current_d_gate_tau_d = L_type_Ca_current_d_gate_alpha_d*L_type_Ca_current_d_gate_beta_d+L_type_Ca_current_d_gate_gamma_d;
   dY[0] = (L_type_Ca_current_d_gate_d_inf-Y[0])/L_type_Ca_current_d_gate_tau_d;
   double L_type_Ca_current_f2_gate_f2_inf = 0.67/(1.0+exp((Y[11]+35.0)/7.0))+0.33;
   double L_type_Ca_current_f2_gate_tau_f2 = 562.0*exp(-pow(Y[11]+27.0, 2.0)/240.0)+31.0/(1.0+exp((25.0-Y[11])/10.0))+80.0/(1.0+exp((Y[11]+30.0)/10.0));
   dY[1] = (L_type_Ca_current_f2_gate_f2_inf-Y[1])/L_type_Ca_current_f2_gate_tau_f2;
   double L_type_Ca_current_fCass_gate_fCass_inf = 0.6/(1.0+pow(Y[6]/0.05, 2.0))+0.4;
   double L_type_Ca_current_fCass_gate_tau_fCass = 80.0/(1.0+pow(Y[6]/0.05, 2.0))+2.0;
   dY[2] = (L_type_Ca_current_fCass_gate_fCass_inf-Y[2])/L_type_Ca_current_fCass_gate_tau_fCass;
   double L_type_Ca_current_f_gate_f_inf = 1.0/(1.0+exp((Y[11]+20.0)/7.0));
   double L_type_Ca_current_f_gate_tau_f = 1102.5*exp(-pow(Y[11]+27.0, 2.0)/225.0)+200.0/(1.0+exp((13.0-Y[11])/10.0))+180.0/(1.0+exp((Y[11]+30.0)/10.0))+20.0;
   dY[3] = (L_type_Ca_current_f_gate_f_inf-Y[3])/L_type_Ca_current_f_gate_tau_f;
   double reversal_potentials_E_Ca = 0.5*membrane_R*membrane_T/membrane_F*log(calcium_dynamics_Ca_o/Y[5]);
   double calcium_background_current_i_b_Ca = calcium_background_current_g_bca*(Y[11]-reversal_potentials_E_Ca);
   double calcium_dynamics_kcasr = calcium_dynamics_max_sr-(calcium_dynamics_max_sr-calcium_dynamics_min_sr)/(1.0+pow(calcium_dynamics_EC/Y[4], 2.0));
   double calcium_dynamics_k1 = calcium_dynamics_k1_prime/calcium_dynamics_kcasr;
   double calcium_dynamics_O = calcium_dynamics_k1*pow(Y[6], 2.0)*Y[7]/(calcium_dynamics_k3+calcium_dynamics_k1*pow(Y[6], 2.0));
   double calcium_dynamics_i_rel = calcium_dynamics_V_rel*calcium_dynamics_O*(Y[4]-Y[6]);
   double calcium_dynamics_i_up = calcium_dynamics_Vmax_up/(1.0+pow(calcium_dynamics_K_up, 2.0)/pow(Y[5], 2.0));
   double calcium_dynamics_i_leak = calcium_dynamics_V_leak*(Y[4]-Y[5]);
   double calcium_dynamics_i_xfer = calcium_dynamics_V_xfer*(Y[6]-Y[5]);
   double calcium_dynamics_k2 = calcium_dynamics_k2_prime*calcium_dynamics_kcasr;
   dY[7] = -calcium_dynamics_k2*Y[6]*Y[7]+calcium_dynamics_k4*(1.0-Y[7]);
   double calcium_dynamics_Ca_i_bufc = 1.0/(1.0+calcium_dynamics_Buf_c*calcium_dynamics_K_buf_c/pow(Y[5]+calcium_dynamics_K_buf_c, 2.0));
   double calcium_dynamics_Ca_sr_bufsr = 1.0/(1.0+calcium_dynamics_Buf_sr*calcium_dynamics_K_buf_sr/pow(Y[4]+calcium_dynamics_K_buf_sr, 2.0));
   double calcium_dynamics_Ca_ss_bufss = 1.0/(1.0+calcium_dynamics_Buf_ss*calcium_dynamics_K_buf_ss/pow(Y[6]+calcium_dynamics_K_buf_ss, 2.0));
   double calcium_pump_current_i_p_Ca = calcium_pump_current_g_pCa*Y[5]/(Y[5]+calcium_pump_current_K_pCa);
   double sodium_calcium_exchanger_current_i_NaCa = sodium_calcium_exchanger_current_K_NaCa*(exp(sodium_calcium_exchanger_current_gamma*Y[11]*membrane_F/(membrane_R*membrane_T))*pow(Y[16], 3.0)*calcium_dynamics_Ca_o-exp((sodium_calcium_exchanger_current_gamma-1.0)*Y[11]*membrane_F/(membrane_R*membrane_T))*pow(sodium_dynamics_Na_o, 3.0)*Y[5]*sodium_calcium_exchanger_current_alpha)/((pow(sodium_calcium_exchanger_current_Km_Nai, 3.0)+pow(sodium_dynamics_Na_o, 3.0))*(sodium_calcium_exchanger_current_Km_Ca+calcium_dynamics_Ca_o)*(1.0+sodium_calcium_exchanger_current_K_sat*exp((sodium_calcium_exchanger_current_gamma-1.0)*Y[11]*membrane_F/(membrane_R*membrane_T))));
   dY[5] = calcium_dynamics_Ca_i_bufc*((calcium_dynamics_i_leak-calcium_dynamics_i_up)*calcium_dynamics_V_sr/membrane_V_c+calcium_dynamics_i_xfer-(calcium_background_current_i_b_Ca+calcium_pump_current_i_p_Ca-2.0*sodium_calcium_exchanger_current_i_NaCa)*membrane_Cm/(2.0*membrane_V_c*membrane_F));
   dY[4] = calcium_dynamics_Ca_sr_bufsr*(calcium_dynamics_i_up-(calcium_dynamics_i_rel+calcium_dynamics_i_leak));
   dY[6] = calcium_dynamics_Ca_ss_bufss*(-L_type_Ca_current_i_CaL*membrane_Cm/(2.0*calcium_dynamics_V_ss*membrane_F)+calcium_dynamics_i_rel*calcium_dynamics_V_sr/calcium_dynamics_V_ss-calcium_dynamics_i_xfer*membrane_V_c/calcium_dynamics_V_ss);
   double reversal_potentials_E_Na = membrane_R*membrane_T/membrane_F*log(sodium_dynamics_Na_o/Y[16]);
   double fast_sodium_current_i_Na = fast_sodium_current_g_Na*pow(Y[10], 3.0)*Y[8]*Y[9]*(Y[11]-reversal_potentials_E_Na);
   double fast_sodium_current_h_gate_h_inf = 1.0/pow(1.0+exp((Y[11]+71.55)/7.43), 2.0);

   double fast_sodium_current_h_gate_alpha_h,fast_sodium_current_h_gate_beta_h;
   if (Y[11] < -40.0)
      fast_sodium_current_h_gate_alpha_h = 0.057*exp(-(Y[11]+80.0)/6.8);
   else
      fast_sodium_current_h_gate_alpha_h = 0.0;

   if (Y[11] < -40.0)
      fast_sodium_current_h_gate_beta_h = 2.7*exp(0.079*Y[11])+310000.0*exp(0.3485*Y[11]);
   else
      fast_sodium_current_h_gate_beta_h = 0.77/(0.13*(1.0+exp((Y[11]+10.66)/-11.1)));

   double fast_sodium_current_h_gate_tau_h = 1.0/(fast_sodium_current_h_gate_alpha_h+fast_sodium_current_h_gate_beta_h);
   dY[8] = (fast_sodium_current_h_gate_h_inf-Y[8])/fast_sodium_current_h_gate_tau_h;
   double fast_sodium_current_j_gate_j_inf = 1.0/pow(1.0+exp((Y[11]+71.55)/7.43), 2.0);

   double fast_sodium_current_j_gate_alpha_j ,fast_sodium_current_j_gate_beta_j;
   if (Y[11] < -40.0)
      fast_sodium_current_j_gate_alpha_j = (-25428.0*exp(0.2444*Y[11])-6.948e-6*exp(-0.04391*Y[11]))*(Y[11]+37.78)/(1.0+exp(0.311*(Y[11]+79.23)));
   else
      fast_sodium_current_j_gate_alpha_j = 0.0;

   if (Y[11] < -40.0)
      fast_sodium_current_j_gate_beta_j = 0.02424*exp(-0.01052*Y[11])/(1.0+exp(-0.1378*(Y[11]+40.14)));
   else
      fast_sodium_current_j_gate_beta_j = 0.6*exp(0.057*Y[11])/(1.0+exp(-0.1*(Y[11]+32.0)));

   double fast_sodium_current_j_gate_tau_j = 1.0/(fast_sodium_current_j_gate_alpha_j+fast_sodium_current_j_gate_beta_j);
   dY[9] = (fast_sodium_current_j_gate_j_inf-Y[9])/fast_sodium_current_j_gate_tau_j;
   double fast_sodium_current_m_gate_m_inf = 1.0/pow(1.0+exp((-56.86-Y[11])/9.03), 2.0);
   double fast_sodium_current_m_gate_alpha_m = 1.0/(1.0+exp((-60.0-Y[11])/5.0));
   double fast_sodium_current_m_gate_beta_m = 0.1/(1.0+exp((Y[11]+35.0)/5.0))+0.1/(1.0+exp((Y[11]-50.0)/200.0));
   double fast_sodium_current_m_gate_tau_m = fast_sodium_current_m_gate_alpha_m*fast_sodium_current_m_gate_beta_m;
   dY[10] = (fast_sodium_current_m_gate_m_inf-Y[10])/fast_sodium_current_m_gate_tau_m;
   double reversal_potentials_E_K = membrane_R*membrane_T/membrane_F*log(potassium_dynamics_K_o/Y[12]);
   double inward_rectifier_potassium_current_alpha_K1 = 0.1/(1.0+exp(0.06*(Y[11]-reversal_potentials_E_K-200.0)));
   double inward_rectifier_potassium_current_beta_K1 = (3.0*exp(0.0002*(Y[11]-reversal_potentials_E_K+100.0))+exp(0.1*(Y[11]-reversal_potentials_E_K-10.0)))/(1.0+exp(-0.5*(Y[11]-reversal_potentials_E_K)));
   double inward_rectifier_potassium_current_xK1_inf = inward_rectifier_potassium_current_alpha_K1/(inward_rectifier_potassium_current_alpha_K1+inward_rectifier_potassium_current_beta_K1);
   double inward_rectifier_potassium_current_i_K1 = inward_rectifier_potassium_current_g_K1*inward_rectifier_potassium_current_xK1_inf*(Y[11]-reversal_potentials_E_K);


   double transient_outward_current_i_to = mScaleFactorIto*transient_outward_current_g_to*Y[17]*Y[18]*(Y[11]-reversal_potentials_E_K);
   double rapid_time_dependent_potassium_current_i_Kr = mScaleFactorGkr*rapid_time_dependent_potassium_current_g_Kr*sqrt(potassium_dynamics_K_o/5.4)*Y[13]*Y[14]*(Y[11]-reversal_potentials_E_K);
   double reversal_potentials_E_Ks = membrane_R*membrane_T/membrane_F*log((potassium_dynamics_K_o+reversal_potentials_P_kna*sodium_dynamics_Na_o)/(Y[12]+reversal_potentials_P_kna*Y[16]));
   double slow_time_dependent_potassium_current_i_Ks = mScaleFactorGks*slow_time_dependent_potassium_current_g_Ks*pow(Y[15], 2.0)*(Y[11]-reversal_potentials_E_Ks);
   double sodium_potassium_pump_current_i_NaK = sodium_potassium_pump_current_P_NaK*potassium_dynamics_K_o/(potassium_dynamics_K_o+sodium_potassium_pump_current_K_mk)*Y[16]/(Y[16]+sodium_potassium_pump_current_K_mNa)/(1.0+0.1245*exp(-0.1*Y[11]*membrane_F/(membrane_R*membrane_T))+0.0353*exp(-Y[11]*membrane_F/(membrane_R*membrane_T)));
   double sodium_background_current_i_b_Na = sodium_background_current_g_bna*(Y[11]-reversal_potentials_E_Na);
   double potassium_pump_current_i_p_K = potassium_pump_current_g_pK*(Y[11]-reversal_potentials_E_K)/(1.0+exp((25.0-Y[11])/5.98));

   //stimulus current
   double i_stim = GetStimulus(time);

   dY[11] = -1.0/1.0*(inward_rectifier_potassium_current_i_K1+transient_outward_current_i_to+rapid_time_dependent_potassium_current_i_Kr+slow_time_dependent_potassium_current_i_Ks+L_type_Ca_current_i_CaL+sodium_potassium_pump_current_i_NaK+fast_sodium_current_i_Na+sodium_background_current_i_b_Na+sodium_calcium_exchanger_current_i_NaCa+calcium_background_current_i_b_Ca+potassium_pump_current_i_p_K+calcium_pump_current_i_p_Ca+i_stim);
   dY[12] = -(inward_rectifier_potassium_current_i_K1+transient_outward_current_i_to+rapid_time_dependent_potassium_current_i_Kr+slow_time_dependent_potassium_current_i_Ks+potassium_pump_current_i_p_K+i_stim-2.0*sodium_potassium_pump_current_i_NaK)/(membrane_V_c*membrane_F)*membrane_Cm;
   double rapid_time_dependent_potassium_current_Xr1_gate_xr1_inf = 1.0/(1.0+exp((-26.0-Y[11])/7.0));
   double rapid_time_dependent_potassium_current_Xr1_gate_alpha_xr1 = 450.0/(1.0+exp((-45.0-Y[11])/10.0));
   double rapid_time_dependent_potassium_current_Xr1_gate_beta_xr1 = 6.0/(1.0+exp((Y[11]+30.0)/11.5));
   double rapid_time_dependent_potassium_current_Xr1_gate_tau_xr1 = rapid_time_dependent_potassium_current_Xr1_gate_alpha_xr1*rapid_time_dependent_potassium_current_Xr1_gate_beta_xr1;
   dY[13] = (rapid_time_dependent_potassium_current_Xr1_gate_xr1_inf-Y[13])/rapid_time_dependent_potassium_current_Xr1_gate_tau_xr1;
   double rapid_time_dependent_potassium_current_Xr2_gate_xr2_inf = 1.0/(1.0+exp((Y[11]+88.0)/24.0));
   double rapid_time_dependent_potassium_current_Xr2_gate_alpha_xr2 = 3.0/(1.0+exp((-60.0-Y[11])/20.0));
   double rapid_time_dependent_potassium_current_Xr2_gate_beta_xr2 = 1.12/(1.0+exp((Y[11]-60.0)/20.0));
   double rapid_time_dependent_potassium_current_Xr2_gate_tau_xr2 = rapid_time_dependent_potassium_current_Xr2_gate_alpha_xr2*rapid_time_dependent_potassium_current_Xr2_gate_beta_xr2;
   dY[14] = (rapid_time_dependent_potassium_current_Xr2_gate_xr2_inf-Y[14])/rapid_time_dependent_potassium_current_Xr2_gate_tau_xr2;
   double slow_time_dependent_potassium_current_Xs_gate_xs_inf = 1.0/(1.0+exp((-5.0-Y[11])/14.0));
   double slow_time_dependent_potassium_current_Xs_gate_alpha_xs = 1400.0/sqrt(1.0+exp((5.0-Y[11])/6.0));
   double slow_time_dependent_potassium_current_Xs_gate_beta_xs = 1.0/(1.0+exp((Y[11]-35.0)/15.0));
   double slow_time_dependent_potassium_current_Xs_gate_tau_xs = slow_time_dependent_potassium_current_Xs_gate_alpha_xs*slow_time_dependent_potassium_current_Xs_gate_beta_xs+80.0;
   dY[15] = (slow_time_dependent_potassium_current_Xs_gate_xs_inf-Y[15])/slow_time_dependent_potassium_current_Xs_gate_tau_xs;
   dY[16] = -(fast_sodium_current_i_Na+sodium_background_current_i_b_Na+3.0*sodium_potassium_pump_current_i_NaK+3.0*sodium_calcium_exchanger_current_i_NaCa)/(membrane_V_c*membrane_F)*membrane_Cm;
   double transient_outward_current_r_gate_r_inf = 1.0/(1.0+exp((20.0-Y[11])/6.0));
   double transient_outward_current_r_gate_tau_r = 9.5*exp(-pow(Y[11]+40.0, 2.0)/1800.0)+0.8;
   dY[17] = (transient_outward_current_r_gate_r_inf-Y[17])/transient_outward_current_r_gate_tau_r;
   double transient_outward_current_s_gate_s_inf = 1.0/(1.0+exp((Y[11]+20.0)/5.0));
   double transient_outward_current_s_gate_tau_s = 85.0*exp(-pow(Y[11]+45.0, 2.0)/320.0)+5.0/(1.0+exp((Y[11]-20.0)/5.0))+3.0;
   dY[18] = (transient_outward_current_s_gate_s_inf-Y[18])/transient_outward_current_s_gate_tau_s;


    // do not update voltage if the mSetVoltageDerivativeToZero flag has been set
    if (mSetVoltageDerivativeToZero)
    {
        dY[11] = 0;
    }

    rDY[0] = dY[0];
    rDY[1] = dY[1];
    rDY[2] = dY[2];
    rDY[3] = dY[3];
    rDY[4] = dY[4];
    rDY[5] = dY[5];
    rDY[6] = dY[6];
    rDY[7] = dY[7];
    rDY[8] = dY[8];
    rDY[9] = dY[9];
    rDY[10]= dY[10];
    rDY[11] = dY[11];
    rDY[12] = dY[12];
    rDY[13] = dY[13];
    rDY[14] = dY[14];
    rDY[15] = dY[15];
    rDY[16] = dY[16];
    rDY[17] = dY[17];
    rDY[18] = dY[18];
}


double TenTusscher2006OdeSystem::GetIIonic()
{
      //Vector for state variables//
      double Y[19];

    Y[0] =  mStateVariables[0];// L_t Ype_Ca_current_d_gate_d (dimensionless)
    Y[1] =  mStateVariables[1];// L_t Ype_Ca_current_f2_gate_f2 (dimensionless)
    Y[2] =  mStateVariables[2];// L_t Ype_Ca_current_fCass_gate_fCass (dimensionless)
    Y[3] =  mStateVariables[3];// L_t Ype_Ca_current_f_gate_f (dimensionless)
    Y[4] =  mStateVariables[4];// calcium_d Ynamics_Ca_SR (millimolar)
    Y[5] =  mStateVariables[5];// calcium_d Ynamics_Ca_i (millimolar)
    Y[6] =  mStateVariables[6];// calcium_d Ynamics_Ca_ss (millimolar)
    Y[7] =  mStateVariables[7];// calcium_d Ynamics_R_prime (dimensionless)
    Y[8] =  mStateVariables[8];// fast_sodium_current_h_gate_h (dimensionless)
    Y[9] =  mStateVariables[9];// fast_sodium_current_j_gate_j (dimensionless)
    Y[10] =  mStateVariables[10];// fast_sodium_current_m_gate_m (dimensionless)
    Y[11] =  mStateVariables[11];// membrane_V (millivolt)
    Y[12] =  mStateVariables[12];// potassium_d Ynamics_K_i (millimolar)
    Y[13] =  mStateVariables[13];// rapid_time_dependent_potassium_current_Xr1_gate_Xr1 (dimensionless)
    Y[14] =  mStateVariables[14];// rapid_time_dependent_potassium_current_Xr2_gate_Xr2 (dimensionless)
    Y[15] =  mStateVariables[15];// slow_time_dependent_potassium_current_Xs_gate_Xs (dimensionless)
    Y[16] =  mStateVariables[16];// sodium_d Ynamics_Na_i (millimolar)
    Y[17] =  mStateVariables[17];// transient_outward_current_r_gate_r (dimensionless)
    Y[18] =  mStateVariables[18];// transient_outward_current_s_gate_s (dimensionless)

   double L_type_Ca_current_i_CaL = L_type_Ca_current_g_CaL*Y[0]*Y[3]*Y[1]*Y[2]*4.0*(Y[11]-15.0)*pow(membrane_F, 2.0)/(membrane_R*membrane_T)*(0.25*Y[6]*exp(2.0*(Y[11]-15.0)*membrane_F/(membrane_R*membrane_T))-calcium_dynamics_Ca_o)/(exp(2.0*(Y[11]-15.0)*membrane_F/(membrane_R*membrane_T))-1.0);
   double reversal_potentials_E_Ca = 0.5*membrane_R*membrane_T/membrane_F*log(calcium_dynamics_Ca_o/Y[5]);
   double calcium_background_current_i_b_Ca = calcium_background_current_g_bca*(Y[11]-reversal_potentials_E_Ca);
   double calcium_pump_current_i_p_Ca = calcium_pump_current_g_pCa*Y[5]/(Y[5]+calcium_pump_current_K_pCa);
   double sodium_calcium_exchanger_current_i_NaCa = sodium_calcium_exchanger_current_K_NaCa*(exp(sodium_calcium_exchanger_current_gamma*Y[11]*membrane_F/(membrane_R*membrane_T))*pow(Y[16], 3.0)*calcium_dynamics_Ca_o-exp((sodium_calcium_exchanger_current_gamma-1.0)*Y[11]*membrane_F/(membrane_R*membrane_T))*pow(sodium_dynamics_Na_o, 3.0)*Y[5]*sodium_calcium_exchanger_current_alpha)/((pow(sodium_calcium_exchanger_current_Km_Nai, 3.0)+pow(sodium_dynamics_Na_o, 3.0))*(sodium_calcium_exchanger_current_Km_Ca+calcium_dynamics_Ca_o)*(1.0+sodium_calcium_exchanger_current_K_sat*exp((sodium_calcium_exchanger_current_gamma-1.0)*Y[11]*membrane_F/(membrane_R*membrane_T))));
   double reversal_potentials_E_Na = membrane_R*membrane_T/membrane_F*log(sodium_dynamics_Na_o/Y[16]);
   double fast_sodium_current_i_Na = fast_sodium_current_g_Na*pow(Y[10], 3.0)*Y[8]*Y[9]*(Y[11]-reversal_potentials_E_Na);

   double reversal_potentials_E_K = membrane_R*membrane_T/membrane_F*log(potassium_dynamics_K_o/Y[12]);
   double inward_rectifier_potassium_current_alpha_K1 = 0.1/(1.0+exp(0.06*(Y[11]-reversal_potentials_E_K-200.0)));
   double inward_rectifier_potassium_current_beta_K1 = (3.0*exp(0.0002*(Y[11]-reversal_potentials_E_K+100.0))+exp(0.1*(Y[11]-reversal_potentials_E_K-10.0)))/(1.0+exp(-0.5*(Y[11]-reversal_potentials_E_K)));
   double inward_rectifier_potassium_current_xK1_inf = inward_rectifier_potassium_current_alpha_K1/(inward_rectifier_potassium_current_alpha_K1+inward_rectifier_potassium_current_beta_K1);
   double inward_rectifier_potassium_current_i_K1 = inward_rectifier_potassium_current_g_K1*inward_rectifier_potassium_current_xK1_inf*(Y[11]-reversal_potentials_E_K);

   double transient_outward_current_i_to = mScaleFactorIto*transient_outward_current_g_to*Y[17]*Y[18]*(Y[11]-reversal_potentials_E_K);
   double rapid_time_dependent_potassium_current_i_Kr = mScaleFactorGkr*rapid_time_dependent_potassium_current_g_Kr*sqrt(potassium_dynamics_K_o/5.4)*Y[13]*Y[14]*(Y[11]-reversal_potentials_E_K);
   double reversal_potentials_E_Ks = membrane_R*membrane_T/membrane_F*log((potassium_dynamics_K_o+reversal_potentials_P_kna*sodium_dynamics_Na_o)/(Y[12]+reversal_potentials_P_kna*Y[16]));
   double slow_time_dependent_potassium_current_i_Ks = mScaleFactorGks*slow_time_dependent_potassium_current_g_Ks*pow(Y[15], 2.0)*(Y[11]-reversal_potentials_E_Ks);
   double sodium_potassium_pump_current_i_NaK = sodium_potassium_pump_current_P_NaK*potassium_dynamics_K_o/(potassium_dynamics_K_o+sodium_potassium_pump_current_K_mk)*Y[16]/(Y[16]+sodium_potassium_pump_current_K_mNa)/(1.0+0.1245*exp(-0.1*Y[11]*membrane_F/(membrane_R*membrane_T))+0.0353*exp(-Y[11]*membrane_F/(membrane_R*membrane_T)));
   double sodium_background_current_i_b_Na = sodium_background_current_g_bna*(Y[11]-reversal_potentials_E_Na);
   double potassium_pump_current_i_p_K = potassium_pump_current_g_pK*(Y[11]-reversal_potentials_E_K)/(1.0+exp((25.0-Y[11])/5.98));

    double i_ionic = inward_rectifier_potassium_current_i_K1+transient_outward_current_i_to+rapid_time_dependent_potassium_current_i_Kr+slow_time_dependent_potassium_current_i_Ks+L_type_Ca_current_i_CaL+sodium_potassium_pump_current_i_NaK+fast_sodium_current_i_Na+sodium_background_current_i_b_Na+sodium_calcium_exchanger_current_i_NaCa+calcium_background_current_i_b_Ca+potassium_pump_current_i_p_K+calcium_pump_current_i_p_Ca; /*this is in nA*/

    assert(!std::isnan(i_ionic));

    double i_ionic_in_microA_per_cm2=i_ionic*1.0;
    return i_ionic_in_microA_per_cm2;

    /*   i_ionic for this model is in pA/pF.
     *    Please note that in the mono/bidomain formulation, i_ionic needs to be in microA/cm2.
     *    We then need to divide by the cell capacitance.
     *    The cell capacitance of the tenTusscher model is
     *    2.0 uF/cm2 in the paper
     *    0.185 uF/cm2 in this code (membrane_C)
     *    1.0 uF/cm2 in the EvaluateRhsDerivatives method above
     *
     *    For consistency, we choose the last option.
     *    i_ion*pow(10,-6) will be in microA/pF.
     *    Cm*pow(10,6) will be in pF/cm2.
     *    i_ion*pow(10,-6)*Cm*pow(10,6) = i_ion*Cm is in microA/cm2, the correct units
     */
}

void TenTusscher2006OdeSystem::VerifyStateVariables()
{
    const std::vector<double>& rY = rGetStateVariables();

      double Y[19];

    Y[0] = rY[0];//  L_t Ype_Ca_current_d_gate_d (dimensionless)
    Y[1] = rY[1];//  L_t Ype_Ca_current_f2_gate_f2 (dimensionless)
    Y[2] = rY[2];//  L_t Ype_Ca_current_fCass_gate_fCass (dimensionless)
    Y[3] = rY[3];//  L_t Ype_Ca_current_f_gate_f (dimensionless)
    Y[4] = rY[4];//  calcium_d Ynamics_Ca_SR (millimolar)
    Y[5] = rY[5];//  calcium_d Ynamics_Ca_i (millimolar)
    Y[6] = rY[6];//  calcium_d Ynamics_Ca_ss (millimolar)
    Y[7] = rY[7];//  calcium_d Ynamics_R_prime (dimensionless)
    Y[8] = rY[8];//  fast_sodium_current_h_gate_h (dimensionless)
    Y[9] = rY[9];//  fast_sodium_current_j_gate_j (dimensionless)
    Y[10] = rY[10];//fast_sodium_current_m_gate_m (dimensionless)
    Y[11] = rY[11];//membrane_V (millivolt)
    Y[12] = rY[12];//potassium_d Ynamics_K_i (millimolar)
    Y[13] = rY[13];//rapid_time_dependent_potassium_current_Xr1_gate_Xr1 (dimensionless)
    Y[14] = rY[14];//rapid_time_dependent_potassium_current_Xr2_gate_Xr2 (dimensionless)
    Y[15] = rY[15];//slow_time_dependent_potassium_current_Xs_gate_Xs (dimensionless)
    Y[16] = rY[16];//sodium_d Ynamics_Na_i (millimolar)
    Y[17] = rY[17];//transient_outward_current_r_gate_r (dimensionless)
    Y[18] = rY[18];//transient_outward_current_s_gate_s (dimensionless)

    #define COVERAGE_IGNORE
    if (!(0<=Y[0] && Y[0]<=1))
    {
        EXCEPTION(DumpState("d gate of L type calcium channel is out of range!"));
    }
    if (!(0<=Y[1] && Y[1]<=1))
    {
        EXCEPTION(DumpState("f2 gate of L type calcium channel is out of range!"));
    }
    if (!(0<=Y[2] && Y[2]<=1))
    {
        EXCEPTION(DumpState("fCa gate of L type calcium channel is out of range!"));
    }
    if (!(0<=Y[3] && Y[3]<=1))
    {
        EXCEPTION(DumpState("f gate of L type calcium channel is out of range!"));
    }
    if (!(0<=Y[4]))
    {
        EXCEPTION(DumpState("CaSR is negative!"));
    }
    if (!(0<=Y[5]))
    {
        EXCEPTION(DumpState("Cai is negative!"));
    }
    if (!(0<=Y[6]))
    {
        EXCEPTION(DumpState("CaSS is negative!"));
    }
    if (!(0<=Y[7] && Y[7]<=1))
    {
        EXCEPTION(DumpState("R gate is out of range!"));
    }
    if (!(0<=Y[8] && Y[8]<=1))
    {
        EXCEPTION(DumpState("h gate of Na channel is out of range!"));
    }
    if (!(0<=Y[9] && Y[9]<=1))
    {
        EXCEPTION(DumpState("j gate of Na channel is out of range!"));
    }
    if (!(0<=Y[10] && Y[10]<=1))
    {
        EXCEPTION(DumpState("m gate of Na channel is out of range!"));
    }
    if (!(-200<=Y[11] && Y[11]<=200))
    {
        EXCEPTION(DumpState("Vm is REALLY out of range!"));
    }
    if (!(0<=Y[12]))
    {
        EXCEPTION(DumpState("Ki is negative!"));
    }
    if (!(0<=Y[13] && Y[13]<=1))
    {
        EXCEPTION(DumpState("xr1 gate of IKR channel is out of range!"));
    }
    if (!(0<=Y[14] && Y[14]<=1))
    {
        EXCEPTION(DumpState("xr2 gate of IKR channel is out of range!"));
    }
    if (!(0<=Y[15] && Y[15]<=1))
    {
        EXCEPTION(DumpState("xs1 gate of IKS channel is out of range!"));
    }
    if (!(0<=Y[16]))
    {
        EXCEPTION(DumpState("Nai is negative!"));
    }
    if (!(0<=Y[17] && Y[17]<=1))
    {
        EXCEPTION(DumpState("r gate of Ito channel is out of range!"));
    }
    if (!(0<=Y[18] && Y[18]<=1))
    {
        EXCEPTION(DumpState("s gate of Ito channel is out of range!"));
    }

    #undef COVERAGE_IGNORE
}

template<>
void OdeSystemInformation<TenTusscher2006OdeSystem>::Initialise(void)
{
    //Initailisation values are for steady state pacing at 500 ms BCL, please note that values may differ for different BCLs
   this->mVariableNames.push_back("d_gate_L");
   this->mVariableUnits.push_back("");
   this->mInitialConditions.push_back(0.00003373);

   this->mVariableNames.push_back("f2_gate_L");
   this->mVariableUnits.push_back("");
   this->mInitialConditions.push_back(0.9755);

   this->mVariableNames.push_back("fca_gate_L");
   this->mVariableUnits.push_back("");
   this->mInitialConditions.push_back(0.9953);

   this->mVariableNames.push_back("f_gate_L");
   this->mVariableUnits.push_back("");
   this->mInitialConditions.push_back(0.7888);

   this->mVariableNames.push_back("CaSR");
   this->mVariableUnits.push_back("mM");
   this->mInitialConditions.push_back(3.64);

   this->mVariableNames.push_back("Cai");
   this->mVariableUnits.push_back("mM");
   this->mInitialConditions.push_back(0.000126);

   this->mVariableNames.push_back("CaSS");
   this->mVariableUnits.push_back("mM");
   this->mInitialConditions.push_back(0.00036);

   this->mVariableNames.push_back("R_gate");
   this->mVariableUnits.push_back("");
   this->mInitialConditions.push_back(0.9073);

   this->mVariableNames.push_back("h_gate_Na");
   this->mVariableUnits.push_back("");
   this->mInitialConditions.push_back(0.7444);

   this->mVariableNames.push_back("j_gate_Na");
   this->mVariableUnits.push_back("");
   this->mInitialConditions.push_back(0.7045);

   this->mVariableNames.push_back("m_gate_Na");
   this->mVariableUnits.push_back("");
   this->mInitialConditions.push_back(0.00172);

   this->mVariableNames.push_back("V");
   this->mVariableUnits.push_back("mV");
   this->mInitialConditions.push_back(-85.23);

   this->mVariableNames.push_back("Ki");
   this->mVariableUnits.push_back("mM");
   this->mInitialConditions.push_back(136.89);

   this->mVariableNames.push_back("xr1_gate");
   this->mVariableUnits.push_back("");
   this->mInitialConditions.push_back(0.00621);

   this->mVariableNames.push_back("xr2_gate");
   this->mVariableUnits.push_back("");
   this->mInitialConditions.push_back(0.4712);

   this->mVariableNames.push_back("xS_gate");
   this->mVariableUnits.push_back("");
   this->mInitialConditions.push_back(0.0095);

   this->mVariableNames.push_back("Nai");
   this->mVariableUnits.push_back("mM");
   this->mInitialConditions.push_back(8.604);

   this->mVariableNames.push_back("r_gate_Ito");
   this->mVariableUnits.push_back("");
   this->mInitialConditions.push_back(0.0000000242);

   this->mVariableNames.push_back("s_gate_Ito");
   this->mVariableUnits.push_back("");
   this->mInitialConditions.push_back(0.999998);

    this->mInitialised = true;
}



