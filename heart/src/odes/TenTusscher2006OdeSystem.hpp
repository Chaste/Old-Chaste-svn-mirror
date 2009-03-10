/*

Copyright (C) University of Oxford, 2005-2009

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

#ifndef _TENTUSSCHER2006ODESYSTEM_HPP_
#define _TENTUSSCHER2006ODESYSTEM_HPP_

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class sets up the equations for the TenTusscher 2006 model. This is for an epicardial cell.
 */
class TenTusscher2006OdeSystem : public AbstractCardiacCell
{
private:

    /** Scale factor for Gks*/ 
    double mScaleFactorGks;
    /** Scale factor for Gto*/
    double mScaleFactorGto;
    
    //////////////////////////////////////////////////////////////
    //Constants for the TenTusscher2006 model, values for epicardial cell.
    //////////////////////////////////////////////////////////////
   static const double L_type_Ca_current_g_CaL = 0.0000398;   // nanoS_per_picoF
   static const double calcium_background_current_g_bca = 0.000592;   // nanoS_per_picoF
   static const double calcium_dynamics_Buf_c = 0.2;   // millimolar
   static const double calcium_dynamics_Buf_sr = 10.0;   // millimolar
   static const double calcium_dynamics_Buf_ss = 0.4;   // millimolar
   static const double calcium_dynamics_Ca_o = 2.0;   // millimolar
   static const double calcium_dynamics_EC = 1.5;   // millimolar
   static const double calcium_dynamics_K_buf_c = 0.001;   // millimolar
   static const double calcium_dynamics_K_buf_sr = 0.3;   // millimolar
   static const double calcium_dynamics_K_buf_ss = 0.00025;   // millimolar
   static const double calcium_dynamics_K_up = 0.00025;   // millimolar
   static const double calcium_dynamics_V_leak = 0.00036;   // millimolar_per_millisecond
   static const double calcium_dynamics_V_rel = 0.102;   // millimolar_per_millisecond
   static const double calcium_dynamics_V_sr = 0.001094;   // micrometre3
   static const double calcium_dynamics_V_ss = 0.00005468;   // micrometre3
   static const double calcium_dynamics_V_xfer = 0.0038;   // millimolar_per_millisecond
   static const double calcium_dynamics_Vmax_up = 0.006375;   // millimolar_per_millisecond
   static const double calcium_dynamics_k1_prime = 0.15;   // per_millimolar2_per_millisecond
   static const double calcium_dynamics_k2_prime = 0.045;   // per_millimolar_per_millisecond
   static const double calcium_dynamics_k3 = 0.06;   // per_millisecond
   static const double calcium_dynamics_k4 = 0.005;   // per_millisecond
   static const double calcium_dynamics_max_sr = 2.5;   // dimensionless
   static const double calcium_dynamics_min_sr = 1.0;   // dimensionless
   static const double calcium_pump_current_K_pCa = 0.0005;   // millimolar
   static const double calcium_pump_current_g_pCa = 0.1238;   // nanoS_per_picoF
   static const double fast_sodium_current_g_Na = 14.838;   // nanoS_per_picoF
   static const double inward_rectifier_potassium_current_g_K1 = 5.405;   // nanoS_per_picoF
   static const double membrane_Cm = 0.185;   // microF_per_cm2
   static const double membrane_F = 96485.3415;   // coulomb_per_millimole
   static const double membrane_R = 8314.472;   // joule_per_mole_kelvin
   static const double membrane_T = 310.0;   // kelvin
   static const double membrane_V_c = 0.016404;   // micrometre3
   static const double potassium_dynamics_K_o = 5.4;   // millimolar
   static const double potassium_pump_current_g_pK = 0.0146;   // nanoS_per_picoF
   static const double rapid_time_dependent_potassium_current_g_Kr = 0.153;   // nanoS_per_picoF
   static const double reversal_potentials_P_kna = 0.03;   // nanoA_per_millimolar
   static const double slow_time_dependent_potassium_current_g_Ks = 0.392;   // nanoS_per_picoF
   static const double sodium_background_current_g_bna = 0.00029;   // nanoS_per_picoF
   static const double sodium_calcium_exchanger_current_K_NaCa = 1000.0;   // picoA_per_picoF
   static const double sodium_calcium_exchanger_current_K_sat = 0.1;   // dimensionless
   static const double sodium_calcium_exchanger_current_Km_Ca = 1.38;   // millimolar
   static const double sodium_calcium_exchanger_current_Km_Nai = 87.5;   // millimolar
   static const double sodium_calcium_exchanger_current_alpha = 2.5;   // dimensionless
   static const double sodium_calcium_exchanger_current_gamma = 0.35;   // dimensionless
   static const double sodium_dynamics_Na_o = 140.0;   // millimolar
   static const double sodium_potassium_pump_current_K_mNa = 40.0;   // millimolar
   static const double sodium_potassium_pump_current_K_mk = 1.0;   // millimolar
   static const double sodium_potassium_pump_current_P_NaK = 2.724;   // picoA_per_picoF
   static const double transient_outward_current_g_to = 0.294;   // nanoS_per_picoF

   ///////////////////////////////////////////////////////////////////////
   //variables that need to be computed
   //////////////////////////////////////////////////////////////////////
      double L_type_Ca_current_d_gate_alpha_d;   // per_millisecond
      double L_type_Ca_current_d_gate_beta_d;   // per_millisecond
      double L_type_Ca_current_d_gate_d_inf;   // dimensionless
      double L_type_Ca_current_d_gate_gamma_d;   // per_millisecond
      double L_type_Ca_current_d_gate_tau_d;   // millisecond
      double L_type_Ca_current_f2_gate_f2_inf;   // dimensionless
      double L_type_Ca_current_f2_gate_tau_f2;   // millisecond
      double L_type_Ca_current_fCass_gate_fCass_inf;   // dimensionless
      double L_type_Ca_current_fCass_gate_tau_fCass;   // millisecond
      double L_type_Ca_current_f_gate_f_inf;   // dimensionless
      double L_type_Ca_current_f_gate_tau_f;   // millisecond
      double L_type_Ca_current_i_CaL;   // picoA_per_picoF
      double calcium_background_current_i_b_Ca;   // picoA_per_picoF
      double calcium_dynamics_Ca_i_bufc;   // millimolar
      double calcium_dynamics_Ca_sr_bufsr;   // millimolar
      double calcium_dynamics_Ca_ss_bufss;   // millimolar
      double calcium_dynamics_O;   // dimensionless
      double calcium_dynamics_i_leak;   // millimolar_per_millisecond
      double calcium_dynamics_i_rel;   // millimolar_per_millisecond
      double calcium_dynamics_i_up;   // millimolar_per_millisecond
      double calcium_dynamics_i_xfer;   // millimolar_per_millisecond
      double calcium_dynamics_k1;   // per_millimolar2_per_millisecond
      double calcium_dynamics_k2;   // per_millimolar_per_millisecond
      double calcium_dynamics_kcasr;   // dimensionless
      double calcium_pump_current_i_p_Ca;   // picoA_per_picoF
      double fast_sodium_current_h_gate_alpha_h;   // per_millisecond
      double fast_sodium_current_h_gate_beta_h;   // per_millisecond
      double fast_sodium_current_h_gate_h_inf;   // dimensionless
      double fast_sodium_current_h_gate_tau_h;   // millisecond
      double fast_sodium_current_i_Na;   // picoA_per_picoF
      double fast_sodium_current_j_gate_alpha_j;   // per_millisecond
      double fast_sodium_current_j_gate_beta_j;   // per_millisecond
      double fast_sodium_current_j_gate_j_inf;   // dimensionless
      double fast_sodium_current_j_gate_tau_j;   // millisecond
      double fast_sodium_current_m_gate_alpha_m;   // per_millisecond
      double fast_sodium_current_m_gate_beta_m;   // per_millisecond
      double fast_sodium_current_m_gate_m_inf;   // dimensionless
      double fast_sodium_current_m_gate_tau_m;   // millisecond
      double inward_rectifier_potassium_current_alpha_K1;   // dimensionless
      double inward_rectifier_potassium_current_beta_K1;   // dimensionless
      double inward_rectifier_potassium_current_i_K1;   // picoA_per_picoF
      double inward_rectifier_potassium_current_xK1_inf;   // dimensionless
      double potassium_pump_current_i_p_K;   // picoA_per_picoF
      double rapid_time_dependent_potassium_current_Xr1_gate_alpha_xr1;   // per_millisecond
      double rapid_time_dependent_potassium_current_Xr1_gate_beta_xr1;   // per_millisecond
      double rapid_time_dependent_potassium_current_Xr1_gate_tau_xr1;   // millisecond
      double rapid_time_dependent_potassium_current_Xr1_gate_xr1_inf;   // dimensionless
      double rapid_time_dependent_potassium_current_Xr2_gate_alpha_xr2;   // per_millisecond
      double rapid_time_dependent_potassium_current_Xr2_gate_beta_xr2;   // per_millisecond
      double rapid_time_dependent_potassium_current_Xr2_gate_tau_xr2;   // millisecond
      double rapid_time_dependent_potassium_current_Xr2_gate_xr2_inf;   // dimensionless
      double rapid_time_dependent_potassium_current_i_Kr;   // picoA_per_picoF
      double reversal_potentials_E_Ca;   // millivolt
      double reversal_potentials_E_K;   // millivolt
      double reversal_potentials_E_Ks;   // millivolt
      double reversal_potentials_E_Na;   // millivolt
      double slow_time_dependent_potassium_current_Xs_gate_alpha_xs;   // per_millisecond
      double slow_time_dependent_potassium_current_Xs_gate_beta_xs;   // per_millisecond
      double slow_time_dependent_potassium_current_Xs_gate_tau_xs;   // millisecond
      double slow_time_dependent_potassium_current_Xs_gate_xs_inf;   // dimensionless
      double slow_time_dependent_potassium_current_i_Ks;   // picoA_per_picoF
      double sodium_background_current_i_b_Na;   // picoA_per_picoF
      double sodium_calcium_exchanger_current_i_NaCa;   // picoA_per_picoF
      double sodium_potassium_pump_current_i_NaK;   // picoA_per_picoF
      double transient_outward_current_i_to;   // picoA_per_picoF
      double transient_outward_current_r_gate_r_inf;   // dimensionless
      double transient_outward_current_r_gate_tau_r;   // millisecond
      double transient_outward_current_s_gate_s_inf;   // dimensionless
      double transient_outward_current_s_gate_tau_s;   // dimensionless

    // This private method will check that gates are within 0 and 1 and concentrations are positive
    void VerifyStateVariables();

public:
    /**
     * Constructor
     */ 
    TenTusscher2006OdeSystem(AbstractIvpOdeSolver *pSolver,
                               AbstractStimulusFunction *pIntracellularStimulus);

    /**
     * Destructor
     */ 
    ~TenTusscher2006OdeSystem();

    /**
     *  This method will compute the RHS of the TenTusscher model
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    
    /**
     * Set the scale factor for Gks in order to differentiate epi M and endo cells
     */
    void SetScaleFactorGks(double sfgks);
    
    /**
     * Set the scale factor for Gks in order to differentiate epi M and endo cells
     */
    void SetScaleFactorGto(double sfgto);
    
     /**
     * Returns the ionic current 
     */
    double GetIIonic();

};

#endif // _TENTUSSCHER2006_HPP_
