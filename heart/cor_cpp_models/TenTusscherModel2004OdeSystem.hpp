#ifndef _TENTUSSCHERMODEL2004ODESYSTEM_HPP_
#define _TENTUSSCHERMODEL2004ODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class sets up the TenTusscherModel2004OdeSystem system of equations.
 */
class TenTusscherModel2004OdeSystem : public AbstractOdeSystem
{
   private:
      // Current and voltage components (objects) of the TenTusscherModel2004OdeSystem model
      AbstractStimulusFunction *mpStimulus;

      // Constants for the TenTusscherModel2004OdeSystem model
      double calcium_background_current_g_bCa;
      double calcium_dynamics_a_rel;
      double calcium_dynamics_b_rel;
      double calcium_dynamics_Bufc;
      double calcium_dynamics_Bufsr;
      double calcium_dynamics_c_rel;
      double calcium_dynamics_K_up;
      double calcium_dynamics_Kbufc;
      double calcium_dynamics_Kbufsr;
      double calcium_dynamics_tau_g;
      double calcium_dynamics_V_leak;
      double calcium_dynamics_Vmax_up;
      double calcium_dynamics_Vsr;
      double calcium_pump_current_g_pCa;
      double calcium_pump_current_K_pCa;
      double fast_sodium_current_g_Na;
      double inward_rectifier_potassium_current_g_K1;
      double L_type_calcium_current_f_Ca_gate_tau_f_Ca;
      double L_type_calcium_current_g_CaL;
      double membrane_Cm;
      double model_parameters_Ca_o;
      double model_parameters_F;
      double model_parameters_K_o;
      double model_parameters_Na_o;
      double model_parameters_R;
      double model_parameters_T;
      double model_parameters_Vc;
      double Na_Ca_exchanger_current_alpha;
      double Na_Ca_exchanger_current_gamma;
      double Na_Ca_exchanger_current_K_mCa;
      double Na_Ca_exchanger_current_K_mNa_i;
      double Na_Ca_exchanger_current_K_NaCa;
      double Na_Ca_exchanger_current_K_sat;
      double potassium_pump_current_g_pK;
      double rapid_delayed_rectifier_current_g_Kr;
      double reversal_potentials_p_KNa;
      double slow_delayed_rectifier_current_g_Ks;
      double sodium_background_current_g_bNa;
      double sodium_potassium_pump_current_K_mK;
      double sodium_potassium_pump_current_K_mNa;
      double sodium_potassium_pump_current_P_NaK;
      double transient_outward_current_g_to;

   public:
      // Constructor
      TenTusscherModel2004OdeSystem(AbstractStimulusFunction *stimulus);
      // Destructor
      ~TenTusscherModel2004OdeSystem();

      // This method will compute the RHS of the TenTusscherModel2004OdeSystem model
      std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY);
};

#endif //
