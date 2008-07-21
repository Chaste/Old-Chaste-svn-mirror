#ifndef _CML_noble_varghese_kohl_noble_1998_basic_backward_
#define _CML_noble_varghese_kohl_noble_1998_basic_backward_

// Model: noble_varghese_kohl_noble_1998_basic
// Processed by pycml - CellML Tools in Python
//     (translate: 4024, pycml: 4024)
// on Wed Jul 16 15:08:14 2008

#include <cmath>
#include <cassert>
#include "AbstractBackwardEulerCardiacCell.hpp"
#include "CardiacNewtonSolver.hpp"
#include "Exception.hpp"
#include "AbstractStimulusFunction.hpp"

class CML_noble_varghese_kohl_noble_1998_basic_backward : public AbstractBackwardEulerCardiacCell<12>
{
public:
    CML_noble_varghese_kohl_noble_1998_basic_backward(AbstractStimulusFunction *pIntracellularStimulus,
                                                      AbstractStimulusFunction *pExtracellularStimulus=NULL)
        : AbstractBackwardEulerCardiacCell<12>(22, 0, pIntracellularStimulus, pExtracellularStimulus)
    {
        // Time units: second
        // 
        mVariableNames.push_back("V");
        mVariableUnits.push_back("millivolt");
        mInitialConditions.push_back(-92.849333);

        mVariableNames.push_back("xr1");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(1.03e-5);

        mVariableNames.push_back("xr2");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(2e-7);

        mVariableNames.push_back("xs");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.001302);

        mVariableNames.push_back("m");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.0016203);

        mVariableNames.push_back("h");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.9944036);

        mVariableNames.push_back("d");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0);

        mVariableNames.push_back("f");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(1);

        mVariableNames.push_back("f2");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.9349197);

        mVariableNames.push_back("f2ds");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.9651958);

        mVariableNames.push_back("s");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.9948645);

        mVariableNames.push_back("r");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0);

        mVariableNames.push_back("ActFrac");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.0042614);

        mVariableNames.push_back("ProdFrac");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.4068154);

        mVariableNames.push_back("Na_i");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(7.3321223);

        mVariableNames.push_back("K_i");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(136.5644281);

        mVariableNames.push_back("Ca_i");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(1.4e-5);

        mVariableNames.push_back("Ca_ds");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(1.88e-5);

        mVariableNames.push_back("Ca_up");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(0.4531889);

        mVariableNames.push_back("Ca_rel");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(0.4481927);

        mVariableNames.push_back("Ca_Calmod");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(0.0005555);

        mVariableNames.push_back("Ca_Trop");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(0.0003542);

        Init();

    }

    ~CML_noble_varghese_kohl_noble_1998_basic_backward(void)
    {
    }

    void VerifyGatingVariables() {}
    void VerifyStateVariables() {}

    double GetIIonic()
    {
        std::vector<double>& rY = rGetStateVariables();
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -92.849333
        double var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1 = rY[1];
        // Units: dimensionless; Initial value: 1.03e-5
        double var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2 = rY[2];
        // Units: dimensionless; Initial value: 2e-7
        double var_slow_delayed_rectifier_potassium_current_xs_gate__xs = rY[3];
        // Units: dimensionless; Initial value: 0.001302
        double var_fast_sodium_current_m_gate__m = rY[4];
        // Units: dimensionless; Initial value: 0.0016203
        double var_fast_sodium_current_h_gate__h = rY[5];
        // Units: dimensionless; Initial value: 0.9944036
        double var_L_type_Ca_channel_d_gate__d = rY[6];
        // Units: dimensionless; Initial value: 0
        double var_L_type_Ca_channel_f_gate__f = rY[7];
        // Units: dimensionless; Initial value: 1
        double var_L_type_Ca_channel_f2_gate__f2 = rY[8];
        // Units: dimensionless; Initial value: 0.9349197
        double var_L_type_Ca_channel_f2ds_gate__f2ds = rY[9];
        // Units: dimensionless; Initial value: 0.9651958
        double var_transient_outward_current_s_gate__s = rY[10];
        // Units: dimensionless; Initial value: 0.9948645
        double var_transient_outward_current_r_gate__r = rY[11];
        // Units: dimensionless; Initial value: 0
        double var_intracellular_sodium_concentration__Na_i = rY[14];
        // Units: millimolar; Initial value: 7.3321223
        double var_intracellular_potassium_concentration__K_i = rY[15];
        // Units: millimolar; Initial value: 136.5644281
        double var_intracellular_calcium_concentration__Ca_i = rY[16];
        // Units: millimolar; Initial value: 1.4e-5
        double var_intracellular_calcium_concentration__Ca_ds = rY[17];
        // Units: millimolar; Initial value: 1.88e-5
        
        const double var_membrane__R = 8314.472;
        const double var_membrane__T = 310.0;
        const double var_membrane__F = 96485.3415;
        double var_reversal_potentials__K_i = var_intracellular_potassium_concentration__K_i;
        double var_reversal_potentials__R = var_membrane__R;
        double var_reversal_potentials__T = var_membrane__T;
        double var_reversal_potentials__F = var_membrane__F;
        const double var_extracellular_potassium_concentration__K_o = 4.0;
        double var_reversal_potentials__K_o = var_extracellular_potassium_concentration__K_o;
        double var_reversal_potentials__E_K = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__K_o / var_reversal_potentials__K_i);
        double var_time_independent_potassium_current__E_K = var_reversal_potentials__E_K;
        double var_time_independent_potassium_current__K_o = var_extracellular_potassium_concentration__K_o;
        double var_time_independent_potassium_current__R = var_membrane__R;
        double var_time_independent_potassium_current__V = var_membrane__V;
        double var_time_independent_potassium_current__T = var_membrane__T;
        const double var_time_independent_potassium_current__K_mk1 = 10.0;
        const double var_time_independent_potassium_current__g_K1 = 0.5;
        double var_time_independent_potassium_current__F = var_membrane__F;
        double var_time_independent_potassium_current__i_K1 = (((var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K_o) / (var_time_independent_potassium_current__K_o + var_time_independent_potassium_current__K_mk1)) * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K)) / (1.0 + exp((((var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K) - 10.0) * var_time_independent_potassium_current__F * 1.25) / (var_time_independent_potassium_current__R * var_time_independent_potassium_current__T)));
        double var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
        double var_transient_outward_current__s = var_transient_outward_current_s_gate__s;
        double var_transient_outward_current__r = var_transient_outward_current_r_gate__r;
        const double var_transient_outward_current__g_to = 0.005;
        double var_transient_outward_current__V = var_membrane__V;
        double var_transient_outward_current__E_K = var_reversal_potentials__E_K;
        const double var_transient_outward_current__g_tos = 0.0;
        double var_transient_outward_current__i_to = var_transient_outward_current__g_to * (var_transient_outward_current__g_tos + (var_transient_outward_current__s * (1.0 - var_transient_outward_current__g_tos))) * var_transient_outward_current__r * (var_transient_outward_current__V - var_transient_outward_current__E_K);
        double var_membrane__i_to = var_transient_outward_current__i_to;
        const double var_rapid_delayed_rectifier_potassium_current__g_Kr2 = 0.0013;
        const double var_rapid_delayed_rectifier_potassium_current__g_Kr1 = 0.0021;
        double var_rapid_delayed_rectifier_potassium_current__xr1 = var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1;
        double var_rapid_delayed_rectifier_potassium_current__xr2 = var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2;
        double var_rapid_delayed_rectifier_potassium_current__V = var_membrane__V;
        double var_rapid_delayed_rectifier_potassium_current__E_K = var_reversal_potentials__E_K;
        double var_rapid_delayed_rectifier_potassium_current__i_Kr = ((((var_rapid_delayed_rectifier_potassium_current__g_Kr1 * var_rapid_delayed_rectifier_potassium_current__xr1) + (var_rapid_delayed_rectifier_potassium_current__g_Kr2 * var_rapid_delayed_rectifier_potassium_current__xr2)) * 1.0) / (1.0 + exp((var_rapid_delayed_rectifier_potassium_current__V + 9.0) / 22.4))) * (var_rapid_delayed_rectifier_potassium_current__V - var_rapid_delayed_rectifier_potassium_current__E_K);
        double var_membrane__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
        double var_slow_delayed_rectifier_potassium_current__xs = var_slow_delayed_rectifier_potassium_current_xs_gate__xs;
        const double var_extracellular_sodium_concentration__Na_o = 140.0;
        double var_reversal_potentials__Na_o = var_extracellular_sodium_concentration__Na_o;
        double var_reversal_potentials__Na_i = var_intracellular_sodium_concentration__Na_i;
        const double var_reversal_potentials__P_kna = 0.03;
        double var_reversal_potentials__E_Ks = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log((var_reversal_potentials__K_o + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_o)) / (var_reversal_potentials__K_i + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_i)));
        double var_slow_delayed_rectifier_potassium_current__E_Ks = var_reversal_potentials__E_Ks;
        const double var_slow_delayed_rectifier_potassium_current__g_Ks = 0.0026;
        double var_slow_delayed_rectifier_potassium_current__V = var_membrane__V;
        double var_slow_delayed_rectifier_potassium_current__i_Ks = var_slow_delayed_rectifier_potassium_current__g_Ks * pow(var_slow_delayed_rectifier_potassium_current__xs, 2.0) * (var_slow_delayed_rectifier_potassium_current__V - var_slow_delayed_rectifier_potassium_current__E_Ks);
        double var_membrane__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
        double var_L_type_Ca_channel__d = var_L_type_Ca_channel_d_gate__d;
        const double var_L_type_Ca_channel__FrICa = 1.0;
        double var_L_type_Ca_channel__f = var_L_type_Ca_channel_f_gate__f;
        double var_L_type_Ca_channel__K_o = var_extracellular_potassium_concentration__K_o;
        double var_L_type_Ca_channel__K_i = var_intracellular_potassium_concentration__K_i;
        double var_L_type_Ca_channel__F = var_membrane__F;
        const double var_L_type_Ca_channel__P_Ca_L = 0.1;
        double var_L_type_Ca_channel__T = var_membrane__T;
        const double var_L_type_Ca_channel__P_CaK = 0.002;
        double var_L_type_Ca_channel__V = var_membrane__V;
        double var_L_type_Ca_channel__f2 = var_L_type_Ca_channel_f2_gate__f2;
        double var_L_type_Ca_channel__R = var_membrane__R;
        double var_L_type_Ca_channel__i_Ca_L_K_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaK * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__K_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__K_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_membrane__i_Ca_L_K_cyt = var_L_type_Ca_channel__i_Ca_L_K_cyt;
        double var_L_type_Ca_channel__f2ds = var_L_type_Ca_channel_f2ds_gate__f2ds;
        double var_L_type_Ca_channel__i_Ca_L_K_ds = (((var_L_type_Ca_channel__FrICa * var_L_type_Ca_channel__P_CaK * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__K_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__K_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_membrane__i_Ca_L_K_ds = var_L_type_Ca_channel__i_Ca_L_K_ds;
        const double var_sodium_potassium_pump__i_NaK_max = 0.7;
        double var_sodium_potassium_pump__Na_i = var_intracellular_sodium_concentration__Na_i;
        double var_sodium_potassium_pump__K_o = var_extracellular_potassium_concentration__K_o;
        const double var_sodium_potassium_pump__K_mNa = 40.0;
        const double var_sodium_potassium_pump__K_mK = 1.0;
        double var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__i_NaK_max * var_sodium_potassium_pump__K_o) / (var_sodium_potassium_pump__K_mK + var_sodium_potassium_pump__K_o)) * var_sodium_potassium_pump__Na_i) / (var_sodium_potassium_pump__K_mNa + var_sodium_potassium_pump__Na_i);
        double var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
        const double var_fast_sodium_current__g_Na = 2.5;
        double var_fast_sodium_current__h = var_fast_sodium_current_h_gate__h;
        double var_fast_sodium_current__V = var_membrane__V;
        double var_reversal_potentials__E_mh = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log((var_reversal_potentials__Na_o + (0.12 * var_reversal_potentials__K_o)) / (var_reversal_potentials__Na_i + (0.12 * var_reversal_potentials__K_i)));
        double var_fast_sodium_current__E_mh = var_reversal_potentials__E_mh;
        double var_fast_sodium_current__m = var_fast_sodium_current_m_gate__m;
        double var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * (var_fast_sodium_current__V - var_fast_sodium_current__E_mh);
        double var_membrane__i_Na = var_fast_sodium_current__i_Na;
        double var_sodium_background_current__V = var_membrane__V;
        double var_reversal_potentials__E_Na = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Na_o / var_reversal_potentials__Na_i);
        double var_sodium_background_current__E_Na = var_reversal_potentials__E_Na;
        const double var_sodium_background_current__g_bna = 0.0006;
        double var_sodium_background_current__i_b_Na = var_sodium_background_current__g_bna * (var_sodium_background_current__V - var_sodium_background_current__E_Na);
        double var_membrane__i_b_Na = var_sodium_background_current__i_b_Na;
        const double var_persistent_sodium_current__g_pna = 0.004;
        double var_persistent_sodium_current__V = var_membrane__V;
        double var_persistent_sodium_current__E_Na = var_reversal_potentials__E_Na;
        double var_persistent_sodium_current__i_p_Na = ((var_persistent_sodium_current__g_pna * 1.0) / (1.0 + exp((-(var_persistent_sodium_current__V + 52.0)) / 8.0))) * (var_persistent_sodium_current__V - var_persistent_sodium_current__E_Na);
        double var_membrane__i_p_Na = var_persistent_sodium_current__i_p_Na;
        const double var_L_type_Ca_channel__P_CaNa = 0.01;
        double var_L_type_Ca_channel__Na_o = var_extracellular_sodium_concentration__Na_o;
        double var_L_type_Ca_channel__Na_i = var_intracellular_sodium_concentration__Na_i;
        double var_L_type_Ca_channel__i_Ca_L_Na_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaNa * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Na_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Na_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_membrane__i_Ca_L_Na_cyt = var_L_type_Ca_channel__i_Ca_L_Na_cyt;
        double var_L_type_Ca_channel__i_Ca_L_Na_ds = (((var_L_type_Ca_channel__FrICa * var_L_type_Ca_channel__P_CaNa * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Na_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Na_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_membrane__i_Ca_L_Na_ds = var_L_type_Ca_channel__i_Ca_L_Na_ds;
        double var_sodium_calcium_exchanger__Na_i = var_intracellular_sodium_concentration__Na_i;
        const double var_sodium_calcium_exchanger__n_NaCa = 3.0;
        const double var_sodium_calcium_exchanger__gamma = 0.5;
        double var_sodium_calcium_exchanger__F = var_membrane__F;
        double var_sodium_calcium_exchanger__Na_o = var_extracellular_sodium_concentration__Na_o;
        const double var_sodium_calcium_exchanger__FRiNaCa = 0.001;
        double var_sodium_calcium_exchanger__R = var_membrane__R;
        double var_sodium_calcium_exchanger__Ca_i = var_intracellular_calcium_concentration__Ca_i;
        double var_sodium_calcium_exchanger__T = var_membrane__T;
        double var_sodium_calcium_exchanger__V = var_membrane__V;
        const double var_sodium_calcium_exchanger__d_NaCa = 0.0;
        const double var_extracellular_calcium_concentration__Ca_o = 2.0;
        double var_sodium_calcium_exchanger__Ca_o = var_extracellular_calcium_concentration__Ca_o;
        const double var_sodium_calcium_exchanger__k_NaCa = 0.0005;
        double var_sodium_calcium_exchanger__i_NaCa_cyt = ((1.0 - var_sodium_calcium_exchanger__FRiNaCa) * var_sodium_calcium_exchanger__k_NaCa * ((exp((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_o) - (exp(((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_i))) / ((1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_sodium_calcium_exchanger__Ca_i * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_sodium_calcium_exchanger__Ca_o * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa))))) * (1.0 + (var_sodium_calcium_exchanger__Ca_i / 0.0069)));
        double var_membrane__i_NaCa_cyt = var_sodium_calcium_exchanger__i_NaCa_cyt;
        double var_sodium_calcium_exchanger__Ca_ds = var_intracellular_calcium_concentration__Ca_ds;
        double var_sodium_calcium_exchanger__i_NaCa_ds = (var_sodium_calcium_exchanger__FRiNaCa * var_sodium_calcium_exchanger__k_NaCa * ((exp((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_o) - (exp(((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_ds))) / ((1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_sodium_calcium_exchanger__Ca_ds * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_sodium_calcium_exchanger__Ca_o * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa))))) * (1.0 + (var_sodium_calcium_exchanger__Ca_ds / 0.0069)));
        double var_membrane__i_NaCa_ds = var_sodium_calcium_exchanger__i_NaCa_ds;
        double var_L_type_Ca_channel__Ca_i = var_intracellular_calcium_concentration__Ca_i;
        double var_L_type_Ca_channel__Ca_o = var_extracellular_calcium_concentration__Ca_o;
        double var_L_type_Ca_channel__i_Ca_L_Ca_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * 4.0 * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Ca_i * exp((100.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Ca_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_membrane__i_Ca_L_Ca_cyt = var_L_type_Ca_channel__i_Ca_L_Ca_cyt;
        double var_L_type_Ca_channel__i_Ca_L_Ca_ds = (((var_L_type_Ca_channel__FrICa * 4.0 * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Ca_i * exp((100.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Ca_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_membrane__i_Ca_L_Ca_ds = var_L_type_Ca_channel__i_Ca_L_Ca_ds;
        double var_reversal_potentials__Ca_o = var_extracellular_calcium_concentration__Ca_o;
        double var_reversal_potentials__Ca_i = var_intracellular_calcium_concentration__Ca_i;
        double var_reversal_potentials__E_Ca = ((0.5 * var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Ca_o / var_reversal_potentials__Ca_i);
        double var_calcium_background_current__E_Ca = var_reversal_potentials__E_Ca;
        const double var_calcium_background_current__g_bca = 0.00025;
        double var_calcium_background_current__V = var_membrane__V;
        double var_calcium_background_current__i_b_Ca = var_calcium_background_current__g_bca * (var_calcium_background_current__V - var_calcium_background_current__E_Ca);
        double var_membrane__i_b_Ca = var_calcium_background_current__i_b_Ca;
        
        return var_membrane__i_K1+var_membrane__i_to+var_membrane__i_Kr+var_membrane__i_Ks+var_membrane__i_Ca_L_K_cyt+var_membrane__i_Ca_L_K_ds+var_membrane__i_NaK+var_membrane__i_Na+var_membrane__i_b_Na+var_membrane__i_p_Na+var_membrane__i_Ca_L_Na_cyt+var_membrane__i_Ca_L_Na_ds+var_membrane__i_NaCa_cyt+var_membrane__i_NaCa_ds+var_membrane__i_Ca_L_Ca_cyt+var_membrane__i_Ca_L_Ca_ds+var_membrane__i_b_Ca;
    }

    void ComputeResidual(const double rCurrentGuess[12], double rResidual[12])
    {
        std::vector<double>& rY = rGetStateVariables();
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -92.849333
        double var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1 = rY[1];
        // Units: dimensionless; Initial value: 1.03e-5
        double var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2 = rY[2];
        // Units: dimensionless; Initial value: 2e-7
        double var_slow_delayed_rectifier_potassium_current_xs_gate__xs = rY[3];
        // Units: dimensionless; Initial value: 0.001302
        double var_fast_sodium_current_m_gate__m = rY[4];
        // Units: dimensionless; Initial value: 0.0016203
        double var_fast_sodium_current_h_gate__h = rY[5];
        // Units: dimensionless; Initial value: 0.9944036
        double var_L_type_Ca_channel_d_gate__d = rY[6];
        // Units: dimensionless; Initial value: 0
        double var_L_type_Ca_channel_f_gate__f = rY[7];
        // Units: dimensionless; Initial value: 1
        double var_transient_outward_current_s_gate__s = rY[10];
        // Units: dimensionless; Initial value: 0.9948645
        double var_transient_outward_current_r_gate__r = rY[11];
        // Units: dimensionless; Initial value: 0
        
        double var_L_type_Ca_channel_f2_gate__f2 = rCurrentGuess[0];
        double var_L_type_Ca_channel_f2ds_gate__f2ds = rCurrentGuess[1];
        double var_calcium_release__ActFrac = rCurrentGuess[2];
        double var_calcium_release__ProdFrac = rCurrentGuess[3];
        double var_intracellular_calcium_concentration__Ca_Calmod = rCurrentGuess[4];
        double var_intracellular_calcium_concentration__Ca_Trop = rCurrentGuess[5];
        double var_intracellular_calcium_concentration__Ca_ds = rCurrentGuess[6];
        double var_intracellular_calcium_concentration__Ca_i = rCurrentGuess[7];
        double var_intracellular_calcium_concentration__Ca_rel = rCurrentGuess[8];
        double var_intracellular_calcium_concentration__Ca_up = rCurrentGuess[9];
        double var_intracellular_potassium_concentration__K_i = rCurrentGuess[10];
        double var_intracellular_sodium_concentration__Na_i = rCurrentGuess[11];
        
        const double var_membrane__R = 8314.472;
        const double var_membrane__T = 310.0;
        const double var_membrane__F = 96485.3415;
        double var_reversal_potentials__K_i = var_intracellular_potassium_concentration__K_i;
        double var_reversal_potentials__R = var_membrane__R;
        double var_reversal_potentials__T = var_membrane__T;
        double var_reversal_potentials__F = var_membrane__F;
        const double var_extracellular_potassium_concentration__K_o = 4.0;
        double var_reversal_potentials__K_o = var_extracellular_potassium_concentration__K_o;
        double var_reversal_potentials__E_K = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__K_o / var_reversal_potentials__K_i);
        double var_time_independent_potassium_current__E_K = var_reversal_potentials__E_K;
        double var_time_independent_potassium_current__K_o = var_extracellular_potassium_concentration__K_o;
        double var_time_independent_potassium_current__R = var_membrane__R;
        double var_time_independent_potassium_current__V = var_membrane__V;
        double var_time_independent_potassium_current__T = var_membrane__T;
        const double var_time_independent_potassium_current__K_mk1 = 10.0;
        const double var_time_independent_potassium_current__g_K1 = 0.5;
        double var_time_independent_potassium_current__F = var_membrane__F;
        double var_time_independent_potassium_current__i_K1 = (((var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K_o) / (var_time_independent_potassium_current__K_o + var_time_independent_potassium_current__K_mk1)) * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K)) / (1.0 + exp((((var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K) - 10.0) * var_time_independent_potassium_current__F * 1.25) / (var_time_independent_potassium_current__R * var_time_independent_potassium_current__T)));
        double var_transient_outward_current__s = var_transient_outward_current_s_gate__s;
        double var_transient_outward_current__r = var_transient_outward_current_r_gate__r;
        const double var_transient_outward_current__g_to = 0.005;
        double var_transient_outward_current__V = var_membrane__V;
        double var_transient_outward_current__E_K = var_reversal_potentials__E_K;
        const double var_transient_outward_current__g_tos = 0.0;
        double var_transient_outward_current__i_to = var_transient_outward_current__g_to * (var_transient_outward_current__g_tos + (var_transient_outward_current__s * (1.0 - var_transient_outward_current__g_tos))) * var_transient_outward_current__r * (var_transient_outward_current__V - var_transient_outward_current__E_K);
        const double var_rapid_delayed_rectifier_potassium_current__g_Kr2 = 0.0013;
        const double var_rapid_delayed_rectifier_potassium_current__g_Kr1 = 0.0021;
        double var_rapid_delayed_rectifier_potassium_current__xr1 = var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1;
        double var_rapid_delayed_rectifier_potassium_current__xr2 = var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2;
        double var_rapid_delayed_rectifier_potassium_current__V = var_membrane__V;
        double var_rapid_delayed_rectifier_potassium_current__E_K = var_reversal_potentials__E_K;
        double var_rapid_delayed_rectifier_potassium_current__i_Kr = ((((var_rapid_delayed_rectifier_potassium_current__g_Kr1 * var_rapid_delayed_rectifier_potassium_current__xr1) + (var_rapid_delayed_rectifier_potassium_current__g_Kr2 * var_rapid_delayed_rectifier_potassium_current__xr2)) * 1.0) / (1.0 + exp((var_rapid_delayed_rectifier_potassium_current__V + 9.0) / 22.4))) * (var_rapid_delayed_rectifier_potassium_current__V - var_rapid_delayed_rectifier_potassium_current__E_K);
        double var_slow_delayed_rectifier_potassium_current__xs = var_slow_delayed_rectifier_potassium_current_xs_gate__xs;
        const double var_extracellular_sodium_concentration__Na_o = 140.0;
        double var_reversal_potentials__Na_o = var_extracellular_sodium_concentration__Na_o;
        double var_reversal_potentials__Na_i = var_intracellular_sodium_concentration__Na_i;
        const double var_reversal_potentials__P_kna = 0.03;
        double var_reversal_potentials__E_Ks = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log((var_reversal_potentials__K_o + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_o)) / (var_reversal_potentials__K_i + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_i)));
        double var_slow_delayed_rectifier_potassium_current__E_Ks = var_reversal_potentials__E_Ks;
        const double var_slow_delayed_rectifier_potassium_current__g_Ks = 0.0026;
        double var_slow_delayed_rectifier_potassium_current__V = var_membrane__V;
        double var_slow_delayed_rectifier_potassium_current__i_Ks = var_slow_delayed_rectifier_potassium_current__g_Ks * pow(var_slow_delayed_rectifier_potassium_current__xs, 2.0) * (var_slow_delayed_rectifier_potassium_current__V - var_slow_delayed_rectifier_potassium_current__E_Ks);
        double var_L_type_Ca_channel__d = var_L_type_Ca_channel_d_gate__d;
        const double var_L_type_Ca_channel__FrICa = 1.0;
        double var_L_type_Ca_channel__f = var_L_type_Ca_channel_f_gate__f;
        double var_L_type_Ca_channel__K_o = var_extracellular_potassium_concentration__K_o;
        double var_L_type_Ca_channel__K_i = var_intracellular_potassium_concentration__K_i;
        double var_L_type_Ca_channel__F = var_membrane__F;
        const double var_L_type_Ca_channel__P_Ca_L = 0.1;
        double var_L_type_Ca_channel__T = var_membrane__T;
        const double var_L_type_Ca_channel__P_CaK = 0.002;
        double var_L_type_Ca_channel__V = var_membrane__V;
        double var_L_type_Ca_channel__f2 = var_L_type_Ca_channel_f2_gate__f2;
        double var_L_type_Ca_channel__R = var_membrane__R;
        double var_L_type_Ca_channel__i_Ca_L_K_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaK * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__K_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__K_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_L_type_Ca_channel__f2ds = var_L_type_Ca_channel_f2ds_gate__f2ds;
        double var_L_type_Ca_channel__i_Ca_L_K_ds = (((var_L_type_Ca_channel__FrICa * var_L_type_Ca_channel__P_CaK * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__K_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__K_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        const double var_sodium_potassium_pump__i_NaK_max = 0.7;
        double var_sodium_potassium_pump__Na_i = var_intracellular_sodium_concentration__Na_i;
        double var_sodium_potassium_pump__K_o = var_extracellular_potassium_concentration__K_o;
        const double var_sodium_potassium_pump__K_mNa = 40.0;
        const double var_sodium_potassium_pump__K_mK = 1.0;
        double var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__i_NaK_max * var_sodium_potassium_pump__K_o) / (var_sodium_potassium_pump__K_mK + var_sodium_potassium_pump__K_o)) * var_sodium_potassium_pump__Na_i) / (var_sodium_potassium_pump__K_mNa + var_sodium_potassium_pump__Na_i);
        const double var_fast_sodium_current__g_Na = 2.5;
        double var_fast_sodium_current__h = var_fast_sodium_current_h_gate__h;
        double var_fast_sodium_current__V = var_membrane__V;
        double var_reversal_potentials__E_mh = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log((var_reversal_potentials__Na_o + (0.12 * var_reversal_potentials__K_o)) / (var_reversal_potentials__Na_i + (0.12 * var_reversal_potentials__K_i)));
        double var_fast_sodium_current__E_mh = var_reversal_potentials__E_mh;
        double var_fast_sodium_current__m = var_fast_sodium_current_m_gate__m;
        double var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * (var_fast_sodium_current__V - var_fast_sodium_current__E_mh);
        double var_sodium_background_current__V = var_membrane__V;
        double var_reversal_potentials__E_Na = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Na_o / var_reversal_potentials__Na_i);
        double var_sodium_background_current__E_Na = var_reversal_potentials__E_Na;
        const double var_sodium_background_current__g_bna = 0.0006;
        double var_sodium_background_current__i_b_Na = var_sodium_background_current__g_bna * (var_sodium_background_current__V - var_sodium_background_current__E_Na);
        const double var_persistent_sodium_current__g_pna = 0.004;
        double var_persistent_sodium_current__V = var_membrane__V;
        double var_persistent_sodium_current__E_Na = var_reversal_potentials__E_Na;
        double var_persistent_sodium_current__i_p_Na = ((var_persistent_sodium_current__g_pna * 1.0) / (1.0 + exp((-(var_persistent_sodium_current__V + 52.0)) / 8.0))) * (var_persistent_sodium_current__V - var_persistent_sodium_current__E_Na);
        const double var_L_type_Ca_channel__P_CaNa = 0.01;
        double var_L_type_Ca_channel__Na_o = var_extracellular_sodium_concentration__Na_o;
        double var_L_type_Ca_channel__Na_i = var_intracellular_sodium_concentration__Na_i;
        double var_L_type_Ca_channel__i_Ca_L_Na_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaNa * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Na_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Na_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_L_type_Ca_channel__i_Ca_L_Na_ds = (((var_L_type_Ca_channel__FrICa * var_L_type_Ca_channel__P_CaNa * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Na_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Na_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_sodium_calcium_exchanger__Na_i = var_intracellular_sodium_concentration__Na_i;
        const double var_sodium_calcium_exchanger__n_NaCa = 3.0;
        const double var_sodium_calcium_exchanger__gamma = 0.5;
        double var_sodium_calcium_exchanger__F = var_membrane__F;
        double var_sodium_calcium_exchanger__Na_o = var_extracellular_sodium_concentration__Na_o;
        const double var_sodium_calcium_exchanger__FRiNaCa = 0.001;
        double var_sodium_calcium_exchanger__R = var_membrane__R;
        double var_sodium_calcium_exchanger__Ca_i = var_intracellular_calcium_concentration__Ca_i;
        double var_sodium_calcium_exchanger__T = var_membrane__T;
        double var_sodium_calcium_exchanger__V = var_membrane__V;
        const double var_sodium_calcium_exchanger__d_NaCa = 0.0;
        const double var_extracellular_calcium_concentration__Ca_o = 2.0;
        double var_sodium_calcium_exchanger__Ca_o = var_extracellular_calcium_concentration__Ca_o;
        const double var_sodium_calcium_exchanger__k_NaCa = 0.0005;
        double var_sodium_calcium_exchanger__i_NaCa_cyt = ((1.0 - var_sodium_calcium_exchanger__FRiNaCa) * var_sodium_calcium_exchanger__k_NaCa * ((exp((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_o) - (exp(((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_i))) / ((1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_sodium_calcium_exchanger__Ca_i * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_sodium_calcium_exchanger__Ca_o * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa))))) * (1.0 + (var_sodium_calcium_exchanger__Ca_i / 0.0069)));
        double var_sodium_calcium_exchanger__Ca_ds = var_intracellular_calcium_concentration__Ca_ds;
        double var_sodium_calcium_exchanger__i_NaCa_ds = (var_sodium_calcium_exchanger__FRiNaCa * var_sodium_calcium_exchanger__k_NaCa * ((exp((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_o) - (exp(((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_ds))) / ((1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_sodium_calcium_exchanger__Ca_ds * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_sodium_calcium_exchanger__Ca_o * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa))))) * (1.0 + (var_sodium_calcium_exchanger__Ca_ds / 0.0069)));
        double var_L_type_Ca_channel__Ca_i = var_intracellular_calcium_concentration__Ca_i;
        double var_L_type_Ca_channel__Ca_o = var_extracellular_calcium_concentration__Ca_o;
        double var_L_type_Ca_channel__i_Ca_L_Ca_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * 4.0 * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Ca_i * exp((100.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Ca_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_L_type_Ca_channel__i_Ca_L_Ca_ds = (((var_L_type_Ca_channel__FrICa * 4.0 * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Ca_i * exp((100.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Ca_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_reversal_potentials__Ca_o = var_extracellular_calcium_concentration__Ca_o;
        double var_reversal_potentials__Ca_i = var_intracellular_calcium_concentration__Ca_i;
        double var_reversal_potentials__E_Ca = ((0.5 * var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Ca_o / var_reversal_potentials__Ca_i);
        double var_calcium_background_current__E_Ca = var_reversal_potentials__E_Ca;
        const double var_calcium_background_current__g_bca = 0.00025;
        double var_calcium_background_current__V = var_membrane__V;
        double var_calcium_background_current__i_b_Ca = var_calcium_background_current__g_bca * (var_calcium_background_current__V - var_calcium_background_current__E_Ca);
        double var_L_type_Ca_channel__Ca_ds = var_intracellular_calcium_concentration__Ca_ds;
        const double var_L_type_Ca_channel__Km_f2 = 100000.0;
        const double var_L_type_Ca_channel__Km_f2ds = 0.001;
        const double var_L_type_Ca_channel__R_decay = 20.0;
        double var_L_type_Ca_channel_f2_gate__Km_f2 = var_L_type_Ca_channel__Km_f2;
        double var_L_type_Ca_channel_f2_gate__Ca_i = var_L_type_Ca_channel__Ca_i;
        double var_L_type_Ca_channel_f2ds_gate__Km_f2ds = var_L_type_Ca_channel__Km_f2ds;
        double var_L_type_Ca_channel_f2ds_gate__R_decay = var_L_type_Ca_channel__R_decay;
        double var_L_type_Ca_channel_f2ds_gate__Ca_ds = var_L_type_Ca_channel__Ca_ds;
        double var_sarcoplasmic_reticulum_calcium_pump__Ca_i = var_intracellular_calcium_concentration__Ca_i;
        double var_sarcoplasmic_reticulum_calcium_pump__Ca_up = var_intracellular_calcium_concentration__Ca_up;
        const double var_sarcoplasmic_reticulum_calcium_pump__alpha_up = 0.4;
        const double var_sarcoplasmic_reticulum_calcium_pump__beta_up = 0.03;
        const double var_sarcoplasmic_reticulum_calcium_pump__K_srca = 0.5;
        const double var_sarcoplasmic_reticulum_calcium_pump__K_xcs = 0.4;
        const double var_sarcoplasmic_reticulum_calcium_pump__K_cyca = 0.0003;
        double var_sarcoplasmic_reticulum_calcium_pump__K_1 = (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca;
        double var_sarcoplasmic_reticulum_calcium_pump__K_2 = var_sarcoplasmic_reticulum_calcium_pump__Ca_i + (var_sarcoplasmic_reticulum_calcium_pump__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_1) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca;
        double var_sarcoplasmic_reticulum_calcium_pump__i_up = ((var_sarcoplasmic_reticulum_calcium_pump__Ca_i / var_sarcoplasmic_reticulum_calcium_pump__K_2) * var_sarcoplasmic_reticulum_calcium_pump__alpha_up) - (((var_sarcoplasmic_reticulum_calcium_pump__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_1) / var_sarcoplasmic_reticulum_calcium_pump__K_2) * var_sarcoplasmic_reticulum_calcium_pump__beta_up);
        double var_calcium_translocation__Ca_rel = var_intracellular_calcium_concentration__Ca_rel;
        double var_calcium_translocation__Ca_up = var_intracellular_calcium_concentration__Ca_up;
        double var_calcium_translocation__i_trans = 50.0 * (var_calcium_translocation__Ca_up - var_calcium_translocation__Ca_rel);
        const double var_calcium_release__K_m_rel = 250.0;
        const double var_calcium_release__K_leak_rate = 0.05;
        double var_calcium_release__Ca_rel = var_intracellular_calcium_concentration__Ca_rel;
        double var_calcium_release__i_rel = ((pow(var_calcium_release__ActFrac / (var_calcium_release__ActFrac + 0.25), 2.0) * var_calcium_release__K_m_rel) + var_calcium_release__K_leak_rate) * var_calcium_release__Ca_rel;
        double var_calcium_release__V = var_membrane__V;
        double var_calcium_release__VoltDep = exp(0.08 * (var_calcium_release__V - 40.0));
        const double var_calcium_release__K_m_Ca_cyt = 0.0005;
        double var_calcium_release__Ca_i = var_intracellular_calcium_concentration__Ca_i;
        double var_calcium_release__CaiReg = var_calcium_release__Ca_i / (var_calcium_release__Ca_i + var_calcium_release__K_m_Ca_cyt);
        double var_calcium_release__Ca_ds = var_intracellular_calcium_concentration__Ca_ds;
        const double var_calcium_release__K_m_Ca_ds = 0.01;
        double var_calcium_release__CadsReg = var_calcium_release__Ca_ds / (var_calcium_release__Ca_ds + var_calcium_release__K_m_Ca_ds);
        double var_calcium_release__RegBindSite = var_calcium_release__CaiReg + ((1.0 - var_calcium_release__CaiReg) * var_calcium_release__CadsReg);
        double var_calcium_release__ActRate = (0.0 * var_calcium_release__VoltDep) + (500.0 * pow(var_calcium_release__RegBindSite, 2.0));
        double var_calcium_release__InactRate = 60.0 + (500.0 * pow(var_calcium_release__RegBindSite, 2.0));
        double var_calcium_release__PrecFrac = (1.0 - var_calcium_release__ActFrac) - var_calcium_release__ProdFrac;
        double var_calcium_release__SpeedRel = (var_calcium_release__V < (-50.0)) ? 5.0 : 1.0;
        const double var_intracellular_calcium_concentration__V_up_ratio = 0.01;
        const double var_intracellular_calcium_concentration__V_rel_ratio = 0.1;
        const double var_intracellular_calcium_concentration__V_e_ratio = 0.4;
        double var_intracellular_calcium_concentration__V_i_ratio = ((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio;
        const double var_intracellular_calcium_concentration__radius = 0.012;
        const double var_intracellular_calcium_concentration__length = 0.074;
        double var_intracellular_calcium_concentration__V_Cell = 3.141592654 * pow(var_intracellular_calcium_concentration__radius, 2.0) * var_intracellular_calcium_concentration__length;
        double var_intracellular_calcium_concentration__V_i = var_intracellular_calcium_concentration__V_Cell * var_intracellular_calcium_concentration__V_i_ratio;
        double var_intracellular_sodium_concentration__V_i = var_intracellular_calcium_concentration__V_i;
        double var_intracellular_sodium_concentration__F = var_membrane__F;
        double var_intracellular_sodium_concentration__i_Na = var_fast_sodium_current__i_Na;
        double var_intracellular_sodium_concentration__i_b_Na = var_sodium_background_current__i_b_Na;
        double var_intracellular_sodium_concentration__i_p_Na = var_persistent_sodium_current__i_p_Na;
        double var_intracellular_sodium_concentration__i_Ca_L_Na_cyt = var_L_type_Ca_channel__i_Ca_L_Na_cyt;
        double var_intracellular_sodium_concentration__i_Ca_L_Na_ds = var_L_type_Ca_channel__i_Ca_L_Na_ds;
        double var_intracellular_sodium_concentration__i_NaK = var_sodium_potassium_pump__i_NaK;
        double var_intracellular_sodium_concentration__i_NaCa_cyt = var_sodium_calcium_exchanger__i_NaCa_cyt;
        double var_intracellular_potassium_concentration__V_i = var_intracellular_calcium_concentration__V_i;
        double var_intracellular_potassium_concentration__i_K1 = var_time_independent_potassium_current__i_K1;
        double var_intracellular_potassium_concentration__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
        double var_intracellular_potassium_concentration__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
        double var_intracellular_potassium_concentration__i_Ca_L_K_cyt = var_L_type_Ca_channel__i_Ca_L_K_cyt;
        double var_intracellular_potassium_concentration__i_Ca_L_K_ds = var_L_type_Ca_channel__i_Ca_L_K_ds;
        double var_intracellular_potassium_concentration__i_to = var_transient_outward_current__i_to;
        double var_intracellular_potassium_concentration__i_NaK = var_sodium_potassium_pump__i_NaK;
        double var_intracellular_potassium_concentration__F = var_membrane__F;
        const double var_intracellular_calcium_concentration__Calmod = 0.02;
        const double var_intracellular_calcium_concentration__Trop = 0.05;
        const double var_intracellular_calcium_concentration__alpha_Calmod = 100000.0;
        const double var_intracellular_calcium_concentration__beta_Calmod = 50.0;
        const double var_intracellular_calcium_concentration__alpha_Trop = 100000.0;
        const double var_intracellular_calcium_concentration__beta_Trop = 200.0;
        const double var_intracellular_calcium_concentration__V_ds_ratio = 0.1;
        const double var_intracellular_calcium_concentration__Kdecay = 10.0;
        double var_intracellular_calcium_concentration__i_up = var_sarcoplasmic_reticulum_calcium_pump__i_up;
        double var_intracellular_calcium_concentration__i_trans = var_calcium_translocation__i_trans;
        double var_intracellular_calcium_concentration__i_rel = var_calcium_release__i_rel;
        double var_intracellular_calcium_concentration__i_NaCa_cyt = var_sodium_calcium_exchanger__i_NaCa_cyt;
        double var_intracellular_calcium_concentration__i_NaCa_ds = var_sodium_calcium_exchanger__i_NaCa_ds;
        double var_intracellular_calcium_concentration__i_Ca_L_Ca_cyt = var_L_type_Ca_channel__i_Ca_L_Ca_cyt;
        double var_intracellular_calcium_concentration__i_Ca_L_Ca_ds = var_L_type_Ca_channel__i_Ca_L_Ca_ds;
        double var_intracellular_calcium_concentration__i_b_Ca = var_calcium_background_current__i_b_Ca;
        double var_intracellular_calcium_concentration__F = var_membrane__F;
        double d_dt_L_type_Ca_channel_f2_gate__f2 = 1.0 - (1.0 * ((var_L_type_Ca_channel_f2_gate__Ca_i / (var_L_type_Ca_channel_f2_gate__Km_f2 + var_L_type_Ca_channel_f2_gate__Ca_i)) + var_L_type_Ca_channel_f2_gate__f2));
        double d_dt_L_type_Ca_channel_f2ds_gate__f2ds = var_L_type_Ca_channel_f2ds_gate__R_decay * (1.0 - ((var_L_type_Ca_channel_f2ds_gate__Ca_ds / (var_L_type_Ca_channel_f2ds_gate__Km_f2ds + var_L_type_Ca_channel_f2ds_gate__Ca_ds)) + var_L_type_Ca_channel_f2ds_gate__f2ds));
        double d_dt_calcium_release__ActFrac = (var_calcium_release__PrecFrac * var_calcium_release__SpeedRel * var_calcium_release__ActRate) - (var_calcium_release__ActFrac * var_calcium_release__SpeedRel * var_calcium_release__InactRate);
        double d_dt_calcium_release__ProdFrac = (var_calcium_release__ActFrac * var_calcium_release__SpeedRel * var_calcium_release__InactRate) - (var_calcium_release__SpeedRel * 1.0 * var_calcium_release__ProdFrac);
        double d_dt_intracellular_sodium_concentration__Na_i = ((-1.0) / (1.0 * var_intracellular_sodium_concentration__V_i * var_intracellular_sodium_concentration__F)) * (var_intracellular_sodium_concentration__i_Na + var_intracellular_sodium_concentration__i_p_Na + var_intracellular_sodium_concentration__i_b_Na + (3.0 * var_intracellular_sodium_concentration__i_NaK) + (3.0 * var_intracellular_sodium_concentration__i_NaCa_cyt) + var_intracellular_sodium_concentration__i_Ca_L_Na_cyt + var_intracellular_sodium_concentration__i_Ca_L_Na_ds);
        double d_dt_intracellular_potassium_concentration__K_i = ((-1.0) / (1.0 * var_intracellular_potassium_concentration__V_i * var_intracellular_potassium_concentration__F)) * ((var_intracellular_potassium_concentration__i_K1 + var_intracellular_potassium_concentration__i_Kr + var_intracellular_potassium_concentration__i_Ks + var_intracellular_potassium_concentration__i_Ca_L_K_cyt + var_intracellular_potassium_concentration__i_Ca_L_K_ds + var_intracellular_potassium_concentration__i_to) - (2.0 * var_intracellular_potassium_concentration__i_NaK));
        double d_dt_intracellular_calcium_concentration__Ca_Trop = (var_intracellular_calcium_concentration__alpha_Trop * var_intracellular_calcium_concentration__Ca_i * (var_intracellular_calcium_concentration__Trop - var_intracellular_calcium_concentration__Ca_Trop)) - (var_intracellular_calcium_concentration__beta_Trop * var_intracellular_calcium_concentration__Ca_Trop);
        double d_dt_intracellular_calcium_concentration__Ca_Calmod = (var_intracellular_calcium_concentration__alpha_Calmod * var_intracellular_calcium_concentration__Ca_i * (var_intracellular_calcium_concentration__Calmod - var_intracellular_calcium_concentration__Ca_Calmod)) - (var_intracellular_calcium_concentration__beta_Calmod * var_intracellular_calcium_concentration__Ca_Calmod);
        double d_dt_intracellular_calcium_concentration__Ca_i = ((((((-1.0) / (2.0 * 1.0 * var_intracellular_calcium_concentration__V_i * var_intracellular_calcium_concentration__F)) * (((var_intracellular_calcium_concentration__i_Ca_L_Ca_cyt + var_intracellular_calcium_concentration__i_b_Ca) - (2.0 * var_intracellular_calcium_concentration__i_NaCa_cyt)) - (2.0 * var_intracellular_calcium_concentration__i_NaCa_ds))) + (var_intracellular_calcium_concentration__Ca_ds * var_intracellular_calcium_concentration__V_ds_ratio * var_intracellular_calcium_concentration__Kdecay) + ((var_intracellular_calcium_concentration__i_rel * var_intracellular_calcium_concentration__V_rel_ratio) / var_intracellular_calcium_concentration__V_i_ratio)) - d_dt_intracellular_calcium_concentration__Ca_Calmod) - d_dt_intracellular_calcium_concentration__Ca_Trop) - var_intracellular_calcium_concentration__i_up;
        double d_dt_intracellular_calcium_concentration__Ca_ds = (((-1.0) * var_intracellular_calcium_concentration__i_Ca_L_Ca_ds) / (2.0 * 1.0 * var_intracellular_calcium_concentration__V_ds_ratio * var_intracellular_calcium_concentration__V_i * var_intracellular_calcium_concentration__F)) - (var_intracellular_calcium_concentration__Ca_ds * var_intracellular_calcium_concentration__Kdecay);
        double d_dt_intracellular_calcium_concentration__Ca_up = ((var_intracellular_calcium_concentration__V_i_ratio / var_intracellular_calcium_concentration__V_up_ratio) * var_intracellular_calcium_concentration__i_up) - var_intracellular_calcium_concentration__i_trans;
        double d_dt_intracellular_calcium_concentration__Ca_rel = ((var_intracellular_calcium_concentration__V_up_ratio / var_intracellular_calcium_concentration__V_rel_ratio) * var_intracellular_calcium_concentration__i_trans) - var_intracellular_calcium_concentration__i_rel;
        
        rResidual[0] = rCurrentGuess[0] - rY[8] - mDt*0.001*d_dt_L_type_Ca_channel_f2_gate__f2;
        rResidual[1] = rCurrentGuess[1] - rY[9] - mDt*0.001*d_dt_L_type_Ca_channel_f2ds_gate__f2ds;
        rResidual[2] = rCurrentGuess[2] - rY[12] - mDt*0.001*d_dt_calcium_release__ActFrac;
        rResidual[3] = rCurrentGuess[3] - rY[13] - mDt*0.001*d_dt_calcium_release__ProdFrac;
        rResidual[11] = rCurrentGuess[11] - rY[14] - mDt*0.001*d_dt_intracellular_sodium_concentration__Na_i;
        rResidual[10] = rCurrentGuess[10] - rY[15] - mDt*0.001*d_dt_intracellular_potassium_concentration__K_i;
        rResidual[7] = rCurrentGuess[7] - rY[16] - mDt*0.001*d_dt_intracellular_calcium_concentration__Ca_i;
        rResidual[6] = rCurrentGuess[6] - rY[17] - mDt*0.001*d_dt_intracellular_calcium_concentration__Ca_ds;
        rResidual[9] = rCurrentGuess[9] - rY[18] - mDt*0.001*d_dt_intracellular_calcium_concentration__Ca_up;
        rResidual[8] = rCurrentGuess[8] - rY[19] - mDt*0.001*d_dt_intracellular_calcium_concentration__Ca_rel;
        rResidual[4] = rCurrentGuess[4] - rY[20] - mDt*0.001*d_dt_intracellular_calcium_concentration__Ca_Calmod;
        rResidual[5] = rCurrentGuess[5] - rY[21] - mDt*0.001*d_dt_intracellular_calcium_concentration__Ca_Trop;
    }

    void ComputeJacobian(const double rCurrentGuess[12], double rJacobian[12][12])
    {
        std::vector<double>& rY = rGetStateVariables();
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -92.849333
        double var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1 = rY[1];
        // Units: dimensionless; Initial value: 1.03e-5
        double var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2 = rY[2];
        // Units: dimensionless; Initial value: 2e-7
        double var_slow_delayed_rectifier_potassium_current_xs_gate__xs = rY[3];
        // Units: dimensionless; Initial value: 0.001302
        double var_fast_sodium_current_m_gate__m = rY[4];
        // Units: dimensionless; Initial value: 0.0016203
        double var_fast_sodium_current_h_gate__h = rY[5];
        // Units: dimensionless; Initial value: 0.9944036
        double var_L_type_Ca_channel_d_gate__d = rY[6];
        // Units: dimensionless; Initial value: 0
        double var_L_type_Ca_channel_f_gate__f = rY[7];
        // Units: dimensionless; Initial value: 1
        double var_transient_outward_current_s_gate__s = rY[10];
        // Units: dimensionless; Initial value: 0.9948645
        double var_transient_outward_current_r_gate__r = rY[11];
        // Units: dimensionless; Initial value: 0
        
        double var_L_type_Ca_channel_f2_gate__f2 = rCurrentGuess[0];
        double var_L_type_Ca_channel_f2ds_gate__f2ds = rCurrentGuess[1];
        double var_calcium_release__ActFrac = rCurrentGuess[2];
        double var_calcium_release__ProdFrac = rCurrentGuess[3];
        double var_intracellular_calcium_concentration__Ca_Calmod = rCurrentGuess[4];
        double var_intracellular_calcium_concentration__Ca_Trop = rCurrentGuess[5];
        double var_intracellular_calcium_concentration__Ca_ds = rCurrentGuess[6];
        double var_intracellular_calcium_concentration__Ca_i = rCurrentGuess[7];
        double var_intracellular_calcium_concentration__Ca_rel = rCurrentGuess[8];
        double var_intracellular_calcium_concentration__Ca_up = rCurrentGuess[9];
        double var_intracellular_potassium_concentration__K_i = rCurrentGuess[10];
        double var_intracellular_sodium_concentration__Na_i = rCurrentGuess[11];
        
        const double var_membrane__R = 8314.472;
        const double var_membrane__T = 310.0;
        const double var_membrane__F = 96485.3415;
        const double var_extracellular_potassium_concentration__K_o = 4.0;
        const double var_time_independent_potassium_current__K_mk1 = 10.0;
        const double var_time_independent_potassium_current__g_K1 = 0.5;
        const double var_transient_outward_current__g_to = 0.005;
        const double var_transient_outward_current__g_tos = 0.0;
        const double var_rapid_delayed_rectifier_potassium_current__g_Kr2 = 0.0013;
        const double var_rapid_delayed_rectifier_potassium_current__g_Kr1 = 0.0021;
        const double var_extracellular_sodium_concentration__Na_o = 140.0;
        const double var_reversal_potentials__P_kna = 0.03;
        const double var_slow_delayed_rectifier_potassium_current__g_Ks = 0.0026;
        const double var_L_type_Ca_channel__FrICa = 1.0;
        const double var_L_type_Ca_channel__P_Ca_L = 0.1;
        const double var_L_type_Ca_channel__P_CaK = 0.002;
        const double var_sodium_potassium_pump__i_NaK_max = 0.7;
        const double var_sodium_potassium_pump__K_mNa = 40.0;
        const double var_sodium_potassium_pump__K_mK = 1.0;
        const double var_fast_sodium_current__g_Na = 2.5;
        const double var_sodium_background_current__g_bna = 0.0006;
        const double var_persistent_sodium_current__g_pna = 0.004;
        const double var_L_type_Ca_channel__P_CaNa = 0.01;
        const double var_sodium_calcium_exchanger__n_NaCa = 3.0;
        const double var_sodium_calcium_exchanger__gamma = 0.5;
        const double var_sodium_calcium_exchanger__FRiNaCa = 0.001;
        const double var_sodium_calcium_exchanger__d_NaCa = 0.0;
        const double var_extracellular_calcium_concentration__Ca_o = 2.0;
        const double var_sodium_calcium_exchanger__k_NaCa = 0.0005;
        const double var_calcium_background_current__g_bca = 0.00025;
        const double var_L_type_Ca_channel__Km_f2 = 100000.0;
        const double var_L_type_Ca_channel__Km_f2ds = 0.001;
        const double var_L_type_Ca_channel__R_decay = 20.0;
        const double var_sarcoplasmic_reticulum_calcium_pump__alpha_up = 0.4;
        const double var_sarcoplasmic_reticulum_calcium_pump__beta_up = 0.03;
        const double var_sarcoplasmic_reticulum_calcium_pump__K_srca = 0.5;
        const double var_sarcoplasmic_reticulum_calcium_pump__K_xcs = 0.4;
        const double var_sarcoplasmic_reticulum_calcium_pump__K_cyca = 0.0003;
        const double var_calcium_release__K_m_rel = 250.0;
        const double var_calcium_release__K_leak_rate = 0.05;
        const double var_calcium_release__K_m_Ca_cyt = 0.0005;
        const double var_calcium_release__K_m_Ca_ds = 0.01;
        const double var_intracellular_calcium_concentration__V_up_ratio = 0.01;
        const double var_intracellular_calcium_concentration__V_rel_ratio = 0.1;
        const double var_intracellular_calcium_concentration__V_e_ratio = 0.4;
        const double var_intracellular_calcium_concentration__radius = 0.012;
        const double var_intracellular_calcium_concentration__length = 0.074;
        const double var_intracellular_calcium_concentration__Calmod = 0.02;
        const double var_intracellular_calcium_concentration__Trop = 0.05;
        const double var_intracellular_calcium_concentration__alpha_Calmod = 100000.0;
        const double var_intracellular_calcium_concentration__beta_Calmod = 50.0;
        const double var_intracellular_calcium_concentration__alpha_Trop = 100000.0;
        const double var_intracellular_calcium_concentration__beta_Trop = 200.0;
        const double var_intracellular_calcium_concentration__V_ds_ratio = 0.1;
        const double var_intracellular_calcium_concentration__Kdecay = 10.0;
        
        rJacobian[0][0] = 1.0 + 0.001*mDt;
        rJacobian[0][1] = 0.0;
        rJacobian[0][2] = 0.0;
        rJacobian[0][3] = 0.0;
        rJacobian[0][4] = 0.0;
        rJacobian[0][5] = 0.0;
        rJacobian[0][6] = 0.0;
        rJacobian[0][7] = (-0.001*mDt) * (((-1.0) / (var_L_type_Ca_channel__Km_f2 + var_intracellular_calcium_concentration__Ca_i)) + (var_intracellular_calcium_concentration__Ca_i / pow(var_L_type_Ca_channel__Km_f2 + var_intracellular_calcium_concentration__Ca_i, 2.0)));
        rJacobian[0][8] = 0.0;
        rJacobian[0][9] = 0.0;
        rJacobian[0][10] = 0.0;
        rJacobian[0][11] = 0.0;
        rJacobian[1][0] = 0.0;
        rJacobian[1][1] = 1.0 + (0.001*mDt * var_L_type_Ca_channel__R_decay);
        rJacobian[1][2] = 0.0;
        rJacobian[1][3] = 0.0;
        rJacobian[1][4] = 0.0;
        rJacobian[1][5] = 0.0;
        rJacobian[1][6] = ((-0.001*mDt) * var_L_type_Ca_channel__R_decay) * (((-1.0) / (var_L_type_Ca_channel__Km_f2ds + var_intracellular_calcium_concentration__Ca_ds)) + (var_intracellular_calcium_concentration__Ca_ds / pow(var_L_type_Ca_channel__Km_f2ds + var_intracellular_calcium_concentration__Ca_ds, 2.0)));
        rJacobian[1][7] = 0.0;
        rJacobian[1][8] = 0.0;
        rJacobian[1][9] = 0.0;
        rJacobian[1][10] = 0.0;
        rJacobian[1][11] = 0.0;
        rJacobian[2][0] = 0.0;
        rJacobian[2][1] = 0.0;
        rJacobian[2][2] = 1.0 - (0.001*mDt * ((((-500.0) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * pow((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)), 2.0)) - (((var_membrane__V < (-50.0)) ? 5.0 : 1.0) * (60.0 + (500.0 * pow((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)), 2.0))))));
        rJacobian[2][3] = ((500.0 * 0.001*mDt) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * pow((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)), 2.0);
        rJacobian[2][4] = 0.0;
        rJacobian[2][5] = 0.0;
        rJacobian[2][6] = (-0.001*mDt) * (((((1000.0 * ((1.0 - var_calcium_release__ActFrac) - var_calcium_release__ProdFrac)) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * ((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)))) * (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)) - (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) * var_intracellular_calcium_concentration__Ca_ds) / pow(var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds, 2.0)))) - ((((1000.0 * var_calcium_release__ActFrac) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * ((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)))) * (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)) - (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) * var_intracellular_calcium_concentration__Ca_ds) / pow(var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds, 2.0)))));
        rJacobian[2][7] = (-0.001*mDt) * (((((1000.0 * ((1.0 - var_calcium_release__ActFrac) - var_calcium_release__ProdFrac)) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * ((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)))) * (((1.0 / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) - (var_intracellular_calcium_concentration__Ca_i / pow(var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt, 2.0))) + (((((-1.0) / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) + (var_intracellular_calcium_concentration__Ca_i / pow(var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt, 2.0))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)))) - ((((1000.0 * var_calcium_release__ActFrac) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * ((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)))) * (((1.0 / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) - (var_intracellular_calcium_concentration__Ca_i / pow(var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt, 2.0))) + (((((-1.0) / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) + (var_intracellular_calcium_concentration__Ca_i / pow(var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt, 2.0))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)))));
        rJacobian[2][8] = 0.0;
        rJacobian[2][9] = 0.0;
        rJacobian[2][10] = 0.0;
        rJacobian[2][11] = 0.0;
        rJacobian[3][0] = 0.0;
        rJacobian[3][1] = 0.0;
        rJacobian[3][2] = ((-0.001*mDt) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * (60.0 + (500.0 * pow((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)), 2.0)));
        rJacobian[3][3] = 1.0 + (0.001*mDt * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0));
        rJacobian[3][4] = 0.0;
        rJacobian[3][5] = 0.0;
        rJacobian[3][6] = (((((-1000.0) * 0.001*mDt) * var_calcium_release__ActFrac) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * ((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)))) * (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)) - (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) * var_intracellular_calcium_concentration__Ca_ds) / pow(var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds, 2.0)));
        rJacobian[3][7] = (((((-1000.0) * 0.001*mDt) * var_calcium_release__ActFrac) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * ((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)))) * (((1.0 / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) - (var_intracellular_calcium_concentration__Ca_i / pow(var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt, 2.0))) + (((((-1.0) / (var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt)) + (var_intracellular_calcium_concentration__Ca_i / pow(var_intracellular_calcium_concentration__Ca_i + var_calcium_release__K_m_Ca_cyt, 2.0))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + var_calcium_release__K_m_Ca_ds)));
        rJacobian[3][8] = 0.0;
        rJacobian[3][9] = 0.0;
        rJacobian[3][10] = 0.0;
        rJacobian[3][11] = 0.0;
        rJacobian[4][0] = 0.0;
        rJacobian[4][1] = 0.0;
        rJacobian[4][2] = 0.0;
        rJacobian[4][3] = 0.0;
        rJacobian[4][4] = 1.0 - (0.001*mDt * (((-var_intracellular_calcium_concentration__alpha_Calmod) * var_intracellular_calcium_concentration__Ca_i) - var_intracellular_calcium_concentration__beta_Calmod));
        rJacobian[4][5] = 0.0;
        rJacobian[4][6] = 0.0;
        rJacobian[4][7] = ((-0.001*mDt) * var_intracellular_calcium_concentration__alpha_Calmod) * (var_intracellular_calcium_concentration__Calmod - var_intracellular_calcium_concentration__Ca_Calmod);
        rJacobian[4][8] = 0.0;
        rJacobian[4][9] = 0.0;
        rJacobian[4][10] = 0.0;
        rJacobian[4][11] = 0.0;
        rJacobian[5][0] = 0.0;
        rJacobian[5][1] = 0.0;
        rJacobian[5][2] = 0.0;
        rJacobian[5][3] = 0.0;
        rJacobian[5][4] = 0.0;
        rJacobian[5][5] = 1.0 - (0.001*mDt * (((-var_intracellular_calcium_concentration__alpha_Trop) * var_intracellular_calcium_concentration__Ca_i) - var_intracellular_calcium_concentration__beta_Trop));
        rJacobian[5][6] = 0.0;
        rJacobian[5][7] = ((-0.001*mDt) * var_intracellular_calcium_concentration__alpha_Trop) * (var_intracellular_calcium_concentration__Trop - var_intracellular_calcium_concentration__Ca_Trop);
        rJacobian[5][8] = 0.0;
        rJacobian[5][9] = 0.0;
        rJacobian[5][10] = 0.0;
        rJacobian[5][11] = 0.0;
        rJacobian[6][0] = 0.0;
        rJacobian[6][1] = (((((((((((((0.6366197725 * 0.001*mDt) * var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_Ca_L) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * (var_membrane__V - 50.0)) / var_membrane__R) / var_membrane__T) / (1.0 - exp((((2.0 * ((-var_membrane__V) + 50.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T))) * ((var_intracellular_calcium_concentration__Ca_i * exp(((100.0 * var_membrane__F) / var_membrane__R) / var_membrane__T)) - (var_extracellular_calcium_concentration__Ca_o * exp((((2.0 * ((-var_membrane__V) + 50.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T)))) / var_intracellular_calcium_concentration__V_ds_ratio) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio);
        rJacobian[6][2] = 0.0;
        rJacobian[6][3] = 0.0;
        rJacobian[6][4] = 0.0;
        rJacobian[6][5] = 0.0;
        rJacobian[6][6] = 1.0 + (0.001*mDt * var_intracellular_calcium_concentration__Kdecay);
        rJacobian[6][7] = ((((((((((((((0.6366197725 * 0.001*mDt) * var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_Ca_L) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * var_L_type_Ca_channel_f2ds_gate__f2ds) * (var_membrane__V - 50.0)) / var_membrane__R) / var_membrane__T) / (1.0 - exp((((2.0 * ((-var_membrane__V) + 50.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T))) * exp(((100.0 * var_membrane__F) / var_membrane__R) / var_membrane__T)) / var_intracellular_calcium_concentration__V_ds_ratio) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio);
        rJacobian[6][8] = 0.0;
        rJacobian[6][9] = 0.0;
        rJacobian[6][10] = 0.0;
        rJacobian[6][11] = 0.0;
        rJacobian[7][0] = ((((((((((((0.636619772 * 0.001*mDt) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) * (1.0 - var_L_type_Ca_channel__FrICa)) * var_L_type_Ca_channel__P_Ca_L) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * (var_membrane__V - 50.0)) / var_membrane__R) / var_membrane__T) / (1.0 - exp((((2.0 * ((-var_membrane__V) + 50.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T))) * ((var_intracellular_calcium_concentration__Ca_i * exp(((100.0 * var_membrane__F) / var_membrane__R) / var_membrane__T)) - (var_extracellular_calcium_concentration__Ca_o * exp((((2.0 * ((-var_membrane__V) + 50.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T)));
        rJacobian[7][1] = 0.0;
        rJacobian[7][2] = ((((-0.001*mDt) * ((((2.0 * var_calcium_release__ActFrac) / pow(var_calcium_release__ActFrac + 0.25, 2.0)) * var_calcium_release__K_m_rel) - (((2.0 * pow(var_calcium_release__ActFrac, 2.0)) / pow(var_calcium_release__ActFrac + 0.25, 3.0)) * var_calcium_release__K_m_rel))) * var_intracellular_calcium_concentration__Ca_rel) * var_intracellular_calcium_concentration__V_rel_ratio) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio);
        rJacobian[7][3] = 0.0;
        rJacobian[7][4] = (-0.001*mDt) * ((var_intracellular_calcium_concentration__alpha_Calmod * var_intracellular_calcium_concentration__Ca_i) + var_intracellular_calcium_concentration__beta_Calmod);
        rJacobian[7][5] = (-0.001*mDt) * ((var_intracellular_calcium_concentration__alpha_Trop * var_intracellular_calcium_concentration__Ca_i) + var_intracellular_calcium_concentration__beta_Trop);
        rJacobian[7][6] = (-0.001*mDt) * (((((((-0.159154943) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) / var_membrane__F) * ((((((((2.0 * var_sodium_calcium_exchanger__FRiNaCa) * var_sodium_calcium_exchanger__k_NaCa) * exp((((((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T)) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) / (1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_ds * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))))) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_ds))) + (((((((2.0 * var_sodium_calcium_exchanger__FRiNaCa) * var_sodium_calcium_exchanger__k_NaCa) * (((exp(((((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_extracellular_calcium_concentration__Ca_o) - ((exp((((((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) * var_intracellular_calcium_concentration__Ca_ds))) / pow(1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_ds * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))), 2.0)) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_ds))) * var_sodium_calcium_exchanger__d_NaCa) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa))) + (((((289.8550724 * var_sodium_calcium_exchanger__FRiNaCa) * var_sodium_calcium_exchanger__k_NaCa) * (((exp(((((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_extracellular_calcium_concentration__Ca_o) - ((exp((((((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) * var_intracellular_calcium_concentration__Ca_ds))) / (1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_ds * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))))) / pow(1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_ds), 2.0)))) + (var_intracellular_calcium_concentration__V_ds_ratio * var_intracellular_calcium_concentration__Kdecay));
        rJacobian[7][7] = 1.0 - (0.001*mDt * (((((((((((-0.159154943) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) / var_membrane__F) * (((((((((((((((4.0 * (1.0 - var_L_type_Ca_channel__FrICa)) * var_L_type_Ca_channel__P_Ca_L) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * var_L_type_Ca_channel_f2_gate__f2) * (var_membrane__V - 50.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T) / (1.0 - exp((((2.0 * ((-var_membrane__V) + 50.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T))) * exp(((100.0 * var_membrane__F) / var_membrane__R) / var_membrane__T)) + (((((0.5 * var_calcium_background_current__g_bca) * var_membrane__R) * var_membrane__T) / var_membrane__F) / var_intracellular_calcium_concentration__Ca_i)) + ((((((2.0 * (1.0 - var_sodium_calcium_exchanger__FRiNaCa)) * var_sodium_calcium_exchanger__k_NaCa) * exp((((((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T)) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) / (1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_i * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))))) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i)))) + (((((((2.0 * (1.0 - var_sodium_calcium_exchanger__FRiNaCa)) * var_sodium_calcium_exchanger__k_NaCa) * (((exp(((((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_extracellular_calcium_concentration__Ca_o) - ((exp((((((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) * var_intracellular_calcium_concentration__Ca_i))) / pow(1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_i * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))), 2.0)) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i))) * var_sodium_calcium_exchanger__d_NaCa) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa))) + (((((289.8550724 * (1.0 - var_sodium_calcium_exchanger__FRiNaCa)) * var_sodium_calcium_exchanger__k_NaCa) * (((exp(((((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_extracellular_calcium_concentration__Ca_o) - ((exp((((((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) * var_intracellular_calcium_concentration__Ca_i))) / (1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_i * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))))) / pow(1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i), 2.0)))) - (var_intracellular_calcium_concentration__alpha_Calmod * (var_intracellular_calcium_concentration__Calmod - var_intracellular_calcium_concentration__Ca_Calmod))) - (var_intracellular_calcium_concentration__alpha_Trop * (var_intracellular_calcium_concentration__Trop - var_intracellular_calcium_concentration__Ca_Trop))) - ((1.0 / (((var_intracellular_calcium_concentration__Ca_i + (((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca)) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs)) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca)) * var_sarcoplasmic_reticulum_calcium_pump__alpha_up)) + ((var_intracellular_calcium_concentration__Ca_i / pow(((var_intracellular_calcium_concentration__Ca_i + (((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca)) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs)) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca, 2.0)) * var_sarcoplasmic_reticulum_calcium_pump__alpha_up)) - (((((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca) / pow(((var_intracellular_calcium_concentration__Ca_i + (((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca)) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs)) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca, 2.0)) * var_sarcoplasmic_reticulum_calcium_pump__beta_up)));
        rJacobian[7][8] = (((-0.001*mDt) * (((pow(var_calcium_release__ActFrac, 2.0) / pow(var_calcium_release__ActFrac + 0.25, 2.0)) * var_calcium_release__K_m_rel) + var_calcium_release__K_leak_rate)) * var_intracellular_calcium_concentration__V_rel_ratio) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio);
        rJacobian[7][9] = (-0.001*mDt) * (((((((var_intracellular_calcium_concentration__Ca_i / pow(((var_intracellular_calcium_concentration__Ca_i + (((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca)) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs)) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca, 2.0)) * var_sarcoplasmic_reticulum_calcium_pump__alpha_up) * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca) + ((((var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca) / (((var_intracellular_calcium_concentration__Ca_i + (((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca)) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs)) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca)) * var_sarcoplasmic_reticulum_calcium_pump__beta_up)) - (((((var_intracellular_calcium_concentration__Ca_up * pow(var_sarcoplasmic_reticulum_calcium_pump__K_cyca, 2.0)) * pow(var_sarcoplasmic_reticulum_calcium_pump__K_xcs, 2.0)) / pow(var_sarcoplasmic_reticulum_calcium_pump__K_srca, 2.0)) / pow(((var_intracellular_calcium_concentration__Ca_i + (((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca)) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs)) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca, 2.0)) * var_sarcoplasmic_reticulum_calcium_pump__beta_up));
        rJacobian[7][10] = 0.0;
        rJacobian[7][11] = (((((0.159154943 * 0.001*mDt) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) / var_membrane__F) * (((((((((((((-2.0) * (1.0 - var_sodium_calcium_exchanger__FRiNaCa)) * var_sodium_calcium_exchanger__k_NaCa) * exp(((((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T)) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_sodium_calcium_exchanger__n_NaCa) / var_intracellular_sodium_concentration__Na_i) * var_extracellular_calcium_concentration__Ca_o) / (1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_i * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))))) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i))) + ((((((((((2.0 * (1.0 - var_sodium_calcium_exchanger__FRiNaCa)) * var_sodium_calcium_exchanger__k_NaCa) * (((exp(((((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_extracellular_calcium_concentration__Ca_o) - ((exp((((((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) * var_intracellular_calcium_concentration__Ca_i))) / pow(1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_i * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))), 2.0)) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i))) * var_sodium_calcium_exchanger__d_NaCa) * var_extracellular_calcium_concentration__Ca_o) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_sodium_calcium_exchanger__n_NaCa) / var_intracellular_sodium_concentration__Na_i)) - (((((((((2.0 * var_sodium_calcium_exchanger__FRiNaCa) * var_sodium_calcium_exchanger__k_NaCa) * exp(((((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T)) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_sodium_calcium_exchanger__n_NaCa) / var_intracellular_sodium_concentration__Na_i) * var_extracellular_calcium_concentration__Ca_o) / (1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_ds * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))))) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_ds)))) + ((((((((((2.0 * var_sodium_calcium_exchanger__FRiNaCa) * var_sodium_calcium_exchanger__k_NaCa) * (((exp(((((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_extracellular_calcium_concentration__Ca_o) - ((exp((((((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) * var_intracellular_calcium_concentration__Ca_ds))) / pow(1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_ds * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))), 2.0)) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_ds))) * var_sodium_calcium_exchanger__d_NaCa) * var_extracellular_calcium_concentration__Ca_o) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_sodium_calcium_exchanger__n_NaCa) / var_intracellular_sodium_concentration__Na_i));
        rJacobian[8][0] = 0.0;
        rJacobian[8][1] = 0.0;
        rJacobian[8][2] = (0.001*mDt * ((((2.0 * var_calcium_release__ActFrac) / pow(var_calcium_release__ActFrac + 0.25, 2.0)) * var_calcium_release__K_m_rel) - (((2.0 * pow(var_calcium_release__ActFrac, 2.0)) / pow(var_calcium_release__ActFrac + 0.25, 3.0)) * var_calcium_release__K_m_rel))) * var_intracellular_calcium_concentration__Ca_rel;
        rJacobian[8][3] = 0.0;
        rJacobian[8][4] = 0.0;
        rJacobian[8][5] = 0.0;
        rJacobian[8][6] = 0.0;
        rJacobian[8][7] = 0.0;
        rJacobian[8][8] = 1.0 - (0.001*mDt * (((((-50.0) * var_intracellular_calcium_concentration__V_up_ratio) / var_intracellular_calcium_concentration__V_rel_ratio) - ((pow(var_calcium_release__ActFrac, 2.0) / pow(var_calcium_release__ActFrac + 0.25, 2.0)) * var_calcium_release__K_m_rel)) - var_calcium_release__K_leak_rate));
        rJacobian[8][9] = (((-50.0) * 0.001*mDt) * var_intracellular_calcium_concentration__V_up_ratio) / var_intracellular_calcium_concentration__V_rel_ratio;
        rJacobian[8][10] = 0.0;
        rJacobian[8][11] = 0.0;
        rJacobian[9][0] = 0.0;
        rJacobian[9][1] = 0.0;
        rJacobian[9][2] = 0.0;
        rJacobian[9][3] = 0.0;
        rJacobian[9][4] = 0.0;
        rJacobian[9][5] = 0.0;
        rJacobian[9][6] = 0.0;
        rJacobian[9][7] = (((-0.001*mDt) * (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) / var_intracellular_calcium_concentration__V_up_ratio) * ((((1.0 / (((var_intracellular_calcium_concentration__Ca_i + (((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca)) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs)) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca)) * var_sarcoplasmic_reticulum_calcium_pump__alpha_up) - ((var_intracellular_calcium_concentration__Ca_i / pow(((var_intracellular_calcium_concentration__Ca_i + (((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca)) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs)) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca, 2.0)) * var_sarcoplasmic_reticulum_calcium_pump__alpha_up)) + (((((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca) / pow(((var_intracellular_calcium_concentration__Ca_i + (((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca)) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs)) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca, 2.0)) * var_sarcoplasmic_reticulum_calcium_pump__beta_up));
        rJacobian[9][8] = (-50.0) * 0.001*mDt;
        rJacobian[9][9] = 1.0 - (0.001*mDt * ((((((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio) / var_intracellular_calcium_concentration__V_up_ratio) * ((((((((-var_intracellular_calcium_concentration__Ca_i) / pow(((var_intracellular_calcium_concentration__Ca_i + (((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca)) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs)) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca, 2.0)) * var_sarcoplasmic_reticulum_calcium_pump__alpha_up) * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca) - ((((var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca) / (((var_intracellular_calcium_concentration__Ca_i + (((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca)) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs)) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca)) * var_sarcoplasmic_reticulum_calcium_pump__beta_up)) + (((((var_intracellular_calcium_concentration__Ca_up * pow(var_sarcoplasmic_reticulum_calcium_pump__K_cyca, 2.0)) * pow(var_sarcoplasmic_reticulum_calcium_pump__K_xcs, 2.0)) / pow(var_sarcoplasmic_reticulum_calcium_pump__K_srca, 2.0)) / pow(((var_intracellular_calcium_concentration__Ca_i + (((var_intracellular_calcium_concentration__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_cyca) * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca)) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs)) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca, 2.0)) * var_sarcoplasmic_reticulum_calcium_pump__beta_up))) - 50.0));
        rJacobian[9][10] = 0.0;
        rJacobian[9][11] = 0.0;
        rJacobian[10][0] = (((((((((((((0.3183098861 * 0.001*mDt) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) * (1.0 - var_L_type_Ca_channel__FrICa)) * var_L_type_Ca_channel__P_CaK) * var_L_type_Ca_channel__P_Ca_L) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * (var_membrane__V - 50.0)) / var_membrane__R) / var_membrane__T) / (1.0 - exp(((((-var_membrane__V) + 50.0) * var_membrane__F) / var_membrane__R) / var_membrane__T))) * ((var_intracellular_potassium_concentration__K_i * exp(((50.0 * var_membrane__F) / var_membrane__R) / var_membrane__T)) - (var_extracellular_potassium_concentration__K_o * exp(((((-var_membrane__V) + 50.0) * var_membrane__F) / var_membrane__R) / var_membrane__T)));
        rJacobian[10][1] = (((((((((((((0.3183098861 * 0.001*mDt) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) * var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaK) * var_L_type_Ca_channel__P_Ca_L) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * (var_membrane__V - 50.0)) / var_membrane__R) / var_membrane__T) / (1.0 - exp(((((-var_membrane__V) + 50.0) * var_membrane__F) / var_membrane__R) / var_membrane__T))) * ((var_intracellular_potassium_concentration__K_i * exp(((50.0 * var_membrane__F) / var_membrane__R) / var_membrane__T)) - (var_extracellular_potassium_concentration__K_o * exp(((((-var_membrane__V) + 50.0) * var_membrane__F) / var_membrane__R) / var_membrane__T)));
        rJacobian[10][2] = 0.0;
        rJacobian[10][3] = 0.0;
        rJacobian[10][4] = 0.0;
        rJacobian[10][5] = 0.0;
        rJacobian[10][6] = 0.0;
        rJacobian[10][7] = 0.0;
        rJacobian[10][8] = 0.0;
        rJacobian[10][9] = 0.0;
        rJacobian[10][10] = 1.0 + ((((((0.3183098861 * 0.001*mDt) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) / var_membrane__F) * (((((((((((((var_time_independent_potassium_current__g_K1 * var_extracellular_potassium_concentration__K_o) / (var_extracellular_potassium_concentration__K_o + var_time_independent_potassium_current__K_mk1)) * var_membrane__R) * var_membrane__T) / var_membrane__F) / var_intracellular_potassium_concentration__K_i) / (1.0 + exp((((1.25 * ((var_membrane__V - (((var_membrane__R * var_membrane__T) / var_membrane__F) * log(var_extracellular_potassium_concentration__K_o / var_intracellular_potassium_concentration__K_i))) - 10.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T))) - (((((((1.25 * var_time_independent_potassium_current__g_K1) * var_extracellular_potassium_concentration__K_o) / (var_extracellular_potassium_concentration__K_o + var_time_independent_potassium_current__K_mk1)) * (var_membrane__V - (((var_membrane__R * var_membrane__T) / var_membrane__F) * log(var_extracellular_potassium_concentration__K_o / var_intracellular_potassium_concentration__K_i)))) / pow(1.0 + exp((((1.25 * ((var_membrane__V - (((var_membrane__R * var_membrane__T) / var_membrane__F) * log(var_extracellular_potassium_concentration__K_o / var_intracellular_potassium_concentration__K_i))) - 10.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T), 2.0)) / var_intracellular_potassium_concentration__K_i) * exp((((1.25 * ((var_membrane__V - (((var_membrane__R * var_membrane__T) / var_membrane__F) * log(var_extracellular_potassium_concentration__K_o / var_intracellular_potassium_concentration__K_i))) - 10.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T))) + (((((((var_rapid_delayed_rectifier_potassium_current__g_Kr1 * var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1) + (var_rapid_delayed_rectifier_potassium_current__g_Kr2 * var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2)) / (1.0 + exp((0.04464285714 * var_membrane__V) + 0.4017857143))) * var_membrane__R) * var_membrane__T) / var_membrane__F) / var_intracellular_potassium_concentration__K_i)) + (((((var_slow_delayed_rectifier_potassium_current__g_Ks * pow(var_slow_delayed_rectifier_potassium_current_xs_gate__xs, 2.0)) * var_membrane__R) * var_membrane__T) / var_membrane__F) / (var_intracellular_potassium_concentration__K_i + (var_reversal_potentials__P_kna * var_intracellular_sodium_concentration__Na_i)))) + ((((((((((((1.0 - var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaK) * var_L_type_Ca_channel__P_Ca_L) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * var_L_type_Ca_channel_f2_gate__f2) * (var_membrane__V - 50.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T) / (1.0 - exp(((((-var_membrane__V) + 50.0) * var_membrane__F) / var_membrane__R) / var_membrane__T))) * exp(((50.0 * var_membrane__F) / var_membrane__R) / var_membrane__T))) + (((((((((((var_L_type_Ca_channel__FrICa * var_L_type_Ca_channel__P_CaK) * var_L_type_Ca_channel__P_Ca_L) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * var_L_type_Ca_channel_f2ds_gate__f2ds) * (var_membrane__V - 50.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T) / (1.0 - exp(((((-var_membrane__V) + 50.0) * var_membrane__F) / var_membrane__R) / var_membrane__T))) * exp(((50.0 * var_membrane__F) / var_membrane__R) / var_membrane__T))) + ((((((var_transient_outward_current__g_to * (var_transient_outward_current__g_tos + (var_transient_outward_current_s_gate__s * (1.0 - var_transient_outward_current__g_tos)))) * var_transient_outward_current_r_gate__r) * var_membrane__R) * var_membrane__T) / var_membrane__F) / var_intracellular_potassium_concentration__K_i)));
        rJacobian[10][11] = (((((0.3183098861 * 0.001*mDt) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) / var_membrane__F) * ((((((((var_slow_delayed_rectifier_potassium_current__g_Ks * pow(var_slow_delayed_rectifier_potassium_current_xs_gate__xs, 2.0)) * var_membrane__R) * var_membrane__T) / var_membrane__F) / (var_intracellular_potassium_concentration__K_i + (var_reversal_potentials__P_kna * var_intracellular_sodium_concentration__Na_i))) * var_reversal_potentials__P_kna) - ((((2.0 * var_sodium_potassium_pump__i_NaK_max) * var_extracellular_potassium_concentration__K_o) / (var_sodium_potassium_pump__K_mK + var_extracellular_potassium_concentration__K_o)) / (var_sodium_potassium_pump__K_mNa + var_intracellular_sodium_concentration__Na_i))) + (((((2.0 * var_sodium_potassium_pump__i_NaK_max) * var_extracellular_potassium_concentration__K_o) / (var_sodium_potassium_pump__K_mK + var_extracellular_potassium_concentration__K_o)) * var_intracellular_sodium_concentration__Na_i) / pow(var_sodium_potassium_pump__K_mNa + var_intracellular_sodium_concentration__Na_i, 2.0)));
        rJacobian[11][0] = (((((((((((((0.3183098861 * 0.001*mDt) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) * (1.0 - var_L_type_Ca_channel__FrICa)) * var_L_type_Ca_channel__P_CaNa) * var_L_type_Ca_channel__P_Ca_L) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * (var_membrane__V - 50.0)) / var_membrane__R) / var_membrane__T) / (1.0 - exp(((((-var_membrane__V) + 50.0) * var_membrane__F) / var_membrane__R) / var_membrane__T))) * ((var_intracellular_sodium_concentration__Na_i * exp(((50.0 * var_membrane__F) / var_membrane__R) / var_membrane__T)) - (var_extracellular_sodium_concentration__Na_o * exp(((((-var_membrane__V) + 50.0) * var_membrane__F) / var_membrane__R) / var_membrane__T)));
        rJacobian[11][1] = (((((((((((((0.3183098861 * 0.001*mDt) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) * var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaNa) * var_L_type_Ca_channel__P_Ca_L) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * (var_membrane__V - 50.0)) / var_membrane__R) / var_membrane__T) / (1.0 - exp(((((-var_membrane__V) + 50.0) * var_membrane__F) / var_membrane__R) / var_membrane__T))) * ((var_intracellular_sodium_concentration__Na_i * exp(((50.0 * var_membrane__F) / var_membrane__R) / var_membrane__T)) - (var_extracellular_sodium_concentration__Na_o * exp(((((-var_membrane__V) + 50.0) * var_membrane__F) / var_membrane__R) / var_membrane__T)));
        rJacobian[11][2] = 0.0;
        rJacobian[11][3] = 0.0;
        rJacobian[11][4] = 0.0;
        rJacobian[11][5] = 0.0;
        rJacobian[11][6] = 0.0;
        rJacobian[11][7] = (((((0.3183098861 * 0.001*mDt) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) / var_membrane__F) * (((((((((-3.0) * (1.0 - var_sodium_calcium_exchanger__FRiNaCa)) * var_sodium_calcium_exchanger__k_NaCa) * exp((((((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T)) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) / (1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_i * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))))) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i))) - (((((((3.0 * (1.0 - var_sodium_calcium_exchanger__FRiNaCa)) * var_sodium_calcium_exchanger__k_NaCa) * (((exp(((((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_extracellular_calcium_concentration__Ca_o) - ((exp((((((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) * var_intracellular_calcium_concentration__Ca_i))) / pow(1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_i * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))), 2.0)) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i))) * var_sodium_calcium_exchanger__d_NaCa) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa))) - (((((434.7826086 * (1.0 - var_sodium_calcium_exchanger__FRiNaCa)) * var_sodium_calcium_exchanger__k_NaCa) * (((exp(((((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_extracellular_calcium_concentration__Ca_o) - ((exp((((((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) * var_intracellular_calcium_concentration__Ca_i))) / (1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_i * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))))) / pow(1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i), 2.0)));
        rJacobian[11][8] = 0.0;
        rJacobian[11][9] = 0.0;
        rJacobian[11][10] = ((((((((((0.03819718633 * 0.001*mDt) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) / pow(var_membrane__F, 2.0)) * var_fast_sodium_current__g_Na) * pow(var_fast_sodium_current_m_gate__m, 3.0)) * var_fast_sodium_current_h_gate__h) * var_membrane__R) * var_membrane__T) / (var_intracellular_sodium_concentration__Na_i + (0.12 * var_intracellular_potassium_concentration__K_i));
        rJacobian[11][11] = 1.0 + ((((((0.3183098861 * 0.001*mDt) / pow(var_intracellular_calcium_concentration__radius, 2.0)) / var_intracellular_calcium_concentration__length) / (((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio)) / var_membrane__F) * ((((((((((((((var_fast_sodium_current__g_Na * pow(var_fast_sodium_current_m_gate__m, 3.0)) * var_fast_sodium_current_h_gate__h) * var_membrane__R) * var_membrane__T) / var_membrane__F) / (var_intracellular_sodium_concentration__Na_i + (0.12 * var_intracellular_potassium_concentration__K_i))) + (((((var_persistent_sodium_current__g_pna / (1.0 + exp((((-1.0) / 8.0) * var_membrane__V) - (13.0 / 2.0)))) * var_membrane__R) * var_membrane__T) / var_membrane__F) / var_intracellular_sodium_concentration__Na_i)) + ((((var_sodium_background_current__g_bna * var_membrane__R) * var_membrane__T) / var_membrane__F) / var_intracellular_sodium_concentration__Na_i)) + ((((3.0 * var_sodium_potassium_pump__i_NaK_max) * var_extracellular_potassium_concentration__K_o) / (var_sodium_potassium_pump__K_mK + var_extracellular_potassium_concentration__K_o)) / (var_sodium_potassium_pump__K_mNa + var_intracellular_sodium_concentration__Na_i))) - (((((3.0 * var_sodium_potassium_pump__i_NaK_max) * var_extracellular_potassium_concentration__K_o) / (var_sodium_potassium_pump__K_mK + var_extracellular_potassium_concentration__K_o)) * var_intracellular_sodium_concentration__Na_i) / pow(var_sodium_potassium_pump__K_mNa + var_intracellular_sodium_concentration__Na_i, 2.0))) + (((((((((3.0 * (1.0 - var_sodium_calcium_exchanger__FRiNaCa)) * var_sodium_calcium_exchanger__k_NaCa) * exp(((((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T)) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_sodium_calcium_exchanger__n_NaCa) / var_intracellular_sodium_concentration__Na_i) * var_extracellular_calcium_concentration__Ca_o) / (1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_i * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))))) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i)))) - ((((((((((3.0 * (1.0 - var_sodium_calcium_exchanger__FRiNaCa)) * var_sodium_calcium_exchanger__k_NaCa) * (((exp(((((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_extracellular_calcium_concentration__Ca_o) - ((exp((((((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0)) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) * var_intracellular_calcium_concentration__Ca_i))) / pow(1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_intracellular_calcium_concentration__Ca_i * pow(var_extracellular_sodium_concentration__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_extracellular_calcium_concentration__Ca_o * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)))), 2.0)) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i))) * var_sodium_calcium_exchanger__d_NaCa) * var_extracellular_calcium_concentration__Ca_o) * pow(var_intracellular_sodium_concentration__Na_i, var_sodium_calcium_exchanger__n_NaCa)) * var_sodium_calcium_exchanger__n_NaCa) / var_intracellular_sodium_concentration__Na_i)) + ((((((((((((1.0 - var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaNa) * var_L_type_Ca_channel__P_Ca_L) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * var_L_type_Ca_channel_f2_gate__f2) * (var_membrane__V - 50.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T) / (1.0 - exp(((((-var_membrane__V) + 50.0) * var_membrane__F) / var_membrane__R) / var_membrane__T))) * exp(((50.0 * var_membrane__F) / var_membrane__R) / var_membrane__T))) + (((((((((((var_L_type_Ca_channel__FrICa * var_L_type_Ca_channel__P_CaNa) * var_L_type_Ca_channel__P_Ca_L) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * var_L_type_Ca_channel_f2ds_gate__f2ds) * (var_membrane__V - 50.0)) * var_membrane__F) / var_membrane__R) / var_membrane__T) / (1.0 - exp(((((-var_membrane__V) + 50.0) * var_membrane__F) / var_membrane__R) / var_membrane__T))) * exp(((50.0 * var_membrane__F) / var_membrane__R) / var_membrane__T))));
    }

protected:
    void UpdateTransmembranePotential(double var_environment__time)
    {
        // Time units: second
        var_environment__time *= 0.001;
        std::vector<double>& rY = rGetStateVariables();
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -92.849333
        double var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1 = rY[1];
        // Units: dimensionless; Initial value: 1.03e-5
        double var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2 = rY[2];
        // Units: dimensionless; Initial value: 2e-7
        double var_slow_delayed_rectifier_potassium_current_xs_gate__xs = rY[3];
        // Units: dimensionless; Initial value: 0.001302
        double var_fast_sodium_current_m_gate__m = rY[4];
        // Units: dimensionless; Initial value: 0.0016203
        double var_fast_sodium_current_h_gate__h = rY[5];
        // Units: dimensionless; Initial value: 0.9944036
        double var_L_type_Ca_channel_d_gate__d = rY[6];
        // Units: dimensionless; Initial value: 0
        double var_L_type_Ca_channel_f_gate__f = rY[7];
        // Units: dimensionless; Initial value: 1
        double var_L_type_Ca_channel_f2_gate__f2 = rY[8];
        // Units: dimensionless; Initial value: 0.9349197
        double var_L_type_Ca_channel_f2ds_gate__f2ds = rY[9];
        // Units: dimensionless; Initial value: 0.9651958
        double var_transient_outward_current_s_gate__s = rY[10];
        // Units: dimensionless; Initial value: 0.9948645
        double var_transient_outward_current_r_gate__r = rY[11];
        // Units: dimensionless; Initial value: 0
        double var_intracellular_sodium_concentration__Na_i = rY[14];
        // Units: millimolar; Initial value: 7.3321223
        double var_intracellular_potassium_concentration__K_i = rY[15];
        // Units: millimolar; Initial value: 136.5644281
        double var_intracellular_calcium_concentration__Ca_i = rY[16];
        // Units: millimolar; Initial value: 1.4e-5
        double var_intracellular_calcium_concentration__Ca_ds = rY[17];
        // Units: millimolar; Initial value: 1.88e-5
        
        const double var_membrane__R = 8314.472;
        const double var_membrane__T = 310.0;
        const double var_membrane__F = 96485.3415;
        const double var_membrane__Cm = 9.5e-05;
        //double var_membrane__time = var_environment__time;
        double var_reversal_potentials__K_i = var_intracellular_potassium_concentration__K_i;
        double var_reversal_potentials__R = var_membrane__R;
        double var_reversal_potentials__T = var_membrane__T;
        double var_reversal_potentials__F = var_membrane__F;
        const double var_extracellular_potassium_concentration__K_o = 4.0;
        double var_reversal_potentials__K_o = var_extracellular_potassium_concentration__K_o;
        double var_reversal_potentials__E_K = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__K_o / var_reversal_potentials__K_i);
        double var_time_independent_potassium_current__E_K = var_reversal_potentials__E_K;
        double var_time_independent_potassium_current__K_o = var_extracellular_potassium_concentration__K_o;
        double var_time_independent_potassium_current__R = var_membrane__R;
        double var_time_independent_potassium_current__V = var_membrane__V;
        double var_time_independent_potassium_current__T = var_membrane__T;
        const double var_time_independent_potassium_current__K_mk1 = 10.0;
        const double var_time_independent_potassium_current__g_K1 = 0.5;
        double var_time_independent_potassium_current__F = var_membrane__F;
        double var_time_independent_potassium_current__i_K1 = (((var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K_o) / (var_time_independent_potassium_current__K_o + var_time_independent_potassium_current__K_mk1)) * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K)) / (1.0 + exp((((var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K) - 10.0) * var_time_independent_potassium_current__F * 1.25) / (var_time_independent_potassium_current__R * var_time_independent_potassium_current__T)));
        double var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
        double var_transient_outward_current__s = var_transient_outward_current_s_gate__s;
        double var_transient_outward_current__r = var_transient_outward_current_r_gate__r;
        const double var_transient_outward_current__g_to = 0.005;
        double var_transient_outward_current__V = var_membrane__V;
        double var_transient_outward_current__E_K = var_reversal_potentials__E_K;
        const double var_transient_outward_current__g_tos = 0.0;
        double var_transient_outward_current__i_to = var_transient_outward_current__g_to * (var_transient_outward_current__g_tos + (var_transient_outward_current__s * (1.0 - var_transient_outward_current__g_tos))) * var_transient_outward_current__r * (var_transient_outward_current__V - var_transient_outward_current__E_K);
        double var_membrane__i_to = var_transient_outward_current__i_to;
        const double var_rapid_delayed_rectifier_potassium_current__g_Kr2 = 0.0013;
        const double var_rapid_delayed_rectifier_potassium_current__g_Kr1 = 0.0021;
        double var_rapid_delayed_rectifier_potassium_current__xr1 = var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1;
        double var_rapid_delayed_rectifier_potassium_current__xr2 = var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2;
        double var_rapid_delayed_rectifier_potassium_current__V = var_membrane__V;
        double var_rapid_delayed_rectifier_potassium_current__E_K = var_reversal_potentials__E_K;
        double var_rapid_delayed_rectifier_potassium_current__i_Kr = ((((var_rapid_delayed_rectifier_potassium_current__g_Kr1 * var_rapid_delayed_rectifier_potassium_current__xr1) + (var_rapid_delayed_rectifier_potassium_current__g_Kr2 * var_rapid_delayed_rectifier_potassium_current__xr2)) * 1.0) / (1.0 + exp((var_rapid_delayed_rectifier_potassium_current__V + 9.0) / 22.4))) * (var_rapid_delayed_rectifier_potassium_current__V - var_rapid_delayed_rectifier_potassium_current__E_K);
        double var_membrane__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
        double var_slow_delayed_rectifier_potassium_current__xs = var_slow_delayed_rectifier_potassium_current_xs_gate__xs;
        const double var_extracellular_sodium_concentration__Na_o = 140.0;
        double var_reversal_potentials__Na_o = var_extracellular_sodium_concentration__Na_o;
        double var_reversal_potentials__Na_i = var_intracellular_sodium_concentration__Na_i;
        const double var_reversal_potentials__P_kna = 0.03;
        double var_reversal_potentials__E_Ks = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log((var_reversal_potentials__K_o + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_o)) / (var_reversal_potentials__K_i + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_i)));
        double var_slow_delayed_rectifier_potassium_current__E_Ks = var_reversal_potentials__E_Ks;
        const double var_slow_delayed_rectifier_potassium_current__g_Ks = 0.0026;
        double var_slow_delayed_rectifier_potassium_current__V = var_membrane__V;
        double var_slow_delayed_rectifier_potassium_current__i_Ks = var_slow_delayed_rectifier_potassium_current__g_Ks * pow(var_slow_delayed_rectifier_potassium_current__xs, 2.0) * (var_slow_delayed_rectifier_potassium_current__V - var_slow_delayed_rectifier_potassium_current__E_Ks);
        double var_membrane__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
        double var_L_type_Ca_channel__d = var_L_type_Ca_channel_d_gate__d;
        const double var_L_type_Ca_channel__FrICa = 1.0;
        double var_L_type_Ca_channel__f = var_L_type_Ca_channel_f_gate__f;
        double var_L_type_Ca_channel__K_o = var_extracellular_potassium_concentration__K_o;
        double var_L_type_Ca_channel__K_i = var_intracellular_potassium_concentration__K_i;
        double var_L_type_Ca_channel__F = var_membrane__F;
        const double var_L_type_Ca_channel__P_Ca_L = 0.1;
        double var_L_type_Ca_channel__T = var_membrane__T;
        const double var_L_type_Ca_channel__P_CaK = 0.002;
        double var_L_type_Ca_channel__V = var_membrane__V;
        double var_L_type_Ca_channel__f2 = var_L_type_Ca_channel_f2_gate__f2;
        double var_L_type_Ca_channel__R = var_membrane__R;
        double var_L_type_Ca_channel__i_Ca_L_K_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaK * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__K_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__K_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_membrane__i_Ca_L_K_cyt = var_L_type_Ca_channel__i_Ca_L_K_cyt;
        double var_L_type_Ca_channel__f2ds = var_L_type_Ca_channel_f2ds_gate__f2ds;
        double var_L_type_Ca_channel__i_Ca_L_K_ds = (((var_L_type_Ca_channel__FrICa * var_L_type_Ca_channel__P_CaK * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__K_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__K_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_membrane__i_Ca_L_K_ds = var_L_type_Ca_channel__i_Ca_L_K_ds;
        const double var_sodium_potassium_pump__i_NaK_max = 0.7;
        double var_sodium_potassium_pump__Na_i = var_intracellular_sodium_concentration__Na_i;
        double var_sodium_potassium_pump__K_o = var_extracellular_potassium_concentration__K_o;
        const double var_sodium_potassium_pump__K_mNa = 40.0;
        const double var_sodium_potassium_pump__K_mK = 1.0;
        double var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__i_NaK_max * var_sodium_potassium_pump__K_o) / (var_sodium_potassium_pump__K_mK + var_sodium_potassium_pump__K_o)) * var_sodium_potassium_pump__Na_i) / (var_sodium_potassium_pump__K_mNa + var_sodium_potassium_pump__Na_i);
        double var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
        const double var_fast_sodium_current__g_Na = 2.5;
        double var_fast_sodium_current__h = var_fast_sodium_current_h_gate__h;
        double var_fast_sodium_current__V = var_membrane__V;
        double var_reversal_potentials__E_mh = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log((var_reversal_potentials__Na_o + (0.12 * var_reversal_potentials__K_o)) / (var_reversal_potentials__Na_i + (0.12 * var_reversal_potentials__K_i)));
        double var_fast_sodium_current__E_mh = var_reversal_potentials__E_mh;
        double var_fast_sodium_current__m = var_fast_sodium_current_m_gate__m;
        double var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * (var_fast_sodium_current__V - var_fast_sodium_current__E_mh);
        double var_membrane__i_Na = var_fast_sodium_current__i_Na;
        double var_sodium_background_current__V = var_membrane__V;
        double var_reversal_potentials__E_Na = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Na_o / var_reversal_potentials__Na_i);
        double var_sodium_background_current__E_Na = var_reversal_potentials__E_Na;
        const double var_sodium_background_current__g_bna = 0.0006;
        double var_sodium_background_current__i_b_Na = var_sodium_background_current__g_bna * (var_sodium_background_current__V - var_sodium_background_current__E_Na);
        double var_membrane__i_b_Na = var_sodium_background_current__i_b_Na;
        const double var_persistent_sodium_current__g_pna = 0.004;
        double var_persistent_sodium_current__V = var_membrane__V;
        double var_persistent_sodium_current__E_Na = var_reversal_potentials__E_Na;
        double var_persistent_sodium_current__i_p_Na = ((var_persistent_sodium_current__g_pna * 1.0) / (1.0 + exp((-(var_persistent_sodium_current__V + 52.0)) / 8.0))) * (var_persistent_sodium_current__V - var_persistent_sodium_current__E_Na);
        double var_membrane__i_p_Na = var_persistent_sodium_current__i_p_Na;
        const double var_L_type_Ca_channel__P_CaNa = 0.01;
        double var_L_type_Ca_channel__Na_o = var_extracellular_sodium_concentration__Na_o;
        double var_L_type_Ca_channel__Na_i = var_intracellular_sodium_concentration__Na_i;
        double var_L_type_Ca_channel__i_Ca_L_Na_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaNa * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Na_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Na_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_membrane__i_Ca_L_Na_cyt = var_L_type_Ca_channel__i_Ca_L_Na_cyt;
        double var_L_type_Ca_channel__i_Ca_L_Na_ds = (((var_L_type_Ca_channel__FrICa * var_L_type_Ca_channel__P_CaNa * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Na_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Na_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_membrane__i_Ca_L_Na_ds = var_L_type_Ca_channel__i_Ca_L_Na_ds;
        double var_sodium_calcium_exchanger__Na_i = var_intracellular_sodium_concentration__Na_i;
        const double var_sodium_calcium_exchanger__n_NaCa = 3.0;
        const double var_sodium_calcium_exchanger__gamma = 0.5;
        double var_sodium_calcium_exchanger__F = var_membrane__F;
        double var_sodium_calcium_exchanger__Na_o = var_extracellular_sodium_concentration__Na_o;
        const double var_sodium_calcium_exchanger__FRiNaCa = 0.001;
        double var_sodium_calcium_exchanger__R = var_membrane__R;
        double var_sodium_calcium_exchanger__Ca_i = var_intracellular_calcium_concentration__Ca_i;
        double var_sodium_calcium_exchanger__T = var_membrane__T;
        double var_sodium_calcium_exchanger__V = var_membrane__V;
        const double var_sodium_calcium_exchanger__d_NaCa = 0.0;
        const double var_extracellular_calcium_concentration__Ca_o = 2.0;
        double var_sodium_calcium_exchanger__Ca_o = var_extracellular_calcium_concentration__Ca_o;
        const double var_sodium_calcium_exchanger__k_NaCa = 0.0005;
        double var_sodium_calcium_exchanger__i_NaCa_cyt = ((1.0 - var_sodium_calcium_exchanger__FRiNaCa) * var_sodium_calcium_exchanger__k_NaCa * ((exp((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_o) - (exp(((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_i))) / ((1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_sodium_calcium_exchanger__Ca_i * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_sodium_calcium_exchanger__Ca_o * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa))))) * (1.0 + (var_sodium_calcium_exchanger__Ca_i / 0.0069)));
        double var_membrane__i_NaCa_cyt = var_sodium_calcium_exchanger__i_NaCa_cyt;
        double var_sodium_calcium_exchanger__Ca_ds = var_intracellular_calcium_concentration__Ca_ds;
        double var_sodium_calcium_exchanger__i_NaCa_ds = (var_sodium_calcium_exchanger__FRiNaCa * var_sodium_calcium_exchanger__k_NaCa * ((exp((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_o) - (exp(((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_ds))) / ((1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_sodium_calcium_exchanger__Ca_ds * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_sodium_calcium_exchanger__Ca_o * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa))))) * (1.0 + (var_sodium_calcium_exchanger__Ca_ds / 0.0069)));
        double var_membrane__i_NaCa_ds = var_sodium_calcium_exchanger__i_NaCa_ds;
        double var_L_type_Ca_channel__Ca_i = var_intracellular_calcium_concentration__Ca_i;
        double var_L_type_Ca_channel__Ca_o = var_extracellular_calcium_concentration__Ca_o;
        double var_L_type_Ca_channel__i_Ca_L_Ca_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * 4.0 * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Ca_i * exp((100.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Ca_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_membrane__i_Ca_L_Ca_cyt = var_L_type_Ca_channel__i_Ca_L_Ca_cyt;
        double var_L_type_Ca_channel__i_Ca_L_Ca_ds = (((var_L_type_Ca_channel__FrICa * 4.0 * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Ca_i * exp((100.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Ca_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
        double var_membrane__i_Ca_L_Ca_ds = var_L_type_Ca_channel__i_Ca_L_Ca_ds;
        double var_reversal_potentials__Ca_o = var_extracellular_calcium_concentration__Ca_o;
        double var_reversal_potentials__Ca_i = var_intracellular_calcium_concentration__Ca_i;
        double var_reversal_potentials__E_Ca = ((0.5 * var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Ca_o / var_reversal_potentials__Ca_i);
        double var_calcium_background_current__E_Ca = var_reversal_potentials__E_Ca;
        const double var_calcium_background_current__g_bca = 0.00025;
        double var_calcium_background_current__V = var_membrane__V;
        double var_calcium_background_current__i_b_Ca = var_calcium_background_current__g_bca * (var_calcium_background_current__V - var_calcium_background_current__E_Ca);
        double var_membrane__i_b_Ca = var_calcium_background_current__i_b_Ca;
        //const double var_membrane__stim_end = 100000.0;
        //const double var_membrane__stim_amplitude =  -3.0;
        //const double var_membrane__stim_duration = 0.003;
        //const double var_membrane__stim_period = 1.0;
        //const double var_membrane__stim_start = 0.1;
        double var_membrane__i_Stim = GetStimulus((1.0/0.001)*var_environment__time);
        double d_dt_membrane__V = ((-1.0) / var_membrane__Cm) * (var_membrane__i_Stim + var_membrane__i_K1 + var_membrane__i_to + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_NaK + var_membrane__i_Na + var_membrane__i_b_Na + var_membrane__i_p_Na + var_membrane__i_Ca_L_Na_cyt + var_membrane__i_Ca_L_Na_ds + var_membrane__i_NaCa_cyt + var_membrane__i_NaCa_ds + var_membrane__i_Ca_L_Ca_cyt + var_membrane__i_Ca_L_Ca_ds + var_membrane__i_Ca_L_K_cyt + var_membrane__i_Ca_L_K_ds + var_membrane__i_b_Ca);
        
        rY[0] += mDt * 0.001*d_dt_membrane__V;
    }

    void ComputeOneStepExceptVoltage(double var_environment__time)
    {
        // Time units: second
        var_environment__time *= 0.001;
        std::vector<double>& rY = rGetStateVariables();
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -92.849333
        
        double var_transient_outward_current__V = var_membrane__V;
        double var_rapid_delayed_rectifier_potassium_current__V = var_membrane__V;
        double var_slow_delayed_rectifier_potassium_current__V = var_membrane__V;
        double var_L_type_Ca_channel__V = var_membrane__V;
        double var_fast_sodium_current__V = var_membrane__V;
        double var_rapid_delayed_rectifier_potassium_current_xr1_gate__V = var_rapid_delayed_rectifier_potassium_current__V;
        double var_rapid_delayed_rectifier_potassium_current_xr1_gate__alpha_xr1 = 50.0 / (1.0 + exp((-(var_rapid_delayed_rectifier_potassium_current_xr1_gate__V - 5.0)) / 9.0));
        double var_rapid_delayed_rectifier_potassium_current_xr1_gate__beta_xr1 = 0.05 * exp((-(var_rapid_delayed_rectifier_potassium_current_xr1_gate__V - 20.0)) / 15.0);
        double var_rapid_delayed_rectifier_potassium_current_xr2_gate__V = var_rapid_delayed_rectifier_potassium_current__V;
        double var_rapid_delayed_rectifier_potassium_current_xr2_gate__alpha_xr2 = 50.0 / (1.0 + exp((-(var_rapid_delayed_rectifier_potassium_current_xr2_gate__V - 5.0)) / 9.0));
        double var_rapid_delayed_rectifier_potassium_current_xr2_gate__beta_xr2 = 0.4 * exp(-pow((var_rapid_delayed_rectifier_potassium_current_xr2_gate__V + 30.0) / 30.0, 3.0));
        double var_slow_delayed_rectifier_potassium_current_xs_gate__V = var_slow_delayed_rectifier_potassium_current__V;
        double var_slow_delayed_rectifier_potassium_current_xs_gate__alpha_xs = 14.0 / (1.0 + exp((-(var_slow_delayed_rectifier_potassium_current_xs_gate__V - 40.0)) / 9.0));
        double var_slow_delayed_rectifier_potassium_current_xs_gate__beta_xs = 1.0 * exp((-var_slow_delayed_rectifier_potassium_current_xs_gate__V) / 45.0);
        double var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V;
        double var_fast_sodium_current_m_gate__E0_m = var_fast_sodium_current_m_gate__V + 41.0;
        const double var_fast_sodium_current_m_gate__delta_m = 1e-05;
        double var_fast_sodium_current_m_gate__alpha_m = (fabs(var_fast_sodium_current_m_gate__E0_m) < var_fast_sodium_current_m_gate__delta_m) ? 2000.0 : ((200.0 * var_fast_sodium_current_m_gate__E0_m) / (1.0 - exp((-0.1) * var_fast_sodium_current_m_gate__E0_m)));
        double var_fast_sodium_current_m_gate__beta_m = 8000.0 * exp((-0.056) * (var_fast_sodium_current_m_gate__V + 66.0));
        double var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V;
        const double var_fast_sodium_current_h_gate__shift_h = 0.0;
        double var_fast_sodium_current_h_gate__alpha_h = 20.0 * exp((-0.125) * ((var_fast_sodium_current_h_gate__V + 75.0) - var_fast_sodium_current_h_gate__shift_h));
        double var_fast_sodium_current_h_gate__beta_h = 2000.0 / (1.0 + (320.0 * exp((-0.1) * ((var_fast_sodium_current_h_gate__V + 75.0) - var_fast_sodium_current_h_gate__shift_h))));
        double var_L_type_Ca_channel_d_gate__V = var_L_type_Ca_channel__V;
        double var_L_type_Ca_channel_d_gate__E0_d = (var_L_type_Ca_channel_d_gate__V + 24.0) - 5.0;
        double var_L_type_Ca_channel_d_gate__alpha_d = (fabs(var_L_type_Ca_channel_d_gate__E0_d) < 0.0001) ? 120.0 : ((30.0 * var_L_type_Ca_channel_d_gate__E0_d) / (1.0 - exp((-var_L_type_Ca_channel_d_gate__E0_d) / 4.0)));
        double var_L_type_Ca_channel_d_gate__beta_d = (fabs(var_L_type_Ca_channel_d_gate__E0_d) < 0.0001) ? 120.0 : ((12.0 * var_L_type_Ca_channel_d_gate__E0_d) / (exp(var_L_type_Ca_channel_d_gate__E0_d / 10.0) - 1.0));
        const double var_L_type_Ca_channel_d_gate__speed_d = 3.0;
        double var_L_type_Ca_channel_f_gate__V = var_L_type_Ca_channel__V;
        double var_L_type_Ca_channel_f_gate__E0_f = var_L_type_Ca_channel_f_gate__V + 34.0;
        const double var_L_type_Ca_channel_f_gate__delta_f = 0.0001;
        double var_L_type_Ca_channel_f_gate__alpha_f = (fabs(var_L_type_Ca_channel_f_gate__E0_f) < var_L_type_Ca_channel_f_gate__delta_f) ? 25.0 : ((6.25 * var_L_type_Ca_channel_f_gate__E0_f) / (exp(var_L_type_Ca_channel_f_gate__E0_f / 4.0) - 1.0));
        double var_L_type_Ca_channel_f_gate__beta_f = 12.0 / (1.0 + exp(((-1.0) * (var_L_type_Ca_channel_f_gate__V + 34.0)) / 4.0));
        const double var_L_type_Ca_channel_f_gate__speed_f = 0.3;
        double var_transient_outward_current_s_gate__V = var_transient_outward_current__V;
        double var_transient_outward_current_s_gate__alpha_s = 0.033 * exp((-var_transient_outward_current_s_gate__V) / 17.0);
        double var_transient_outward_current_s_gate__beta_s = 33.0 / (1.0 + exp((-0.125) * (var_transient_outward_current_s_gate__V + 10.0)));
        
        const double _g_0 = var_L_type_Ca_channel_d_gate__speed_d * (var_L_type_Ca_channel_d_gate__alpha_d * 1.0);
        const double _h_0 = var_L_type_Ca_channel_d_gate__speed_d * ((var_L_type_Ca_channel_d_gate__alpha_d * (-1.0)) - (var_L_type_Ca_channel_d_gate__beta_d * 1.0));
        const double _g_1 = var_L_type_Ca_channel_f_gate__speed_f * (var_L_type_Ca_channel_f_gate__alpha_f * 1.0);
        const double _h_1 = var_L_type_Ca_channel_f_gate__speed_f * ((var_L_type_Ca_channel_f_gate__alpha_f * (-1.0)) - (var_L_type_Ca_channel_f_gate__beta_f * 1.0));
        const double _g_2 = var_fast_sodium_current_h_gate__alpha_h * 1.0;
        const double _h_2 = (var_fast_sodium_current_h_gate__alpha_h * (-1.0)) - (var_fast_sodium_current_h_gate__beta_h * 1.0);
        const double _g_3 = var_fast_sodium_current_m_gate__alpha_m * 1.0;
        const double _h_3 = (var_fast_sodium_current_m_gate__alpha_m * (-1.0)) - (var_fast_sodium_current_m_gate__beta_m * 1.0);
        const double _g_4 = var_rapid_delayed_rectifier_potassium_current_xr1_gate__alpha_xr1 * 1.0;
        const double _h_4 = (var_rapid_delayed_rectifier_potassium_current_xr1_gate__alpha_xr1 * (-1.0)) - (var_rapid_delayed_rectifier_potassium_current_xr1_gate__beta_xr1 * 1.0);
        const double _g_5 = var_rapid_delayed_rectifier_potassium_current_xr2_gate__alpha_xr2 * 1.0;
        const double _h_5 = (var_rapid_delayed_rectifier_potassium_current_xr2_gate__alpha_xr2 * (-1.0)) - (var_rapid_delayed_rectifier_potassium_current_xr2_gate__beta_xr2 * 1.0);
        const double _g_6 = var_slow_delayed_rectifier_potassium_current_xs_gate__alpha_xs * 1.0;
        const double _h_6 = (var_slow_delayed_rectifier_potassium_current_xs_gate__alpha_xs * (-1.0)) - (var_slow_delayed_rectifier_potassium_current_xs_gate__beta_xs * 1.0);
        const double _g_7 = 333.0 * (1.0 / (1.0 + exp((-(var_membrane__V + 4.0)) / 5.0)));
        const double _h_7 = 333.0 * (-1.0);
        const double _g_8 = var_transient_outward_current_s_gate__alpha_s * 1.0;
        const double _h_8 = (var_transient_outward_current_s_gate__alpha_s * (-1.0)) - (var_transient_outward_current_s_gate__beta_s * 1.0);
        
        const double dt = mDt*0.001;
        rY[1] = (rY[1] + _g_4*dt) / (1 - _h_4*dt);
        rY[2] = (rY[2] + _g_5*dt) / (1 - _h_5*dt);
        rY[3] = (rY[3] + _g_6*dt) / (1 - _h_6*dt);
        rY[4] = (rY[4] + _g_3*dt) / (1 - _h_3*dt);
        rY[5] = (rY[5] + _g_2*dt) / (1 - _h_2*dt);
        rY[6] = (rY[6] + _g_0*dt) / (1 - _h_0*dt);
        rY[7] = (rY[7] + _g_1*dt) / (1 - _h_1*dt);
        rY[10] = (rY[10] + _g_8*dt) / (1 - _h_8*dt);
        rY[11] = (rY[11] + _g_7*dt) / (1 - _h_7*dt);
        
        double _guess[12] = {rY[8],rY[9],rY[12],rY[13],rY[14],rY[15],rY[16],rY[17],rY[18],rY[19],rY[20],rY[21]};
        CardiacNewtonSolver<12> *_solver = CardiacNewtonSolver<12>::Instance();
        _solver->Solve(*this, _guess);
        rY[8] = _guess[0];
        rY[9] = _guess[1];
        rY[12] = _guess[2];
        rY[13] = _guess[3];
        rY[20] = _guess[4];
        rY[21] = _guess[5];
        rY[17] = _guess[6];
        rY[16] = _guess[7];
        rY[19] = _guess[8];
        rY[18] = _guess[9];
        rY[15] = _guess[10];
        rY[14] = _guess[11];
    }

};

#endif
