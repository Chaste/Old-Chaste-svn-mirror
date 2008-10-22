#ifndef _BackwardEulerFoxModel2002Modified_
#define _BackwardEulerFoxModel2002Modified_

// Model: fox_model_2002
// Processed by pycml - CellML Tools in Python
//     (translate: $Revision: 698 $, pycml: $Revision: 693 $)
// on Wed May  9 09:21:04 2007

#include <cmath>
#include <cassert>
#include "AbstractBackwardEulerCardiacCell.hpp"
#include "CardiacNewtonSolver.hpp"
#include "Exception.hpp"
#include "AbstractStimulusFunction.hpp"

/**
 * Generated from CellML, and P_Ca parameter modified.
 */
class BackwardEulerFoxModel2002Modified : public AbstractBackwardEulerCardiacCell<3>
{
public:
    BackwardEulerFoxModel2002Modified( AbstractStimulusFunction *pIntracellularStimulus,
                       				   AbstractStimulusFunction *pExtracellularStimulus=NULL)
        : AbstractBackwardEulerCardiacCell<3>(13, 0, pIntracellularStimulus, pExtracellularStimulus)
    {
        // Time units: millisecond

        mVariableNames.push_back("V");
        mVariableUnits.push_back("millivolt");
        mInitialConditions.push_back(-94.7);

        mVariableNames.push_back("m");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.00024676);

        mVariableNames.push_back("h");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.99869);

        mVariableNames.push_back("j");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.99887);

        mVariableNames.push_back("X_kr");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.229);

        mVariableNames.push_back("X_ks");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.0001);

        mVariableNames.push_back("X_to");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.00003742);

        mVariableNames.push_back("Y_to");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(1);

        mVariableNames.push_back("f");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.983);

        mVariableNames.push_back("d");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.0001);

        mVariableNames.push_back("f_Ca");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.942);

        mVariableNames.push_back("CaI");
        mVariableUnits.push_back("micromolar");
        mInitialConditions.push_back(0.0472);

        mVariableNames.push_back("Ca_SR");
        mVariableUnits.push_back("micromolar");
        mInitialConditions.push_back(320);

        Init();

    }

    ~BackwardEulerFoxModel2002Modified(void)
    {
    }

    double GetIIonic()
    {
        std::vector<double>& rY = rGetStateVariables();
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -94.7
        double var_fast_sodium_current_m_gate__m = rY[1];
        // Units: dimensionless; Initial value: 0.00024676
        double var_fast_sodium_current_h_gate__h = rY[2];
        // Units: dimensionless; Initial value: 0.99869
        double var_fast_sodium_current_j_gate__j = rY[3];
        // Units: dimensionless; Initial value: 0.99887
        double var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr = rY[4];
        // Units: dimensionless; Initial value: 0.229
        double var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks = rY[5];
        // Units: dimensionless; Initial value: 0.0001
        double var_transient_outward_potassium_current_X_to_gate__X_to = rY[6];
        // Units: dimensionless; Initial value: 0.00003742
        double var_transient_outward_potassium_current_Y_to_gate__Y_to = rY[7];
        // Units: dimensionless; Initial value: 1
        double var_L_type_Ca_current_f_gate__f = rY[8];
        // Units: dimensionless; Initial value: 0.983
        double var_L_type_Ca_current_d_gate__d = rY[9];
        // Units: dimensionless; Initial value: 0.0001
        double var_L_type_Ca_current_f_Ca_gate__f_Ca = rY[10];
        // Units: dimensionless; Initial value: 0.942
        double var_calcium_dynamics__Ca_i = rY[11];
        // Units: micromolar; Initial value: 0.0472

        const double var_membrane__R = 8.314;
        const double var_membrane__T = 310.0;
        const double var_membrane__F = 96.5;
        double var_fast_sodium_current__h = var_fast_sodium_current_h_gate__h;
        const double var_fast_sodium_current__g_Na = 12.8;
        double var_fast_sodium_current__j = var_fast_sodium_current_j_gate__j;
        double var_fast_sodium_current__T = var_membrane__T;
        double var_fast_sodium_current__R = var_membrane__R;
        const double var_standard_ionic_concentrations__Na_i = 10.0;
        double var_fast_sodium_current__Na_i = var_standard_ionic_concentrations__Na_i;
        const double var_standard_ionic_concentrations__Na_o = 138.0;
        double var_fast_sodium_current__Na_o = var_standard_ionic_concentrations__Na_o;
        double var_fast_sodium_current__F = var_membrane__F;
        double var_fast_sodium_current__E_Na = ((var_fast_sodium_current__R * var_fast_sodium_current__T) / var_fast_sodium_current__F) * log(var_fast_sodium_current__Na_o / var_fast_sodium_current__Na_i);
        double var_fast_sodium_current__V = var_membrane__V;
        double var_fast_sodium_current__m = var_fast_sodium_current_m_gate__m;
        double var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na);
        double var_membrane__i_Na = var_fast_sodium_current__i_Na;
        const double var_L_type_Ca_current__P_Ca = 1.26e-05; //2.26e-05;
        double var_L_type_Ca_current__V = var_membrane__V;
        const double var_standard_ionic_concentrations__Ca_o = 2000.0;
        double var_L_type_Ca_current__Ca_o = var_standard_ionic_concentrations__Ca_o;
        double var_L_type_Ca_current__T = var_membrane__T;
        const double var_L_type_Ca_current__C_sc = 1.0;
        double var_L_type_Ca_current__Ca_i = var_calcium_dynamics__Ca_i;
        double var_L_type_Ca_current__R = var_membrane__R;
        double var_L_type_Ca_current__F = var_membrane__F;
        double var_L_type_Ca_current__i_Ca_max = ((((var_L_type_Ca_current__P_Ca / var_L_type_Ca_current__C_sc) * 4.0 * var_L_type_Ca_current__V * pow(var_L_type_Ca_current__F, 2.0)) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) * ((var_L_type_Ca_current__Ca_i * exp((2.0 * var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T))) - (0.341 * var_L_type_Ca_current__Ca_o))) / (exp((2.0 * var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) - 1.0);
        double var_L_type_Ca_current__d = var_L_type_Ca_current_d_gate__d;
        double var_L_type_Ca_current__f = var_L_type_Ca_current_f_gate__f;
        double var_L_type_Ca_current__f_Ca = var_L_type_Ca_current_f_Ca_gate__f_Ca;
        double var_L_type_Ca_current__i_Ca = var_L_type_Ca_current__i_Ca_max * var_L_type_Ca_current__f * var_L_type_Ca_current__d * var_L_type_Ca_current__f_Ca;
        double var_membrane__i_Ca = var_L_type_Ca_current__i_Ca;
        const double var_standard_ionic_concentrations__K_o = 4.0;
        double var_L_type_Ca_current__K_o = var_standard_ionic_concentrations__K_o;
        const double var_L_type_Ca_current__i_Ca_half =  -0.265;
        const double var_L_type_Ca_current__P_CaK = 5.79e-07;
        const double var_standard_ionic_concentrations__K_i = 149.4;
        double var_L_type_Ca_current__K_i = var_standard_ionic_concentrations__K_i;
        double var_L_type_Ca_current__i_CaK = ((((((var_L_type_Ca_current__P_CaK / var_L_type_Ca_current__C_sc) * var_L_type_Ca_current__f * var_L_type_Ca_current__d * var_L_type_Ca_current__f_Ca) / (1.0 + (var_L_type_Ca_current__i_Ca_max / var_L_type_Ca_current__i_Ca_half))) * 1000.0 * var_L_type_Ca_current__V * pow(var_L_type_Ca_current__F, 2.0)) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) * ((var_L_type_Ca_current__K_i * exp((var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T))) - var_L_type_Ca_current__K_o)) / (exp((var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) - 1.0);
        double var_membrane__i_CaK = var_L_type_Ca_current__i_CaK;
        double var_rapid_activating_delayed_rectifiyer_K_current__V = var_membrane__V;
        double var_rapid_activating_delayed_rectifiyer_K_current__R_V = 1.0 / (1.0 + (2.5 * exp(0.1 * (var_rapid_activating_delayed_rectifiyer_K_current__V + 28.0))));
        double var_rapid_activating_delayed_rectifiyer_K_current__X_kr = var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr;
        double var_rapid_activating_delayed_rectifiyer_K_current__R = var_membrane__R;
        double var_rapid_activating_delayed_rectifiyer_K_current__K_o = var_standard_ionic_concentrations__K_o;
        double var_rapid_activating_delayed_rectifiyer_K_current__T = var_membrane__T;
        double var_rapid_activating_delayed_rectifiyer_K_current__K_i = var_standard_ionic_concentrations__K_i;
        double var_rapid_activating_delayed_rectifiyer_K_current__F = var_membrane__F;
        double var_rapid_activating_delayed_rectifiyer_K_current__E_K = ((var_rapid_activating_delayed_rectifiyer_K_current__R * var_rapid_activating_delayed_rectifiyer_K_current__T) / var_rapid_activating_delayed_rectifiyer_K_current__F) * log(var_rapid_activating_delayed_rectifiyer_K_current__K_o / var_rapid_activating_delayed_rectifiyer_K_current__K_i);
        const double var_rapid_activating_delayed_rectifiyer_K_current__g_Kr = 0.0136;
        double var_rapid_activating_delayed_rectifiyer_K_current__i_Kr = var_rapid_activating_delayed_rectifiyer_K_current__g_Kr * var_rapid_activating_delayed_rectifiyer_K_current__R_V * var_rapid_activating_delayed_rectifiyer_K_current__X_kr * sqrt(var_rapid_activating_delayed_rectifiyer_K_current__K_o / 4.0) * (var_rapid_activating_delayed_rectifiyer_K_current__V - var_rapid_activating_delayed_rectifiyer_K_current__E_K);
        double var_membrane__i_Kr = var_rapid_activating_delayed_rectifiyer_K_current__i_Kr;
        const double var_slow_activating_delayed_rectifiyer_K_current__g_Ks = 0.0245;
        double var_slow_activating_delayed_rectifiyer_K_current__V = var_membrane__V;
        double var_slow_activating_delayed_rectifiyer_K_current__K_i = var_standard_ionic_concentrations__K_i;
        double var_slow_activating_delayed_rectifiyer_K_current__F = var_membrane__F;
        double var_slow_activating_delayed_rectifiyer_K_current__K_o = var_standard_ionic_concentrations__K_o;
        double var_slow_activating_delayed_rectifiyer_K_current__T = var_membrane__T;
        double var_slow_activating_delayed_rectifiyer_K_current__Na_o = var_standard_ionic_concentrations__Na_o;
        double var_slow_activating_delayed_rectifiyer_K_current__R = var_membrane__R;
        double var_slow_activating_delayed_rectifiyer_K_current__Na_i = var_standard_ionic_concentrations__Na_i;
        double var_slow_activating_delayed_rectifiyer_K_current__E_Ks = ((var_slow_activating_delayed_rectifiyer_K_current__R * var_slow_activating_delayed_rectifiyer_K_current__T) / var_slow_activating_delayed_rectifiyer_K_current__F) * log((var_slow_activating_delayed_rectifiyer_K_current__K_o + (0.01833 * var_slow_activating_delayed_rectifiyer_K_current__Na_o)) / (var_slow_activating_delayed_rectifiyer_K_current__K_i + (0.01833 * var_slow_activating_delayed_rectifiyer_K_current__Na_i)));
        double var_slow_activating_delayed_rectifiyer_K_current__X_ks = var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks;
        double var_slow_activating_delayed_rectifiyer_K_current__i_Ks = var_slow_activating_delayed_rectifiyer_K_current__g_Ks * pow(var_slow_activating_delayed_rectifiyer_K_current__X_ks, 2.0) * (var_slow_activating_delayed_rectifiyer_K_current__V - var_slow_activating_delayed_rectifiyer_K_current__E_Ks);
        double var_membrane__i_Ks = var_slow_activating_delayed_rectifiyer_K_current__i_Ks;
        const double var_transient_outward_potassium_current__g_to = 0.23815;
        double var_transient_outward_potassium_current__V = var_membrane__V;
        double var_transient_outward_potassium_current__X_to = var_transient_outward_potassium_current_X_to_gate__X_to;
        double var_transient_outward_potassium_current__E_K = var_rapid_activating_delayed_rectifiyer_K_current__E_K;
        double var_transient_outward_potassium_current__Y_to = var_transient_outward_potassium_current_Y_to_gate__Y_to;
        double var_transient_outward_potassium_current__i_to = var_transient_outward_potassium_current__g_to * var_transient_outward_potassium_current__X_to * var_transient_outward_potassium_current__Y_to * (var_transient_outward_potassium_current__V - var_transient_outward_potassium_current__E_K);
        double var_membrane__i_to = var_transient_outward_potassium_current__i_to;
        const double var_time_independent_potassium_current__g_K1 = 2.8;
        double var_time_independent_potassium_current__F = var_membrane__F;
        double var_time_independent_potassium_current_K1_gate__F = var_time_independent_potassium_current__F;
        double var_time_independent_potassium_current__V = var_membrane__V;
        double var_time_independent_potassium_current_K1_gate__V = var_time_independent_potassium_current__V;
        double var_time_independent_potassium_current__T = var_membrane__T;
        double var_time_independent_potassium_current_K1_gate__T = var_time_independent_potassium_current__T;
        double var_time_independent_potassium_current__R = var_membrane__R;
        double var_time_independent_potassium_current_K1_gate__R = var_time_independent_potassium_current__R;
        double var_time_independent_potassium_current__E_K = var_rapid_activating_delayed_rectifiyer_K_current__E_K;
        double var_time_independent_potassium_current_K1_gate__E_K = var_time_independent_potassium_current__E_K;
        double var_time_independent_potassium_current_K1_gate__K1_infinity = 1.0 / (2.0 + exp(((1.62 * var_time_independent_potassium_current_K1_gate__F) / (var_time_independent_potassium_current_K1_gate__R * var_time_independent_potassium_current_K1_gate__T)) * (var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K)));
        double var_time_independent_potassium_current__K1_infinity = var_time_independent_potassium_current_K1_gate__K1_infinity;
        double var_time_independent_potassium_current__K_o = var_standard_ionic_concentrations__K_o;
        const double var_time_independent_potassium_current__K_mK1 = 13.0;
        double var_time_independent_potassium_current__i_K1 = ((var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K1_infinity * var_time_independent_potassium_current__K_o) / (var_time_independent_potassium_current__K_o + var_time_independent_potassium_current__K_mK1)) * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
        double var_plateau_potassium_current__V = var_membrane__V;
        double var_plateau_potassium_current__E_K = var_rapid_activating_delayed_rectifiyer_K_current__E_K;
        const double var_plateau_potassium_current__g_Kp = 0.002216;
        double var_plateau_potassium_current_Kp_gate__V = var_plateau_potassium_current__V;
        double var_plateau_potassium_current_Kp_gate__Kp_V = 1.0 / (1.0 + exp((7.488 - var_plateau_potassium_current_Kp_gate__V) / 5.98));
        double var_plateau_potassium_current__Kp_V = var_plateau_potassium_current_Kp_gate__Kp_V;
        double var_plateau_potassium_current__i_Kp = var_plateau_potassium_current__g_Kp * var_plateau_potassium_current__Kp_V * (var_plateau_potassium_current__V - var_plateau_potassium_current__E_K);
        double var_membrane__i_Kp = var_plateau_potassium_current__i_Kp;
        double var_Na_Ca_exchanger__Na_o = var_standard_ionic_concentrations__Na_o;
        const double var_Na_Ca_exchanger__K_NaCa = 1500.0;
        double var_Na_Ca_exchanger__V = var_membrane__V;
        double var_Na_Ca_exchanger__Ca_o = var_standard_ionic_concentrations__Ca_o;
        const double var_Na_Ca_exchanger__eta = 0.35;
        double var_Na_Ca_exchanger__Na_i = var_standard_ionic_concentrations__Na_i;
        double var_Na_Ca_exchanger__T = var_membrane__T;
        double var_Na_Ca_exchanger__R = var_membrane__R;
        const double var_Na_Ca_exchanger__K_sat = 0.2;
        const double var_Na_Ca_exchanger__K_mCa = 1380.0;
        double var_Na_Ca_exchanger__F = var_membrane__F;
        const double var_Na_Ca_exchanger__K_mNa = 87.5;
        double var_Na_Ca_exchanger__Ca_i = var_calcium_dynamics__Ca_i;
        double var_Na_Ca_exchanger__i_NaCa = (var_Na_Ca_exchanger__K_NaCa / ((pow(var_Na_Ca_exchanger__K_mNa, 3.0) + pow(var_Na_Ca_exchanger__Na_o, 3.0)) * (var_Na_Ca_exchanger__K_mCa + var_Na_Ca_exchanger__Ca_o) * (1.0 + (var_Na_Ca_exchanger__K_sat * exp(((var_Na_Ca_exchanger__eta - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)))))) * ((exp((var_Na_Ca_exchanger__eta * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Na_i, 3.0) * var_Na_Ca_exchanger__Ca_o) - (exp(((var_Na_Ca_exchanger__eta - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Na_o, 3.0) * var_Na_Ca_exchanger__Ca_i));
        double var_membrane__i_NaCa = var_Na_Ca_exchanger__i_NaCa;
        double var_sodium_potassium_pump__Na_i = var_standard_ionic_concentrations__Na_i;
        double var_sodium_potassium_pump__K_o = var_standard_ionic_concentrations__K_o;
        const double var_sodium_potassium_pump__i_NaK_max = 0.693;
        const double var_sodium_potassium_pump__K_mNai = 10.0;
        double var_sodium_potassium_pump__T = var_membrane__T;
        double var_sodium_potassium_pump__R = var_membrane__R;
        double var_sodium_potassium_pump__Na_o = var_standard_ionic_concentrations__Na_o;
        double var_sodium_potassium_pump__sigma = (1.0 / 7.0) * (exp(var_sodium_potassium_pump__Na_o / 67.3) - 1.0);
        double var_sodium_potassium_pump__V = var_membrane__V;
        double var_sodium_potassium_pump__F = var_membrane__F;
        double var_sodium_potassium_pump__f_NaK = 1.0 / (1.0 + (0.1245 * exp(((-0.1) * var_sodium_potassium_pump__V * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))) + (0.0365 * var_sodium_potassium_pump__sigma * exp(((-var_sodium_potassium_pump__V) * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))));
        const double var_sodium_potassium_pump__K_mKo = 1.5;
        double var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__i_NaK_max * var_sodium_potassium_pump__f_NaK) / (1.0 + pow(var_sodium_potassium_pump__K_mNai / var_sodium_potassium_pump__Na_i, 1.5))) * var_sodium_potassium_pump__K_o) / (var_sodium_potassium_pump__K_o + var_sodium_potassium_pump__K_mKo);
        double var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
        const double var_sarcolemmal_calcium_pump__i_pCa_max = 0.05;
        double var_sarcolemmal_calcium_pump__Ca_i = var_calcium_dynamics__Ca_i;
        const double var_sarcolemmal_calcium_pump__K_mpCa = 0.05;
        double var_sarcolemmal_calcium_pump__i_p_Ca = (var_sarcolemmal_calcium_pump__i_pCa_max * var_sarcolemmal_calcium_pump__Ca_i) / (var_sarcolemmal_calcium_pump__K_mpCa + var_sarcolemmal_calcium_pump__Ca_i);
        double var_membrane__i_p_Ca = var_sarcolemmal_calcium_pump__i_p_Ca;
        const double var_calcium_background_current__g_Cab = 0.0003842;
        double var_calcium_background_current__R = var_membrane__R;
        double var_calcium_background_current__Ca_i = var_calcium_dynamics__Ca_i;
        double var_calcium_background_current__F = var_membrane__F;
        double var_calcium_background_current__Ca_o = var_standard_ionic_concentrations__Ca_o;
        double var_calcium_background_current__T = var_membrane__T;
        double var_calcium_background_current__E_Ca = ((var_calcium_background_current__R * var_calcium_background_current__T) / (2.0 * var_calcium_background_current__F)) * log(var_calcium_background_current__Ca_o / var_calcium_background_current__Ca_i);
        double var_calcium_background_current__V = var_membrane__V;
        double var_calcium_background_current__i_Ca_b = var_calcium_background_current__g_Cab * (var_calcium_background_current__V - var_calcium_background_current__E_Ca);
        double var_membrane__i_Ca_b = var_calcium_background_current__i_Ca_b;
        double var_sodium_background_current__V = var_membrane__V;
        const double var_sodium_background_current__g_Nab = 0.0031;
        double var_sodium_background_current__E_Na = var_fast_sodium_current__E_Na;
        double var_sodium_background_current__i_Na_b = var_sodium_background_current__g_Nab * (var_sodium_background_current__V - var_sodium_background_current__E_Na);
        double var_membrane__i_Na_b = var_sodium_background_current__i_Na_b;

        return var_membrane__i_Na+var_membrane__i_Ca+var_membrane__i_CaK+var_membrane__i_Kr+var_membrane__i_Ks+var_membrane__i_to+var_membrane__i_K1+var_membrane__i_Kp+var_membrane__i_NaCa+var_membrane__i_NaK+var_membrane__i_p_Ca+var_membrane__i_Ca_b+var_membrane__i_Na_b;
    }

    void ComputeResidual(const double rCurrentGuess[3], double rResidual[3])
    {
        std::vector<double>& rY = rGetStateVariables();
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -94.7
        // Units: dimensionless; Initial value: 1
        double var_L_type_Ca_current_f_gate__f = rY[8];
        // Units: dimensionless; Initial value: 0.983
        double var_L_type_Ca_current_d_gate__d = rY[9];
        // Units: dimensionless; Initial value: 0.0001

        double var_L_type_Ca_current_f_Ca_gate__f_Ca = rCurrentGuess[0];
        double var_calcium_dynamics__Ca_SR = rCurrentGuess[1];
        double var_calcium_dynamics__Ca_i = rCurrentGuess[2];

        const double var_membrane__R = 8.314;
        const double var_membrane__T = 310.0;
        const double var_membrane__F = 96.5;
        const double var_standard_ionic_concentrations__Na_i = 10.0;
        const double var_standard_ionic_concentrations__Na_o = 138.0;
        const double var_L_type_Ca_current__P_Ca = 1.26e-05; //2.26e-05;
        double var_L_type_Ca_current__V = var_membrane__V;
        const double var_standard_ionic_concentrations__Ca_o = 2000.0;
        double var_L_type_Ca_current__Ca_o = var_standard_ionic_concentrations__Ca_o;
        double var_L_type_Ca_current__T = var_membrane__T;
        const double var_L_type_Ca_current__C_sc = 1.0;
        double var_L_type_Ca_current__Ca_i = var_calcium_dynamics__Ca_i;
        double var_L_type_Ca_current__R = var_membrane__R;
        double var_L_type_Ca_current__F = var_membrane__F;
        double var_L_type_Ca_current__i_Ca_max = ((((var_L_type_Ca_current__P_Ca / var_L_type_Ca_current__C_sc) * 4.0 * var_L_type_Ca_current__V * pow(var_L_type_Ca_current__F, 2.0)) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) * ((var_L_type_Ca_current__Ca_i * exp((2.0 * var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T))) - (0.341 * var_L_type_Ca_current__Ca_o))) / (exp((2.0 * var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) - 1.0);
        double var_L_type_Ca_current__d = var_L_type_Ca_current_d_gate__d;
        double var_L_type_Ca_current__f = var_L_type_Ca_current_f_gate__f;
        double var_L_type_Ca_current__f_Ca = var_L_type_Ca_current_f_Ca_gate__f_Ca;
        double var_L_type_Ca_current__i_Ca = var_L_type_Ca_current__i_Ca_max * var_L_type_Ca_current__f * var_L_type_Ca_current__d * var_L_type_Ca_current__f_Ca;
        double var_Na_Ca_exchanger__Na_o = var_standard_ionic_concentrations__Na_o;
        const double var_Na_Ca_exchanger__K_NaCa = 1500.0;
        double var_Na_Ca_exchanger__V = var_membrane__V;
        double var_Na_Ca_exchanger__Ca_o = var_standard_ionic_concentrations__Ca_o;
        const double var_Na_Ca_exchanger__eta = 0.35;
        double var_Na_Ca_exchanger__Na_i = var_standard_ionic_concentrations__Na_i;
        double var_Na_Ca_exchanger__T = var_membrane__T;
        double var_Na_Ca_exchanger__R = var_membrane__R;
        const double var_Na_Ca_exchanger__K_sat = 0.2;
        const double var_Na_Ca_exchanger__K_mCa = 1380.0;
        double var_Na_Ca_exchanger__F = var_membrane__F;
        const double var_Na_Ca_exchanger__K_mNa = 87.5;
        double var_Na_Ca_exchanger__Ca_i = var_calcium_dynamics__Ca_i;
        double var_Na_Ca_exchanger__i_NaCa = (var_Na_Ca_exchanger__K_NaCa / ((pow(var_Na_Ca_exchanger__K_mNa, 3.0) + pow(var_Na_Ca_exchanger__Na_o, 3.0)) * (var_Na_Ca_exchanger__K_mCa + var_Na_Ca_exchanger__Ca_o) * (1.0 + (var_Na_Ca_exchanger__K_sat * exp(((var_Na_Ca_exchanger__eta - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)))))) * ((exp((var_Na_Ca_exchanger__eta * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Na_i, 3.0) * var_Na_Ca_exchanger__Ca_o) - (exp(((var_Na_Ca_exchanger__eta - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Na_o, 3.0) * var_Na_Ca_exchanger__Ca_i));
        const double var_sarcolemmal_calcium_pump__i_pCa_max = 0.05;
        double var_sarcolemmal_calcium_pump__Ca_i = var_calcium_dynamics__Ca_i;
        const double var_sarcolemmal_calcium_pump__K_mpCa = 0.05;
        double var_sarcolemmal_calcium_pump__i_p_Ca = (var_sarcolemmal_calcium_pump__i_pCa_max * var_sarcolemmal_calcium_pump__Ca_i) / (var_sarcolemmal_calcium_pump__K_mpCa + var_sarcolemmal_calcium_pump__Ca_i);
        const double var_calcium_background_current__g_Cab = 0.0003842;
        double var_calcium_background_current__R = var_membrane__R;
        double var_calcium_background_current__Ca_i = var_calcium_dynamics__Ca_i;
        double var_calcium_background_current__F = var_membrane__F;
        double var_calcium_background_current__Ca_o = var_standard_ionic_concentrations__Ca_o;
        double var_calcium_background_current__T = var_membrane__T;
        double var_calcium_background_current__E_Ca = ((var_calcium_background_current__R * var_calcium_background_current__T) / (2.0 * var_calcium_background_current__F)) * log(var_calcium_background_current__Ca_o / var_calcium_background_current__Ca_i);
        double var_calcium_background_current__V = var_membrane__V;
        double var_calcium_background_current__i_Ca_b = var_calcium_background_current__g_Cab * (var_calcium_background_current__V - var_calcium_background_current__E_Ca);
        double var_L_type_Ca_current_f_Ca_gate__tau_f_Ca = 30.0;
        double var_L_type_Ca_current_f_Ca_gate__Ca_i = var_L_type_Ca_current__Ca_i;
        const double var_L_type_Ca_current_f_Ca_gate__K_mfCa = 0.18;
        double var_L_type_Ca_current_f_Ca_gate__f_Ca_infinity = 1.0 / (1.0 + pow(var_L_type_Ca_current_f_Ca_gate__Ca_i / var_L_type_Ca_current_f_Ca_gate__K_mfCa, 3.0));
        const double var_calcium_dynamics__K_mCMDN = 2.0;
        const double var_calcium_dynamics__CMDN_tot = 10.0;
        double var_calcium_dynamics__beta_i = 1.0 / (1.0 + ((var_calcium_dynamics__CMDN_tot * var_calcium_dynamics__K_mCMDN) / pow(var_calcium_dynamics__K_mCMDN + var_calcium_dynamics__Ca_i, 2.0)));
        const double var_calcium_dynamics__V_myo = 2.584e-05;
        double var_calcium_dynamics__F = var_membrane__F;
        double var_calcium_dynamics__C_sc = var_L_type_Ca_current__C_sc;
        const double var_calcium_dynamics__A_Cap = 0.0001534;
        double var_calcium_dynamics__d = var_L_type_Ca_current__d;
        double var_calcium_dynamics__V = var_membrane__V;
        const double var_calcium_dynamics__P_rel = 6.0;
        double var_calcium_dynamics__gamma = 1.0 / (1.0 + pow(2000.0 / var_calcium_dynamics__Ca_SR, 3.0));
        double var_calcium_dynamics__f = var_L_type_Ca_current__f;
        double var_calcium_dynamics__f_Ca = var_L_type_Ca_current__f_Ca;
        double var_calcium_dynamics__J_rel = (var_calcium_dynamics__P_rel * var_calcium_dynamics__f * var_calcium_dynamics__d * var_calcium_dynamics__f_Ca * ((var_calcium_dynamics__gamma * var_calcium_dynamics__Ca_SR) - var_calcium_dynamics__Ca_i)) / (1.0 + (1.65 * exp(var_calcium_dynamics__V / 20.0)));
        const double var_calcium_dynamics__P_leak = 1e-06;
        double var_calcium_dynamics__J_leak = var_calcium_dynamics__P_leak * (var_calcium_dynamics__Ca_SR - var_calcium_dynamics__Ca_i);
        const double var_calcium_dynamics__V_up = 0.1;
        const double var_calcium_dynamics__K_mup = 0.32;
        double var_calcium_dynamics__J_up = var_calcium_dynamics__V_up / (1.0 + pow(var_calcium_dynamics__K_mup / var_calcium_dynamics__Ca_i, 2.0));
        double var_calcium_dynamics__i_Ca = var_L_type_Ca_current__i_Ca;
        double var_calcium_dynamics__i_Ca_b = var_calcium_background_current__i_Ca_b;
        double var_calcium_dynamics__i_p_Ca = var_sarcolemmal_calcium_pump__i_p_Ca;
        double var_calcium_dynamics__i_NaCa = var_Na_Ca_exchanger__i_NaCa;
        const double var_calcium_dynamics__K_mCSQN = 600.0;
        const double var_calcium_dynamics__CSQN_tot = 10000.0;
        const double var_calcium_dynamics__V_SR = 2e-06;
        double var_calcium_dynamics__beta_SR = 1.0 / (1.0 + ((var_calcium_dynamics__CSQN_tot * var_calcium_dynamics__K_mCSQN) / pow(var_calcium_dynamics__K_mCSQN + var_calcium_dynamics__Ca_SR, 2.0)));
        double d_dt_L_type_Ca_current_f_Ca_gate__f_Ca = (var_L_type_Ca_current_f_Ca_gate__f_Ca_infinity - var_L_type_Ca_current_f_Ca_gate__f_Ca) / var_L_type_Ca_current_f_Ca_gate__tau_f_Ca;
        double d_dt_calcium_dynamics__Ca_SR = (var_calcium_dynamics__beta_SR * ((var_calcium_dynamics__J_up - var_calcium_dynamics__J_leak) - var_calcium_dynamics__J_rel) * var_calcium_dynamics__V_myo) / var_calcium_dynamics__V_SR;
        double d_dt_calcium_dynamics__Ca_i = var_calcium_dynamics__beta_i * (((var_calcium_dynamics__J_rel + var_calcium_dynamics__J_leak) - var_calcium_dynamics__J_up) - (((var_calcium_dynamics__A_Cap * var_calcium_dynamics__C_sc) / (2.0 * var_calcium_dynamics__F * var_calcium_dynamics__V_myo)) * ((var_calcium_dynamics__i_Ca + var_calcium_dynamics__i_Ca_b + var_calcium_dynamics__i_p_Ca) - (2.0 * var_calcium_dynamics__i_NaCa))));

        rResidual[0] = rCurrentGuess[0] - rY[10] - mDt*d_dt_L_type_Ca_current_f_Ca_gate__f_Ca;
        rResidual[2] = rCurrentGuess[2] - rY[11] - mDt*d_dt_calcium_dynamics__Ca_i;
        rResidual[1] = rCurrentGuess[1] - rY[12] - mDt*d_dt_calcium_dynamics__Ca_SR;
    }

    void ComputeJacobian(const double rCurrentGuess[3], double rJacobian[3][3])
    {
        std::vector<double>& rY = rGetStateVariables();
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -94.7
        double var_L_type_Ca_current_f_gate__f = rY[8];
        // Units: dimensionless; Initial value: 0.983
        double var_L_type_Ca_current_d_gate__d = rY[9];
        // Units: dimensionless; Initial value: 0.0001

        double var_L_type_Ca_current_f_Ca_gate__f_Ca = rCurrentGuess[0];
        double var_calcium_dynamics__Ca_SR = rCurrentGuess[1];
        double var_calcium_dynamics__Ca_i = rCurrentGuess[2];

        const double var_membrane__R = 8.314;
        const double var_membrane__T = 310.0;
        const double var_membrane__F = 96.5;
        const double var_standard_ionic_concentrations__Na_i = 10.0;
        const double var_standard_ionic_concentrations__Na_o = 138.0;
        const double var_L_type_Ca_current__P_Ca = 1.26e-05; //2.26e-05;
        const double var_standard_ionic_concentrations__Ca_o = 2000.0;
        const double var_L_type_Ca_current__C_sc = 1.0;
        const double var_Na_Ca_exchanger__K_NaCa = 1500.0;
        const double var_Na_Ca_exchanger__eta = 0.35;
        const double var_Na_Ca_exchanger__K_sat = 0.2;
        const double var_Na_Ca_exchanger__K_mCa = 1380.0;
        const double var_Na_Ca_exchanger__K_mNa = 87.5;
        const double var_sarcolemmal_calcium_pump__i_pCa_max = 0.05;
        const double var_sarcolemmal_calcium_pump__K_mpCa = 0.05;
        const double var_calcium_background_current__g_Cab = 0.0003842;
        const double var_L_type_Ca_current_f_Ca_gate__K_mfCa = 0.18;
        const double var_calcium_dynamics__K_mCMDN = 2.0;
        const double var_calcium_dynamics__CMDN_tot = 10.0;
        const double var_calcium_dynamics__V_myo = 2.584e-05;
        const double var_calcium_dynamics__A_Cap = 0.0001534;
        const double var_calcium_dynamics__P_rel = 6.0;
        const double var_calcium_dynamics__P_leak = 1e-06;
        const double var_calcium_dynamics__V_up = 0.1;
        const double var_calcium_dynamics__K_mup = 0.32;
        const double var_calcium_dynamics__K_mCSQN = 600.0;
        const double var_calcium_dynamics__CSQN_tot = 10000.0;
        const double var_calcium_dynamics__V_SR = 2e-06;

        rJacobian[0][0] = 1.0 + ((1.0 / 30.0) * mDt);
        rJacobian[0][1] = 0.0;
        rJacobian[0][2] = ((((1.0 / 10.0) * mDt) / pow(1.0 + (pow(var_calcium_dynamics__Ca_i, 3.0) / pow(var_L_type_Ca_current_f_Ca_gate__K_mfCa, 3.0)), 2.0)) * pow(var_calcium_dynamics__Ca_i, 2.0)) / pow(var_L_type_Ca_current_f_Ca_gate__K_mfCa, 3.0);
        rJacobian[1][0] = (((((((mDt / (1.0 + ((var_calcium_dynamics__CSQN_tot * var_calcium_dynamics__K_mCSQN) / pow(var_calcium_dynamics__K_mCSQN + var_calcium_dynamics__Ca_SR, 2.0)))) * var_calcium_dynamics__P_rel) * var_L_type_Ca_current_f_gate__f) * var_L_type_Ca_current_d_gate__d) * (((1.0 / (1.0 + (8000000000.0 / pow(var_calcium_dynamics__Ca_SR, 3.0)))) * var_calcium_dynamics__Ca_SR) - var_calcium_dynamics__Ca_i)) / (1.0 + (1.65 * exp((1.0 / 20.0) * var_membrane__V)))) * var_calcium_dynamics__V_myo) / var_calcium_dynamics__V_SR;
        rJacobian[1][1] = (1.0 - ((((((((2.0 * mDt) / pow(1.0 + ((var_calcium_dynamics__CSQN_tot * var_calcium_dynamics__K_mCSQN) / pow(var_calcium_dynamics__K_mCSQN + var_calcium_dynamics__Ca_SR, 2.0)), 2.0)) * (((var_calcium_dynamics__V_up / (1.0 + (pow(var_calcium_dynamics__K_mup, 2.0) / pow(var_calcium_dynamics__Ca_i, 2.0)))) - (var_calcium_dynamics__P_leak * (var_calcium_dynamics__Ca_SR - var_calcium_dynamics__Ca_i))) - (((((var_calcium_dynamics__P_rel * var_L_type_Ca_current_f_gate__f) * var_L_type_Ca_current_d_gate__d) * var_L_type_Ca_current_f_Ca_gate__f_Ca) * (((1.0 / (1.0 + (8000000000.0 / pow(var_calcium_dynamics__Ca_SR, 3.0)))) * var_calcium_dynamics__Ca_SR) - var_calcium_dynamics__Ca_i)) / (1.0 + (1.65 * exp((1.0 / 20.0) * var_membrane__V)))))) * var_calcium_dynamics__V_myo) / var_calcium_dynamics__V_SR) * var_calcium_dynamics__CSQN_tot) * var_calcium_dynamics__K_mCSQN) / pow(var_calcium_dynamics__K_mCSQN + var_calcium_dynamics__Ca_SR, 3.0))) - ((((mDt / (1.0 + ((var_calcium_dynamics__CSQN_tot * var_calcium_dynamics__K_mCSQN) / pow(var_calcium_dynamics__K_mCSQN + var_calcium_dynamics__Ca_SR, 2.0)))) * ((-var_calcium_dynamics__P_leak) - (((((var_calcium_dynamics__P_rel * var_L_type_Ca_current_f_gate__f) * var_L_type_Ca_current_d_gate__d) * var_L_type_Ca_current_f_Ca_gate__f_Ca) * (((24000000000.0 / pow(1.0 + (8000000000.0 / pow(var_calcium_dynamics__Ca_SR, 3.0)), 2.0)) / pow(var_calcium_dynamics__Ca_SR, 3.0)) + (1.0 / (1.0 + (8000000000.0 / pow(var_calcium_dynamics__Ca_SR, 3.0)))))) / (1.0 + (1.65 * exp((1.0 / 20.0) * var_membrane__V)))))) * var_calcium_dynamics__V_myo) / var_calcium_dynamics__V_SR);
        rJacobian[1][2] = ((((-mDt) / (1.0 + ((var_calcium_dynamics__CSQN_tot * var_calcium_dynamics__K_mCSQN) / pow(var_calcium_dynamics__K_mCSQN + var_calcium_dynamics__Ca_SR, 2.0)))) * ((((((2.0 * var_calcium_dynamics__V_up) / pow(1.0 + (pow(var_calcium_dynamics__K_mup, 2.0) / pow(var_calcium_dynamics__Ca_i, 2.0)), 2.0)) * pow(var_calcium_dynamics__K_mup, 2.0)) / pow(var_calcium_dynamics__Ca_i, 3.0)) + var_calcium_dynamics__P_leak) + ((((var_calcium_dynamics__P_rel * var_L_type_Ca_current_f_gate__f) * var_L_type_Ca_current_d_gate__d) * var_L_type_Ca_current_f_Ca_gate__f_Ca) / (1.0 + (1.65 * exp((1.0 / 20.0) * var_membrane__V)))))) * var_calcium_dynamics__V_myo) / var_calcium_dynamics__V_SR;
        rJacobian[2][0] = ((-mDt) / (1.0 + ((var_calcium_dynamics__CMDN_tot * var_calcium_dynamics__K_mCMDN) / pow(var_calcium_dynamics__K_mCMDN + var_calcium_dynamics__Ca_i, 2.0)))) * (((((var_calcium_dynamics__P_rel * var_L_type_Ca_current_f_gate__f) * var_L_type_Ca_current_d_gate__d) * (((1.0 / (1.0 + (8000000000.0 / pow(var_calcium_dynamics__Ca_SR, 3.0)))) * var_calcium_dynamics__Ca_SR) - var_calcium_dynamics__Ca_i)) / (1.0 + (1.65 * exp((1.0 / 20.0) * var_membrane__V)))) - (((((((((((2.0 * var_calcium_dynamics__A_Cap) * var_membrane__F) / var_calcium_dynamics__V_myo) * var_L_type_Ca_current__P_Ca) * var_membrane__V) / var_membrane__R) / var_membrane__T) * ((var_calcium_dynamics__Ca_i * exp((((2.0 * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T)) - (0.341 * var_standard_ionic_concentrations__Ca_o))) / (exp((((2.0 * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) - 1.0)) * var_L_type_Ca_current_f_gate__f) * var_L_type_Ca_current_d_gate__d));
        rJacobian[2][1] = ((-mDt) / (1.0 + ((var_calcium_dynamics__CMDN_tot * var_calcium_dynamics__K_mCMDN) / pow(var_calcium_dynamics__K_mCMDN + var_calcium_dynamics__Ca_i, 2.0)))) * ((((((var_calcium_dynamics__P_rel * var_L_type_Ca_current_f_gate__f) * var_L_type_Ca_current_d_gate__d) * var_L_type_Ca_current_f_Ca_gate__f_Ca) * (((24000000000.0 / pow(1.0 + (8000000000.0 / pow(var_calcium_dynamics__Ca_SR, 3.0)), 2.0)) / pow(var_calcium_dynamics__Ca_SR, 3.0)) + (1.0 / (1.0 + (8000000000.0 / pow(var_calcium_dynamics__Ca_SR, 3.0)))))) / (1.0 + (1.65 * exp((1.0 / 20.0) * var_membrane__V)))) + var_calcium_dynamics__P_leak);
        rJacobian[2][2] = (1.0 - ((((((2.0 * mDt) / pow(1.0 + ((var_calcium_dynamics__CMDN_tot * var_calcium_dynamics__K_mCMDN) / pow(var_calcium_dynamics__K_mCMDN + var_calcium_dynamics__Ca_i, 2.0)), 2.0)) * ((((((((var_calcium_dynamics__P_rel * var_L_type_Ca_current_f_gate__f) * var_L_type_Ca_current_d_gate__d) * var_L_type_Ca_current_f_Ca_gate__f_Ca) * (((1.0 / (1.0 + (8000000000.0 / pow(var_calcium_dynamics__Ca_SR, 3.0)))) * var_calcium_dynamics__Ca_SR) - var_calcium_dynamics__Ca_i)) / (1.0 + (1.65 * exp((1.0 / 20.0) * var_membrane__V)))) + (var_calcium_dynamics__P_leak * (var_calcium_dynamics__Ca_SR - var_calcium_dynamics__Ca_i))) - (var_calcium_dynamics__V_up / (1.0 + (pow(var_calcium_dynamics__K_mup, 2.0) / pow(var_calcium_dynamics__Ca_i, 2.0))))) - ((((((1.0 / 2.0) * var_calcium_dynamics__A_Cap) * var_L_type_Ca_current__C_sc) / var_membrane__F) / var_calcium_dynamics__V_myo) * ((((((((((((((4.0 * var_L_type_Ca_current__P_Ca) / var_L_type_Ca_current__C_sc) * var_membrane__V) * pow(var_membrane__F, 2.0)) / var_membrane__R) / var_membrane__T) * ((var_calcium_dynamics__Ca_i * exp((((2.0 * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T)) - (0.341 * var_standard_ionic_concentrations__Ca_o))) / (exp((((2.0 * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) - 1.0)) * var_L_type_Ca_current_f_gate__f) * var_L_type_Ca_current_d_gate__d) * var_L_type_Ca_current_f_Ca_gate__f_Ca) + (var_calcium_background_current__g_Cab * (var_membrane__V - (((((1.0 / 2.0) * var_membrane__R) * var_membrane__T) / var_membrane__F) * log(var_standard_ionic_concentrations__Ca_o / var_calcium_dynamics__Ca_i))))) + ((var_sarcolemmal_calcium_pump__i_pCa_max * var_calcium_dynamics__Ca_i) / (var_sarcolemmal_calcium_pump__K_mpCa + var_calcium_dynamics__Ca_i))) - (((((2.0 * var_Na_Ca_exchanger__K_NaCa) / (pow(var_Na_Ca_exchanger__K_mNa, 3.0) + pow(var_standard_ionic_concentrations__Na_o, 3.0))) / (var_Na_Ca_exchanger__K_mCa + var_standard_ionic_concentrations__Ca_o)) / (1.0 + (var_Na_Ca_exchanger__K_sat * exp(((((var_Na_Ca_exchanger__eta - 1.0) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T)))) * (((exp((((var_Na_Ca_exchanger__eta * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_standard_ionic_concentrations__Na_i, 3.0)) * var_standard_ionic_concentrations__Ca_o) - ((exp(((((var_Na_Ca_exchanger__eta - 1.0) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) * pow(var_standard_ionic_concentrations__Na_o, 3.0)) * var_calcium_dynamics__Ca_i))))))) * var_calcium_dynamics__CMDN_tot) * var_calcium_dynamics__K_mCMDN) / pow(var_calcium_dynamics__K_mCMDN + var_calcium_dynamics__Ca_i, 3.0))) - ((mDt / (1.0 + ((var_calcium_dynamics__CMDN_tot * var_calcium_dynamics__K_mCMDN) / pow(var_calcium_dynamics__K_mCMDN + var_calcium_dynamics__Ca_i, 2.0)))) * ((((((((-var_calcium_dynamics__P_rel) * var_L_type_Ca_current_f_gate__f) * var_L_type_Ca_current_d_gate__d) * var_L_type_Ca_current_f_Ca_gate__f_Ca) / (1.0 + (1.65 * exp((1.0 / 20.0) * var_membrane__V)))) - var_calcium_dynamics__P_leak) - ((((2.0 * var_calcium_dynamics__V_up) / pow(1.0 + (pow(var_calcium_dynamics__K_mup, 2.0) / pow(var_calcium_dynamics__Ca_i, 2.0)), 2.0)) * pow(var_calcium_dynamics__K_mup, 2.0)) / pow(var_calcium_dynamics__Ca_i, 3.0))) - ((((((1.0 / 2.0) * var_calcium_dynamics__A_Cap) * var_L_type_Ca_current__C_sc) / var_membrane__F) / var_calcium_dynamics__V_myo) * (((((((((((((((4.0 * var_L_type_Ca_current__P_Ca) / var_L_type_Ca_current__C_sc) * var_membrane__V) * pow(var_membrane__F, 2.0)) / var_membrane__R) / var_membrane__T) * exp((((2.0 * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T)) / (exp((((2.0 * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T) - 1.0)) * var_L_type_Ca_current_f_gate__f) * var_L_type_Ca_current_d_gate__d) * var_L_type_Ca_current_f_Ca_gate__f_Ca) + ((((((1.0 / 2.0) * var_calcium_background_current__g_Cab) * var_membrane__R) * var_membrane__T) / var_membrane__F) / var_calcium_dynamics__Ca_i)) + (var_sarcolemmal_calcium_pump__i_pCa_max / (var_sarcolemmal_calcium_pump__K_mpCa + var_calcium_dynamics__Ca_i))) - ((var_sarcolemmal_calcium_pump__i_pCa_max * var_calcium_dynamics__Ca_i) / pow(var_sarcolemmal_calcium_pump__K_mpCa + var_calcium_dynamics__Ca_i, 2.0))) + ((((((2.0 * var_Na_Ca_exchanger__K_NaCa) / (pow(var_Na_Ca_exchanger__K_mNa, 3.0) + pow(var_standard_ionic_concentrations__Na_o, 3.0))) / (var_Na_Ca_exchanger__K_mCa + var_standard_ionic_concentrations__Ca_o)) / (1.0 + (var_Na_Ca_exchanger__K_sat * exp(((((var_Na_Ca_exchanger__eta - 1.0) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T)))) * exp(((((var_Na_Ca_exchanger__eta - 1.0) * var_membrane__V) * var_membrane__F) / var_membrane__R) / var_membrane__T)) * pow(var_standard_ionic_concentrations__Na_o, 3.0))))));
    }

protected:
    void UpdateTransmembranePotential(double var_environment__time)
    {
        // Time units: millisecond
        std::vector<double>& rY = rGetStateVariables();
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -94.7
        double var_fast_sodium_current_m_gate__m = rY[1];
        // Units: dimensionless; Initial value: 0.00024676
        double var_fast_sodium_current_h_gate__h = rY[2];
        // Units: dimensionless; Initial value: 0.99869
        double var_fast_sodium_current_j_gate__j = rY[3];
        // Units: dimensionless; Initial value: 0.99887
        double var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr = rY[4];
        // Units: dimensionless; Initial value: 0.229
        double var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks = rY[5];
        // Units: dimensionless; Initial value: 0.0001
        double var_transient_outward_potassium_current_X_to_gate__X_to = rY[6];
        // Units: dimensionless; Initial value: 0.00003742
        double var_transient_outward_potassium_current_Y_to_gate__Y_to = rY[7];
        // Units: dimensionless; Initial value: 1
        double var_L_type_Ca_current_f_gate__f = rY[8];
        // Units: dimensionless; Initial value: 0.983
        double var_L_type_Ca_current_d_gate__d = rY[9];
        // Units: dimensionless; Initial value: 0.0001
        double var_L_type_Ca_current_f_Ca_gate__f_Ca = rY[10];
        // Units: dimensionless; Initial value: 0.942
        double var_calcium_dynamics__Ca_i = rY[11];
        // Units: micromolar; Initial value: 0.0472

        const double var_membrane__R = 8.314;
        const double var_membrane__T = 310.0;
        const double var_membrane__F = 96.5;
        double var_membrane__time = var_environment__time;
        double var_fast_sodium_current__h = var_fast_sodium_current_h_gate__h;
        const double var_fast_sodium_current__g_Na = 12.8;
        double var_fast_sodium_current__j = var_fast_sodium_current_j_gate__j;
        double var_fast_sodium_current__T = var_membrane__T;
        double var_fast_sodium_current__R = var_membrane__R;
        const double var_standard_ionic_concentrations__Na_i = 10.0;
        double var_fast_sodium_current__Na_i = var_standard_ionic_concentrations__Na_i;
        const double var_standard_ionic_concentrations__Na_o = 138.0;
        double var_fast_sodium_current__Na_o = var_standard_ionic_concentrations__Na_o;
        double var_fast_sodium_current__F = var_membrane__F;
        double var_fast_sodium_current__E_Na = ((var_fast_sodium_current__R * var_fast_sodium_current__T) / var_fast_sodium_current__F) * log(var_fast_sodium_current__Na_o / var_fast_sodium_current__Na_i);
        double var_fast_sodium_current__V = var_membrane__V;
        double var_fast_sodium_current__m = var_fast_sodium_current_m_gate__m;
        double var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na);
        double var_membrane__i_Na = var_fast_sodium_current__i_Na;
        const double var_L_type_Ca_current__P_Ca = 1.26e-05; //2.26e-05;
        double var_L_type_Ca_current__V = var_membrane__V;
        const double var_standard_ionic_concentrations__Ca_o = 2000.0;
        double var_L_type_Ca_current__Ca_o = var_standard_ionic_concentrations__Ca_o;
        double var_L_type_Ca_current__T = var_membrane__T;
        const double var_L_type_Ca_current__C_sc = 1.0;
        double var_L_type_Ca_current__Ca_i = var_calcium_dynamics__Ca_i;
        double var_L_type_Ca_current__R = var_membrane__R;
        double var_L_type_Ca_current__F = var_membrane__F;
        double var_L_type_Ca_current__i_Ca_max = ((((var_L_type_Ca_current__P_Ca / var_L_type_Ca_current__C_sc) * 4.0 * var_L_type_Ca_current__V * pow(var_L_type_Ca_current__F, 2.0)) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) * ((var_L_type_Ca_current__Ca_i * exp((2.0 * var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T))) - (0.341 * var_L_type_Ca_current__Ca_o))) / (exp((2.0 * var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) - 1.0);
        double var_L_type_Ca_current__d = var_L_type_Ca_current_d_gate__d;
        double var_L_type_Ca_current__f = var_L_type_Ca_current_f_gate__f;
        double var_L_type_Ca_current__f_Ca = var_L_type_Ca_current_f_Ca_gate__f_Ca;
        double var_L_type_Ca_current__i_Ca = var_L_type_Ca_current__i_Ca_max * var_L_type_Ca_current__f * var_L_type_Ca_current__d * var_L_type_Ca_current__f_Ca;
        double var_membrane__i_Ca = var_L_type_Ca_current__i_Ca;
        const double var_standard_ionic_concentrations__K_o = 4.0;
        double var_L_type_Ca_current__K_o = var_standard_ionic_concentrations__K_o;
        const double var_L_type_Ca_current__i_Ca_half =  -0.265;
        const double var_L_type_Ca_current__P_CaK = 5.79e-07;
        const double var_standard_ionic_concentrations__K_i = 149.4;
        double var_L_type_Ca_current__K_i = var_standard_ionic_concentrations__K_i;
        double var_L_type_Ca_current__i_CaK = ((((((var_L_type_Ca_current__P_CaK / var_L_type_Ca_current__C_sc) * var_L_type_Ca_current__f * var_L_type_Ca_current__d * var_L_type_Ca_current__f_Ca) / (1.0 + (var_L_type_Ca_current__i_Ca_max / var_L_type_Ca_current__i_Ca_half))) * 1000.0 * var_L_type_Ca_current__V * pow(var_L_type_Ca_current__F, 2.0)) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) * ((var_L_type_Ca_current__K_i * exp((var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T))) - var_L_type_Ca_current__K_o)) / (exp((var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) - 1.0);
        double var_membrane__i_CaK = var_L_type_Ca_current__i_CaK;
        double var_rapid_activating_delayed_rectifiyer_K_current__V = var_membrane__V;
        double var_rapid_activating_delayed_rectifiyer_K_current__R_V = 1.0 / (1.0 + (2.5 * exp(0.1 * (var_rapid_activating_delayed_rectifiyer_K_current__V + 28.0))));
        double var_rapid_activating_delayed_rectifiyer_K_current__X_kr = var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr;
        double var_rapid_activating_delayed_rectifiyer_K_current__R = var_membrane__R;
        double var_rapid_activating_delayed_rectifiyer_K_current__K_o = var_standard_ionic_concentrations__K_o;
        double var_rapid_activating_delayed_rectifiyer_K_current__T = var_membrane__T;
        double var_rapid_activating_delayed_rectifiyer_K_current__K_i = var_standard_ionic_concentrations__K_i;
        double var_rapid_activating_delayed_rectifiyer_K_current__F = var_membrane__F;
        double var_rapid_activating_delayed_rectifiyer_K_current__E_K = ((var_rapid_activating_delayed_rectifiyer_K_current__R * var_rapid_activating_delayed_rectifiyer_K_current__T) / var_rapid_activating_delayed_rectifiyer_K_current__F) * log(var_rapid_activating_delayed_rectifiyer_K_current__K_o / var_rapid_activating_delayed_rectifiyer_K_current__K_i);
        const double var_rapid_activating_delayed_rectifiyer_K_current__g_Kr = 0.0136;
        double var_rapid_activating_delayed_rectifiyer_K_current__i_Kr = var_rapid_activating_delayed_rectifiyer_K_current__g_Kr * var_rapid_activating_delayed_rectifiyer_K_current__R_V * var_rapid_activating_delayed_rectifiyer_K_current__X_kr * sqrt(var_rapid_activating_delayed_rectifiyer_K_current__K_o / 4.0) * (var_rapid_activating_delayed_rectifiyer_K_current__V - var_rapid_activating_delayed_rectifiyer_K_current__E_K);
        double var_membrane__i_Kr = var_rapid_activating_delayed_rectifiyer_K_current__i_Kr;
        const double var_slow_activating_delayed_rectifiyer_K_current__g_Ks = 0.0245;
        double var_slow_activating_delayed_rectifiyer_K_current__V = var_membrane__V;
        double var_slow_activating_delayed_rectifiyer_K_current__K_i = var_standard_ionic_concentrations__K_i;
        double var_slow_activating_delayed_rectifiyer_K_current__F = var_membrane__F;
        double var_slow_activating_delayed_rectifiyer_K_current__K_o = var_standard_ionic_concentrations__K_o;
        double var_slow_activating_delayed_rectifiyer_K_current__T = var_membrane__T;
        double var_slow_activating_delayed_rectifiyer_K_current__Na_o = var_standard_ionic_concentrations__Na_o;
        double var_slow_activating_delayed_rectifiyer_K_current__R = var_membrane__R;
        double var_slow_activating_delayed_rectifiyer_K_current__Na_i = var_standard_ionic_concentrations__Na_i;
        double var_slow_activating_delayed_rectifiyer_K_current__E_Ks = ((var_slow_activating_delayed_rectifiyer_K_current__R * var_slow_activating_delayed_rectifiyer_K_current__T) / var_slow_activating_delayed_rectifiyer_K_current__F) * log((var_slow_activating_delayed_rectifiyer_K_current__K_o + (0.01833 * var_slow_activating_delayed_rectifiyer_K_current__Na_o)) / (var_slow_activating_delayed_rectifiyer_K_current__K_i + (0.01833 * var_slow_activating_delayed_rectifiyer_K_current__Na_i)));
        double var_slow_activating_delayed_rectifiyer_K_current__X_ks = var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks;
        double var_slow_activating_delayed_rectifiyer_K_current__i_Ks = var_slow_activating_delayed_rectifiyer_K_current__g_Ks * pow(var_slow_activating_delayed_rectifiyer_K_current__X_ks, 2.0) * (var_slow_activating_delayed_rectifiyer_K_current__V - var_slow_activating_delayed_rectifiyer_K_current__E_Ks);
        double var_membrane__i_Ks = var_slow_activating_delayed_rectifiyer_K_current__i_Ks;
        const double var_transient_outward_potassium_current__g_to = 0.23815;
        double var_transient_outward_potassium_current__V = var_membrane__V;
        double var_transient_outward_potassium_current__X_to = var_transient_outward_potassium_current_X_to_gate__X_to;
        double var_transient_outward_potassium_current__E_K = var_rapid_activating_delayed_rectifiyer_K_current__E_K;
        double var_transient_outward_potassium_current__Y_to = var_transient_outward_potassium_current_Y_to_gate__Y_to;
        double var_transient_outward_potassium_current__i_to = var_transient_outward_potassium_current__g_to * var_transient_outward_potassium_current__X_to * var_transient_outward_potassium_current__Y_to * (var_transient_outward_potassium_current__V - var_transient_outward_potassium_current__E_K);
        double var_membrane__i_to = var_transient_outward_potassium_current__i_to;
        const double var_time_independent_potassium_current__g_K1 = 2.8;
        double var_time_independent_potassium_current__F = var_membrane__F;
        double var_time_independent_potassium_current_K1_gate__F = var_time_independent_potassium_current__F;
        double var_time_independent_potassium_current__V = var_membrane__V;
        double var_time_independent_potassium_current_K1_gate__V = var_time_independent_potassium_current__V;
        double var_time_independent_potassium_current__T = var_membrane__T;
        double var_time_independent_potassium_current_K1_gate__T = var_time_independent_potassium_current__T;
        double var_time_independent_potassium_current__R = var_membrane__R;
        double var_time_independent_potassium_current_K1_gate__R = var_time_independent_potassium_current__R;
        double var_time_independent_potassium_current__E_K = var_rapid_activating_delayed_rectifiyer_K_current__E_K;
        double var_time_independent_potassium_current_K1_gate__E_K = var_time_independent_potassium_current__E_K;
        double var_time_independent_potassium_current_K1_gate__K1_infinity = 1.0 / (2.0 + exp(((1.62 * var_time_independent_potassium_current_K1_gate__F) / (var_time_independent_potassium_current_K1_gate__R * var_time_independent_potassium_current_K1_gate__T)) * (var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K)));
        double var_time_independent_potassium_current__K1_infinity = var_time_independent_potassium_current_K1_gate__K1_infinity;
        double var_time_independent_potassium_current__K_o = var_standard_ionic_concentrations__K_o;
        const double var_time_independent_potassium_current__K_mK1 = 13.0;
        double var_time_independent_potassium_current__i_K1 = ((var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K1_infinity * var_time_independent_potassium_current__K_o) / (var_time_independent_potassium_current__K_o + var_time_independent_potassium_current__K_mK1)) * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
        double var_plateau_potassium_current__V = var_membrane__V;
        double var_plateau_potassium_current__E_K = var_rapid_activating_delayed_rectifiyer_K_current__E_K;
        const double var_plateau_potassium_current__g_Kp = 0.002216;
        double var_plateau_potassium_current_Kp_gate__V = var_plateau_potassium_current__V;
        double var_plateau_potassium_current_Kp_gate__Kp_V = 1.0 / (1.0 + exp((7.488 - var_plateau_potassium_current_Kp_gate__V) / 5.98));
        double var_plateau_potassium_current__Kp_V = var_plateau_potassium_current_Kp_gate__Kp_V;
        double var_plateau_potassium_current__i_Kp = var_plateau_potassium_current__g_Kp * var_plateau_potassium_current__Kp_V * (var_plateau_potassium_current__V - var_plateau_potassium_current__E_K);
        double var_membrane__i_Kp = var_plateau_potassium_current__i_Kp;
        double var_Na_Ca_exchanger__Na_o = var_standard_ionic_concentrations__Na_o;
        const double var_Na_Ca_exchanger__K_NaCa = 1500.0;
        double var_Na_Ca_exchanger__V = var_membrane__V;
        double var_Na_Ca_exchanger__Ca_o = var_standard_ionic_concentrations__Ca_o;
        const double var_Na_Ca_exchanger__eta = 0.35;
        double var_Na_Ca_exchanger__Na_i = var_standard_ionic_concentrations__Na_i;
        double var_Na_Ca_exchanger__T = var_membrane__T;
        double var_Na_Ca_exchanger__R = var_membrane__R;
        const double var_Na_Ca_exchanger__K_sat = 0.2;
        const double var_Na_Ca_exchanger__K_mCa = 1380.0;
        double var_Na_Ca_exchanger__F = var_membrane__F;
        const double var_Na_Ca_exchanger__K_mNa = 87.5;
        double var_Na_Ca_exchanger__Ca_i = var_calcium_dynamics__Ca_i;
        double var_Na_Ca_exchanger__i_NaCa = (var_Na_Ca_exchanger__K_NaCa / ((pow(var_Na_Ca_exchanger__K_mNa, 3.0) + pow(var_Na_Ca_exchanger__Na_o, 3.0)) * (var_Na_Ca_exchanger__K_mCa + var_Na_Ca_exchanger__Ca_o) * (1.0 + (var_Na_Ca_exchanger__K_sat * exp(((var_Na_Ca_exchanger__eta - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)))))) * ((exp((var_Na_Ca_exchanger__eta * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Na_i, 3.0) * var_Na_Ca_exchanger__Ca_o) - (exp(((var_Na_Ca_exchanger__eta - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Na_o, 3.0) * var_Na_Ca_exchanger__Ca_i));
        double var_membrane__i_NaCa = var_Na_Ca_exchanger__i_NaCa;
        double var_sodium_potassium_pump__Na_i = var_standard_ionic_concentrations__Na_i;
        double var_sodium_potassium_pump__K_o = var_standard_ionic_concentrations__K_o;
        const double var_sodium_potassium_pump__i_NaK_max = 0.693;
        const double var_sodium_potassium_pump__K_mNai = 10.0;
        double var_sodium_potassium_pump__T = var_membrane__T;
        double var_sodium_potassium_pump__R = var_membrane__R;
        double var_sodium_potassium_pump__Na_o = var_standard_ionic_concentrations__Na_o;
        double var_sodium_potassium_pump__sigma = (1.0 / 7.0) * (exp(var_sodium_potassium_pump__Na_o / 67.3) - 1.0);
        double var_sodium_potassium_pump__V = var_membrane__V;
        double var_sodium_potassium_pump__F = var_membrane__F;
        double var_sodium_potassium_pump__f_NaK = 1.0 / (1.0 + (0.1245 * exp(((-0.1) * var_sodium_potassium_pump__V * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))) + (0.0365 * var_sodium_potassium_pump__sigma * exp(((-var_sodium_potassium_pump__V) * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))));
        const double var_sodium_potassium_pump__K_mKo = 1.5;
        double var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__i_NaK_max * var_sodium_potassium_pump__f_NaK) / (1.0 + pow(var_sodium_potassium_pump__K_mNai / var_sodium_potassium_pump__Na_i, 1.5))) * var_sodium_potassium_pump__K_o) / (var_sodium_potassium_pump__K_o + var_sodium_potassium_pump__K_mKo);
        double var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
        const double var_sarcolemmal_calcium_pump__i_pCa_max = 0.05;
        double var_sarcolemmal_calcium_pump__Ca_i = var_calcium_dynamics__Ca_i;
        const double var_sarcolemmal_calcium_pump__K_mpCa = 0.05;
        double var_sarcolemmal_calcium_pump__i_p_Ca = (var_sarcolemmal_calcium_pump__i_pCa_max * var_sarcolemmal_calcium_pump__Ca_i) / (var_sarcolemmal_calcium_pump__K_mpCa + var_sarcolemmal_calcium_pump__Ca_i);
        double var_membrane__i_p_Ca = var_sarcolemmal_calcium_pump__i_p_Ca;
        const double var_calcium_background_current__g_Cab = 0.0003842;
        double var_calcium_background_current__R = var_membrane__R;
        double var_calcium_background_current__Ca_i = var_calcium_dynamics__Ca_i;
        double var_calcium_background_current__F = var_membrane__F;
        double var_calcium_background_current__Ca_o = var_standard_ionic_concentrations__Ca_o;
        double var_calcium_background_current__T = var_membrane__T;
        double var_calcium_background_current__E_Ca = ((var_calcium_background_current__R * var_calcium_background_current__T) / (2.0 * var_calcium_background_current__F)) * log(var_calcium_background_current__Ca_o / var_calcium_background_current__Ca_i);
        double var_calcium_background_current__V = var_membrane__V;
        double var_calcium_background_current__i_Ca_b = var_calcium_background_current__g_Cab * (var_calcium_background_current__V - var_calcium_background_current__E_Ca);
        double var_membrane__i_Ca_b = var_calcium_background_current__i_Ca_b;
        double var_sodium_background_current__V = var_membrane__V;
        const double var_sodium_background_current__g_Nab = 0.0031;
        double var_sodium_background_current__E_Na = var_fast_sodium_current__E_Na;
        double var_sodium_background_current__i_Na_b = var_sodium_background_current__g_Nab * (var_sodium_background_current__V - var_sodium_background_current__E_Na);
        double var_membrane__i_Na_b = var_sodium_background_current__i_Na_b;
        double var_membrane__i_Stim = GetStimulus(var_membrane__time);
        double d_dt_membrane__V = -(var_membrane__i_Na + var_membrane__i_Ca + var_membrane__i_CaK + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_to + var_membrane__i_K1 + var_membrane__i_Kp + var_membrane__i_NaCa + var_membrane__i_NaK + var_membrane__i_p_Ca + var_membrane__i_Na_b + var_membrane__i_Ca_b + var_membrane__i_Stim);

        rY[0] += mDt * d_dt_membrane__V;
    }

    void ComputeOneStepExceptVoltage(double var_environment__time)
    {
        // Time units: millisecond
        std::vector<double>& rY = rGetStateVariables();
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -94.7

        double var_fast_sodium_current__V = var_membrane__V;
        double var_L_type_Ca_current__V = var_membrane__V;
        double var_rapid_activating_delayed_rectifiyer_K_current__V = var_membrane__V;
        double var_slow_activating_delayed_rectifiyer_K_current__V = var_membrane__V;
        double var_transient_outward_potassium_current__V = var_membrane__V;
        double var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V;
        double var_fast_sodium_current_m_gate__E0_m = var_fast_sodium_current_m_gate__V + 47.13;
        double var_fast_sodium_current_m_gate__alpha_m = (0.32 * var_fast_sodium_current_m_gate__E0_m) / (1.0 - exp((-0.1) * var_fast_sodium_current_m_gate__E0_m));
        double var_fast_sodium_current_m_gate__beta_m = 0.08 * exp((-var_fast_sodium_current_m_gate__V) / 11.0);
        double var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V;
        double var_fast_sodium_current_h_gate__alpha_h = 0.135 * exp((var_fast_sodium_current_h_gate__V + 80.0) / (-6.8));
        double var_fast_sodium_current_h_gate__beta_h = 7.5 / (1.0 + exp((-0.1) * (var_fast_sodium_current_h_gate__V + 11.0)));
        double var_fast_sodium_current_j_gate__V = var_fast_sodium_current__V;
        double var_fast_sodium_current_j_gate__alpha_j = (0.175 * exp((var_fast_sodium_current_j_gate__V + 100.0) / (-23.0))) / (1.0 + exp(0.15 * (var_fast_sodium_current_j_gate__V + 79.0)));
        double var_fast_sodium_current_j_gate__beta_j = 0.3 / (1.0 + exp((-0.1) * (var_fast_sodium_current_j_gate__V + 32.0)));
        double var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__V = var_rapid_activating_delayed_rectifiyer_K_current__V;
        double var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr_inf = 1.0 / (1.0 + exp((-2.182) - (0.1819 * var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__V)));
        double var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__tau_X_kr = 43.0 + (1.0 / (exp((-5.495) + (0.1691 * var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__V)) + exp((-7.677) - (0.0128 * var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__V))));
        double var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V = var_slow_activating_delayed_rectifiyer_K_current__V;
        double var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__tau_X_ks = 1.0 / (((7.19e-05 * (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 10.0)) / (1.0 - exp((-0.148) * (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 10.0)))) + ((0.000131 * (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 10.0)) / (exp(0.0687 * (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 10.0)) - 1.0)));
        double var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks_infinity = 1.0 / (1.0 + exp((var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 16.0) / (-13.6)));
        double var_transient_outward_potassium_current_X_to_gate__V = var_transient_outward_potassium_current__V;
        double var_transient_outward_potassium_current_X_to_gate__alpha_X_to = 0.04516 * exp(0.03577 * var_transient_outward_potassium_current_X_to_gate__V);
        double var_transient_outward_potassium_current_X_to_gate__beta_X_to = 0.0989 * exp((-0.06237) * var_transient_outward_potassium_current_X_to_gate__V);
        double var_transient_outward_potassium_current_Y_to_gate__V = var_transient_outward_potassium_current__V;
        double var_transient_outward_potassium_current_Y_to_gate__alpha_Y_to = (0.005415 * exp((var_transient_outward_potassium_current_Y_to_gate__V + 33.5) / (-5.0))) / (1.0 + (0.051335 * exp((var_transient_outward_potassium_current_Y_to_gate__V + 33.5) / (-5.0))));
        double var_transient_outward_potassium_current_Y_to_gate__beta_Y_to = (0.005415 * exp((var_transient_outward_potassium_current_Y_to_gate__V + 33.5) / 5.0)) / (1.0 + (0.051335 * exp((var_transient_outward_potassium_current_Y_to_gate__V + 33.5) / 5.0)));
        double var_L_type_Ca_current_f_gate__V = var_L_type_Ca_current__V;
        double var_L_type_Ca_current_f_gate__f_infinity = 1.0 / (1.0 + exp((var_L_type_Ca_current_f_gate__V + 12.5) / 5.0));
        double var_L_type_Ca_current_f_gate__tau_f = 30.0 + (200.0 / (1.0 + exp((var_L_type_Ca_current_f_gate__V + 20.0) / 9.5)));
        double var_L_type_Ca_current_d_gate__V = var_L_type_Ca_current__V;
        double var_L_type_Ca_current_d_gate__d_infinity = 1.0 / (1.0 + exp((var_L_type_Ca_current_d_gate__V + 10.0) / (-6.24)));
        double var_L_type_Ca_current_d_gate__E0_m = var_L_type_Ca_current_d_gate__V + 40.0;
        double var_L_type_Ca_current_d_gate__tau_d = 1.0 / (((0.25 * exp((-0.01) * var_L_type_Ca_current_d_gate__V)) / (1.0 + exp((-0.07) * var_L_type_Ca_current_d_gate__V))) + ((0.07 * exp((-0.05) * var_L_type_Ca_current_d_gate__E0_m)) / (1.0 + exp(0.05 * var_L_type_Ca_current_d_gate__E0_m))));

        const double _g_0 = var_L_type_Ca_current_d_gate__d_infinity / var_L_type_Ca_current_d_gate__tau_d;
        const double _h_0 = (-1.0) / var_L_type_Ca_current_d_gate__tau_d;
        const double _g_1 = var_L_type_Ca_current_f_gate__f_infinity / var_L_type_Ca_current_f_gate__tau_f;
        const double _h_1 = (-1.0) / var_L_type_Ca_current_f_gate__tau_f;
        const double _g_2 = var_fast_sodium_current_h_gate__alpha_h * 1.0;
        const double _h_2 = (var_fast_sodium_current_h_gate__alpha_h * (-1.0)) - (var_fast_sodium_current_h_gate__beta_h * 1.0);
        const double _g_3 = var_fast_sodium_current_j_gate__alpha_j * 1.0;
        const double _h_3 = (var_fast_sodium_current_j_gate__alpha_j * (-1.0)) - (var_fast_sodium_current_j_gate__beta_j * 1.0);
        const double _g_4 = var_fast_sodium_current_m_gate__alpha_m * 1.0;
        const double _h_4 = (var_fast_sodium_current_m_gate__alpha_m * (-1.0)) - (var_fast_sodium_current_m_gate__beta_m * 1.0);
        const double _g_5 = var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr_inf / var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__tau_X_kr;
        const double _h_5 = (-1.0) / var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__tau_X_kr;
        const double _g_6 = var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks_infinity / var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__tau_X_ks;
        const double _h_6 = (-1.0) / var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__tau_X_ks;
        const double _g_7 = var_transient_outward_potassium_current_X_to_gate__alpha_X_to * 1.0;
        const double _h_7 = (var_transient_outward_potassium_current_X_to_gate__alpha_X_to * (-1.0)) - (var_transient_outward_potassium_current_X_to_gate__beta_X_to * 1.0);
        const double _g_8 = var_transient_outward_potassium_current_Y_to_gate__alpha_Y_to * 1.0;
        const double _h_8 = (var_transient_outward_potassium_current_Y_to_gate__alpha_Y_to * (-1.0)) - (var_transient_outward_potassium_current_Y_to_gate__beta_Y_to * 1.0);

        rY[1] = (rY[1] + _g_4*mDt) / (1 - _h_4*mDt);
        rY[2] = (rY[2] + _g_2*mDt) / (1 - _h_2*mDt);
        rY[3] = (rY[3] + _g_3*mDt) / (1 - _h_3*mDt);
        rY[4] = (rY[4] + _g_5*mDt) / (1 - _h_5*mDt);
        rY[5] = (rY[5] + _g_6*mDt) / (1 - _h_6*mDt);
        rY[6] = (rY[6] + _g_7*mDt) / (1 - _h_7*mDt);
        rY[7] = (rY[7] + _g_8*mDt) / (1 - _h_8*mDt);
        rY[8] = (rY[8] + _g_1*mDt) / (1 - _h_1*mDt);
        rY[9] = (rY[9] + _g_0*mDt) / (1 - _h_0*mDt);

        double _guess[3] = {rY[10],rY[11],rY[12]};
        CardiacNewtonSolver<3> *_solver = CardiacNewtonSolver<3>::Instance();
        _solver->Solve(*this, _guess);
        rY[10] = _guess[0];
        rY[12] = _guess[1];
        rY[11] = _guess[2];
    }

};

#endif
