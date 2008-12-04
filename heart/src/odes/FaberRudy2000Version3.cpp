#ifndef _FaberRudy2000Version3_
#define _FaberRudy2000Version3_

// Model: LR_Dynamic_model_2000
// Processed by pycml - CellML Tools in Python
//     (translate: 1004, pycml: 896)
// on Thu Dec 13 17:50:34 2007

#include <cmath>
#include <cassert>
#include "AbstractCardiacCell.hpp"
#include "Exception.hpp"
#include "AbstractStimulusFunction.hpp"
#include "OdeSystemInformation.hpp"

class FaberRudy2000Version3 : public AbstractCardiacCell
{
private:
    // for heterogeneities
    double mScaleFactorGks;
    double mScaleFactorIto;
public:
    FaberRudy2000Version3(AbstractIvpOdeSolver *pSolver,
                              AbstractStimulusFunction *pIntracellularStimulus,
                              AbstractStimulusFunction *pExtracellularStimulus=NULL)
        : AbstractCardiacCell(pSolver, 25, 0, pIntracellularStimulus, pExtracellularStimulus)
    {
        mpSystemInfo = OdeSystemInformation<FaberRudy2000Version3>::Instance();

        mScaleFactorGks=1.0;
        mScaleFactorIto=0.0;
        
        Init();
    }

    ~FaberRudy2000Version3(void)
    {
    }

    void SetScaleFactorGks(double sfgks)
    {
        assert(sfgks>=0.0);
        mScaleFactorGks=sfgks;
    }

    void SetScaleFactorIto(double sfito)
    {
        assert(sfito>=0.0);
        mScaleFactorIto=sfito;
    }

    double GetIIonic()
    {
        std::vector<double>& rY = rGetStateVariables();
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -90
        double var_fast_sodium_current_m_gate__m = rY[1];
        // Units: dimensionless; Initial value: 0.0008
        double var_fast_sodium_current_h_gate__h = rY[2];
        // Units: dimensionless; Initial value: 0.993771
        double var_fast_sodium_current_j_gate__j = rY[3];
        // Units: dimensionless; Initial value: 0.995727
        double var_L_type_Ca_channel_d_gate__d = rY[4];
        // Units: dimensionless; Initial value: 3.210618e-6
        double var_L_type_Ca_channel_f_gate__f = rY[5];
        // Units: dimensionless; Initial value: 0.999837
        double var_T_type_Ca_channel_b_gate__b = rY[6];
        // Units: dimensionless; Initial value: 0.000970231
        double var_T_type_Ca_channel_g_gate__g = rY[7];
        // Units: dimensionless; Initial value: 0.994305
        double var_rapid_delayed_rectifier_potassium_current_xr_gate__xr = rY[8];
        // Units: dimensionless; Initial value: 0.000124042
        double var_slow_delayed_rectifier_potassium_current_xs1_gate__xs1 = rY[9];
        // Units: dimensionless; Initial value: 0.00445683
        double var_slow_delayed_rectifier_potassium_current_xs2_gate__xs2 = rY[10];
        // Units: dimensionless; Initial value: 0.00445683
        double var_transient_outward_current_zdv_gate__zdv = rY[11];
        // Units: dimensionless; Initial value: 0.5
        double var_transient_outward_current_ydv_gate__ydv = rY[12];
        // Units: dimensionless; Initial value: 0.5
        double var_calcium_dynamics__Cai = rY[13];
        // Units: millimolar; Initial value: 6e-5
        double var_ionic_concentrations__Nai = rY[23];
        // Units: millimolar; Initial value: 9
        double var_ionic_concentrations__Ki = rY[24];
        // Units: millimolar; Initial value: 141.2

        const double var_membrane__R = 8314.0;
        const double var_membrane__T = 310.0;
        const double var_membrane__F = 96485.0;
        double var_fast_sodium_current__j = var_fast_sodium_current_j_gate__j;
        double var_fast_sodium_current__h = var_fast_sodium_current_h_gate__h;
        const double var_fast_sodium_current__g_Na = 16.0;
        double var_fast_sodium_current__m = var_fast_sodium_current_m_gate__m;
        double var_fast_sodium_current__V = var_membrane__V;
        double var_fast_sodium_current__R = var_membrane__R;
        double var_fast_sodium_current__F = var_membrane__F;
        const double var_ionic_concentrations__Nao = 132.0;
        double var_fast_sodium_current__Nao = var_ionic_concentrations__Nao;
        double var_fast_sodium_current__Nai = var_ionic_concentrations__Nai;
        double var_fast_sodium_current__T = var_membrane__T;
        double var_fast_sodium_current__E_Na = ((var_fast_sodium_current__R * var_fast_sodium_current__T) / var_fast_sodium_current__F) * log(var_fast_sodium_current__Nao / var_fast_sodium_current__Nai);
        double var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na);
        double var_membrane__i_Na = var_fast_sodium_current__i_Na;
        double var_L_type_Ca_channel__d = var_L_type_Ca_channel_d_gate__d;
        double var_L_type_Ca_channel__f = var_L_type_Ca_channel_f_gate__f;
        const double var_L_type_Ca_channel_f_Ca_gate__Km_Ca = 0.0006;
        double var_L_type_Ca_channel__Cai = var_calcium_dynamics__Cai;
        double var_L_type_Ca_channel_f_Ca_gate__Cai = var_L_type_Ca_channel__Cai;
        double var_L_type_Ca_channel_f_Ca_gate__f_Ca = 1.0 / (1.0 + (var_L_type_Ca_channel_f_Ca_gate__Cai / var_L_type_Ca_channel_f_Ca_gate__Km_Ca));
        double var_L_type_Ca_channel__f_Ca = var_L_type_Ca_channel_f_Ca_gate__f_Ca;
        const double var_L_type_Ca_channel__gamma_Cai = 1.0;
        const double var_L_type_Ca_channel__gamma_Cao = 0.341;
        const double var_calcium_dynamics__Cao = 1.8;
        double var_L_type_Ca_channel__Cao = var_calcium_dynamics__Cao;
        double var_L_type_Ca_channel__F = var_membrane__F;
        const double var_L_type_Ca_channel__P_Ca = 0.00054;
        double var_L_type_Ca_channel__T = var_membrane__T;
        double var_L_type_Ca_channel__V = var_membrane__V;
        double var_L_type_Ca_channel__R = var_membrane__R;
        double var_L_type_Ca_channel__I_CaCa = (((var_L_type_Ca_channel__P_Ca * pow(2.0, 2.0) * var_L_type_Ca_channel__V * pow(var_L_type_Ca_channel__F, 2.0)) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) * ((var_L_type_Ca_channel__gamma_Cai * var_L_type_Ca_channel__Cai * exp((2.0 * var_L_type_Ca_channel__V * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__gamma_Cao * var_L_type_Ca_channel__Cao))) / (exp((2.0 * var_L_type_Ca_channel__V * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) - 1.0);
        double var_L_type_Ca_channel__i_CaCa = var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f_Ca * var_L_type_Ca_channel__I_CaCa;
        const double var_L_type_Ca_channel__P_Na = 6.75e-07;
        const double var_L_type_Ca_channel__gamma_Nao = 0.75;
        double var_L_type_Ca_channel__Nao = var_ionic_concentrations__Nao;
        double var_L_type_Ca_channel__Nai = var_ionic_concentrations__Nai;
        const double var_L_type_Ca_channel__gamma_Nai = 0.75;
        double var_L_type_Ca_channel__I_CaNa = (((var_L_type_Ca_channel__P_Na * pow(1.0, 2.0) * var_L_type_Ca_channel__V * pow(var_L_type_Ca_channel__F, 2.0)) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) * ((var_L_type_Ca_channel__gamma_Nai * var_L_type_Ca_channel__Nai * exp((1.0 * var_L_type_Ca_channel__V * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__gamma_Nao * var_L_type_Ca_channel__Nao))) / (exp((1.0 * var_L_type_Ca_channel__V * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) - 1.0);
        double var_L_type_Ca_channel__i_CaNa = var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f_Ca * var_L_type_Ca_channel__I_CaNa;
        const double var_ionic_concentrations__Ko = 4.5;
        double var_L_type_Ca_channel__Ko = var_ionic_concentrations__Ko;
        double var_L_type_Ca_channel__Ki = var_ionic_concentrations__Ki;
        const double var_L_type_Ca_channel__P_K = 1.93e-07;
        const double var_L_type_Ca_channel__gamma_Ki = 0.75;
        const double var_L_type_Ca_channel__gamma_Ko = 0.75;
        double var_L_type_Ca_channel__I_CaK = (((var_L_type_Ca_channel__P_K * pow(1.0, 2.0) * var_L_type_Ca_channel__V * pow(var_L_type_Ca_channel__F, 2.0)) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) * ((var_L_type_Ca_channel__gamma_Ki * var_L_type_Ca_channel__Ki * exp((1.0 * var_L_type_Ca_channel__V * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__gamma_Ko * var_L_type_Ca_channel__Ko))) / (exp((1.0 * var_L_type_Ca_channel__V * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) - 1.0);
        double var_L_type_Ca_channel__i_CaK = var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f_Ca * var_L_type_Ca_channel__I_CaK;
        double var_L_type_Ca_channel__i_Ca_L = var_L_type_Ca_channel__i_CaCa + var_L_type_Ca_channel__i_CaK + var_L_type_Ca_channel__i_CaNa;
        double var_membrane__i_Ca_L = var_L_type_Ca_channel__i_Ca_L;
        double var_T_type_Ca_channel__V = var_membrane__V;
        double var_T_type_Ca_channel__g = var_T_type_Ca_channel_g_gate__g;
        const double var_T_type_Ca_channel__g_CaT = 0.05;
        double var_T_type_Ca_channel__b = var_T_type_Ca_channel_b_gate__b;
        double var_calcium_background_current__R = var_membrane__R;
        double var_calcium_background_current__Cai = var_calcium_dynamics__Cai;
        double var_calcium_background_current__F = var_membrane__F;
        double var_calcium_background_current__T = var_membrane__T;
        double var_calcium_background_current__Cao = var_calcium_dynamics__Cao;
        double var_calcium_background_current__E_Ca = ((var_calcium_background_current__R * var_calcium_background_current__T) / (2.0 * var_calcium_background_current__F)) * log(var_calcium_background_current__Cao / var_calcium_background_current__Cai);
        double var_T_type_Ca_channel__E_Ca = var_calcium_background_current__E_Ca;
        double var_T_type_Ca_channel__i_Ca_T = var_T_type_Ca_channel__g_CaT * var_T_type_Ca_channel__b * var_T_type_Ca_channel__b * var_T_type_Ca_channel__g * (var_T_type_Ca_channel__V - var_T_type_Ca_channel__E_Ca);
        double var_membrane__i_Ca_T = var_T_type_Ca_channel__i_Ca_T;
        double var_rapid_delayed_rectifier_potassium_current__Ko = var_ionic_concentrations__Ko;
        double var_rapid_delayed_rectifier_potassium_current__g_Kr = 0.02614 * sqrt(var_rapid_delayed_rectifier_potassium_current__Ko / 5.4);
        double var_rapid_delayed_rectifier_potassium_current__V = var_membrane__V;
        double var_rapid_delayed_rectifier_potassium_current__xr = var_rapid_delayed_rectifier_potassium_current_xr_gate__xr;
        double var_rapid_delayed_rectifier_potassium_current__Rect = 1.0 / (1.0 + exp((var_rapid_delayed_rectifier_potassium_current__V + 9.0) / 22.4));
        double var_time_independent_potassium_current__Ki = var_ionic_concentrations__Ki;
        double var_time_independent_potassium_current__R = var_membrane__R;
        double var_time_independent_potassium_current__F = var_membrane__F;
        double var_time_independent_potassium_current__Ko = var_ionic_concentrations__Ko;
        double var_time_independent_potassium_current__T = var_membrane__T;
        double var_time_independent_potassium_current__E_K = ((var_time_independent_potassium_current__R * var_time_independent_potassium_current__T) / var_time_independent_potassium_current__F) * log(var_time_independent_potassium_current__Ko / var_time_independent_potassium_current__Ki);
        double var_rapid_delayed_rectifier_potassium_current__E_K = var_time_independent_potassium_current__E_K;
        double var_rapid_delayed_rectifier_potassium_current__i_Kr = var_rapid_delayed_rectifier_potassium_current__g_Kr * var_rapid_delayed_rectifier_potassium_current__xr * var_rapid_delayed_rectifier_potassium_current__Rect * (var_rapid_delayed_rectifier_potassium_current__V - var_rapid_delayed_rectifier_potassium_current__E_K);
        double var_membrane__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
        double var_slow_delayed_rectifier_potassium_current__xs2 = var_slow_delayed_rectifier_potassium_current_xs2_gate__xs2;
        double var_slow_delayed_rectifier_potassium_current__Nai = var_ionic_concentrations__Nai;
        double var_slow_delayed_rectifier_potassium_current__Nao = var_ionic_concentrations__Nao;
        double var_slow_delayed_rectifier_potassium_current__F = var_membrane__F;
        const double var_slow_delayed_rectifier_potassium_current__PNaK = 0.01833;
        double var_slow_delayed_rectifier_potassium_current__Ko = var_ionic_concentrations__Ko;
        double var_slow_delayed_rectifier_potassium_current__R = var_membrane__R;
        double var_slow_delayed_rectifier_potassium_current__T = var_membrane__T;
        double var_slow_delayed_rectifier_potassium_current__Ki = var_ionic_concentrations__Ki;
        double var_slow_delayed_rectifier_potassium_current__E_Ks = ((var_slow_delayed_rectifier_potassium_current__R * var_slow_delayed_rectifier_potassium_current__T) / var_slow_delayed_rectifier_potassium_current__F) * log((var_slow_delayed_rectifier_potassium_current__Ko + (var_slow_delayed_rectifier_potassium_current__PNaK * var_slow_delayed_rectifier_potassium_current__Nao)) / (var_slow_delayed_rectifier_potassium_current__Ki + (var_slow_delayed_rectifier_potassium_current__PNaK * var_slow_delayed_rectifier_potassium_current__Nai)));
        double var_slow_delayed_rectifier_potassium_current__xs1 = var_slow_delayed_rectifier_potassium_current_xs1_gate__xs1;
        double var_slow_delayed_rectifier_potassium_current__Cai = var_calcium_dynamics__Cai;
        double var_slow_delayed_rectifier_potassium_current__g_Ks = 0.433 * (1.0 + (0.6 / (1.0 + pow(3.8e-05 / var_slow_delayed_rectifier_potassium_current__Cai, 1.4))));
        double var_slow_delayed_rectifier_potassium_current__V = var_membrane__V;
        double var_slow_delayed_rectifier_potassium_current__i_Ks = mScaleFactorGks*var_slow_delayed_rectifier_potassium_current__g_Ks * var_slow_delayed_rectifier_potassium_current__xs1 * var_slow_delayed_rectifier_potassium_current__xs2 * (var_slow_delayed_rectifier_potassium_current__V - var_slow_delayed_rectifier_potassium_current__E_Ks);
        double var_membrane__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
        double var_sodium_activated_potassium_current__V = var_membrane__V;
        double var_sodium_activated_potassium_current__pov = 0.8 - (0.65 / (1.0 + exp((var_sodium_activated_potassium_current__V + 125.0) / 15.0)));
        double var_sodium_activated_potassium_current__g_K_Na = 0.0 * 0.12848;
        const double var_sodium_activated_potassium_current__nKNa = 2.8;
        double var_sodium_activated_potassium_current__Nai = var_ionic_concentrations__Nai;
        const double var_sodium_activated_potassium_current__kdKNa = 66.0;
        double var_sodium_activated_potassium_current__pona = 0.85 / (1.0 + pow(var_sodium_activated_potassium_current__kdKNa / var_sodium_activated_potassium_current__Nai, var_sodium_activated_potassium_current__nKNa));
        double var_sodium_activated_potassium_current__E_K = var_time_independent_potassium_current__E_K;
        double var_sodium_activated_potassium_current__i_K_Na = var_sodium_activated_potassium_current__g_K_Na * var_sodium_activated_potassium_current__pona * var_sodium_activated_potassium_current__pov * (var_sodium_activated_potassium_current__V - var_sodium_activated_potassium_current__E_K);
        double var_membrane__i_K_Na = var_sodium_activated_potassium_current__i_K_Na;
        double var_ATP_sensitive_potassium_current__V = var_membrane__V;
        const double var_ATP_sensitive_potassium_current__nATP = 0.24;
        const double var_ATP_sensitive_potassium_current__hATP = 2.0;
        const double var_ATP_sensitive_potassium_current__ATPi = 3.0;
        const double var_ATP_sensitive_potassium_current__kATP = 0.00025;
        double var_ATP_sensitive_potassium_current__pATP = 1.0 / (1.0 + pow(var_ATP_sensitive_potassium_current__ATPi / var_ATP_sensitive_potassium_current__kATP, var_ATP_sensitive_potassium_current__hATP));
        const double var_ATP_sensitive_potassium_current__i_K_ATP_on = 1.0;
        const double var_ATP_sensitive_potassium_current__nicholsarea = 5e-05;
        double var_ATP_sensitive_potassium_current__g_K_ATP = (var_ATP_sensitive_potassium_current__i_K_ATP_on * 0.000193) / var_ATP_sensitive_potassium_current__nicholsarea;
        double var_ATP_sensitive_potassium_current__Ko = var_ionic_concentrations__Ko;
        double var_ATP_sensitive_potassium_current__GKbaraATP = var_ATP_sensitive_potassium_current__g_K_ATP * var_ATP_sensitive_potassium_current__pATP * pow(var_ATP_sensitive_potassium_current__Ko / 4.0, var_ATP_sensitive_potassium_current__nATP);
        double var_ATP_sensitive_potassium_current__E_K = var_time_independent_potassium_current__E_K;
        double var_ATP_sensitive_potassium_current__i_K_ATP = var_ATP_sensitive_potassium_current__GKbaraATP * (var_ATP_sensitive_potassium_current__V - var_ATP_sensitive_potassium_current__E_K);
        double var_membrane__i_K_ATP = var_ATP_sensitive_potassium_current__i_K_ATP;
        double var_transient_outward_current__g_to = mScaleFactorIto* 0.5;
        double var_transient_outward_current__V = var_membrane__V;
        double var_transient_outward_current__E_K = var_time_independent_potassium_current__E_K;
        double var_transient_outward_current__zdv = var_transient_outward_current_zdv_gate__zdv;
        double var_transient_outward_current__ydv = var_transient_outward_current_ydv_gate__ydv;
        double var_transient_outward_current__rvdv = exp(var_transient_outward_current__V / 100.0);
        double var_transient_outward_current__i_to = var_transient_outward_current__g_to * pow(var_transient_outward_current__zdv, 3.0) * var_transient_outward_current__ydv * var_transient_outward_current__rvdv * (var_transient_outward_current__V - var_transient_outward_current__E_K);
        double var_membrane__i_to = var_transient_outward_current__i_to;
        double var_Na_Ca_exchanger__Nao = var_ionic_concentrations__Nao;
        double var_Na_Ca_exchanger__F = var_membrane__F;
        const double var_Na_Ca_exchanger__gamma = 0.15;
        double var_Na_Ca_exchanger__Nai = var_ionic_concentrations__Nai;
        const double var_Na_Ca_exchanger__c1 = 0.00025;
        const double var_Na_Ca_exchanger__c2 = 0.0001;
        double var_Na_Ca_exchanger__T = var_membrane__T;
        double var_Na_Ca_exchanger__Cao = var_calcium_dynamics__Cao;
        double var_Na_Ca_exchanger__V = var_membrane__V;
        double var_Na_Ca_exchanger__R = var_membrane__R;
        double var_Na_Ca_exchanger__Cai = var_calcium_dynamics__Cai;
        double var_Na_Ca_exchanger__i_NaCa = (var_Na_Ca_exchanger__c1 * exp(((var_Na_Ca_exchanger__gamma - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * ((exp((var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Nai, 3.0) * var_Na_Ca_exchanger__Cao) - (pow(var_Na_Ca_exchanger__Nao, 3.0) * var_Na_Ca_exchanger__Cai))) / (1.0 + (var_Na_Ca_exchanger__c2 * exp(((var_Na_Ca_exchanger__gamma - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * ((exp((var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Nai, 3.0) * var_Na_Ca_exchanger__Cao) + (pow(var_Na_Ca_exchanger__Nao, 3.0) * var_Na_Ca_exchanger__Cai))));
        double var_membrane__i_NaCa = var_Na_Ca_exchanger__i_NaCa;
        double var_time_independent_potassium_current__V = var_membrane__V;
        double var_time_independent_potassium_current_K1_gate__V = var_time_independent_potassium_current__V;
        double var_time_independent_potassium_current_K1_gate__E_K = var_time_independent_potassium_current__E_K;
        double var_time_independent_potassium_current_K1_gate__beta_K1 = (1000.0 * ((0.49124 * exp(0.08032 * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K) + 5.476))) + exp(0.06175 * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K) - 594.31)))) / (1.0 + exp((-0.5143) * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K) + 4.753)));
        double var_time_independent_potassium_current_K1_gate__alpha_K1 = 1020.0 / (1.0 + exp(0.2385 * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K) - 59.215)));
        double var_time_independent_potassium_current_K1_gate__K1_infinity = var_time_independent_potassium_current_K1_gate__alpha_K1 / (var_time_independent_potassium_current_K1_gate__alpha_K1 + var_time_independent_potassium_current_K1_gate__beta_K1);
        double var_time_independent_potassium_current__K1_infinity = var_time_independent_potassium_current_K1_gate__K1_infinity;
        double var_time_independent_potassium_current__g_K1 = 0.75 * sqrt(var_time_independent_potassium_current__Ko / 5.4);
        double var_time_independent_potassium_current__i_K1 = var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K1_infinity * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
        const double var_plateau_potassium_current__g_Kp = 0.00552;
        double var_plateau_potassium_current__V = var_membrane__V;
        double var_plateau_potassium_current__Kp = 1.0 / (1.0 + exp((7.488 - var_plateau_potassium_current__V) / 5.98));
        double var_plateau_potassium_current__E_K = var_time_independent_potassium_current__E_K;
        double var_plateau_potassium_current__i_Kp = var_plateau_potassium_current__g_Kp * var_plateau_potassium_current__Kp * (var_plateau_potassium_current__V - var_plateau_potassium_current__E_K);
        double var_membrane__i_Kp = var_plateau_potassium_current__i_Kp;
        const double var_sarcolemmal_calcium_pump__I_pCa = 1.15;
        double var_sarcolemmal_calcium_pump__Cai = var_calcium_dynamics__Cai;
        const double var_sarcolemmal_calcium_pump__K_mpCa = 0.0005;
        double var_sarcolemmal_calcium_pump__i_p_Ca = (var_sarcolemmal_calcium_pump__I_pCa * var_sarcolemmal_calcium_pump__Cai) / (var_sarcolemmal_calcium_pump__K_mpCa + var_sarcolemmal_calcium_pump__Cai);
        double var_membrane__i_p_Ca = var_sarcolemmal_calcium_pump__i_p_Ca;
        const double var_sodium_background_current__g_Nab = 0.004;
        double var_sodium_background_current__V = var_membrane__V;
        double var_sodium_background_current__E_Na = var_fast_sodium_current__E_Na;
        double var_sodium_background_current__i_Na_b = var_sodium_background_current__g_Nab * (var_sodium_background_current__V - var_sodium_background_current__E_Na);
        double var_membrane__i_Na_b = var_sodium_background_current__i_Na_b;
        const double var_calcium_background_current__g_Cab = 0.003016;
        double var_calcium_background_current__V = var_membrane__V;
        double var_calcium_background_current__i_Ca_b = var_calcium_background_current__g_Cab * (var_calcium_background_current__V - var_calcium_background_current__E_Ca);
        double var_membrane__i_Ca_b = var_calcium_background_current__i_Ca_b;
        double var_sodium_potassium_pump__Ko = var_ionic_concentrations__Ko;
        const double var_sodium_potassium_pump__K_mNai = 10.0;
        double var_sodium_potassium_pump__Nai = var_ionic_concentrations__Nai;
        double var_sodium_potassium_pump__V = var_membrane__V;
        double var_sodium_potassium_pump__F = var_membrane__F;
        double var_sodium_potassium_pump__T = var_membrane__T;
        double var_sodium_potassium_pump__Nao = var_ionic_concentrations__Nao;
        double var_sodium_potassium_pump__sigma = (1.0 / 7.0) * (exp(var_sodium_potassium_pump__Nao / 67.3) - 1.0);
        double var_sodium_potassium_pump__R = var_membrane__R;
        double var_sodium_potassium_pump__f_NaK = 1.0 / (1.0 + (0.1245 * exp(((-0.1) * var_sodium_potassium_pump__V * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))) + (0.0365 * var_sodium_potassium_pump__sigma * exp(((-var_sodium_potassium_pump__V) * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))));
        const double var_sodium_potassium_pump__I_NaK = 2.25;
        const double var_sodium_potassium_pump__K_mKo = 1.5;
        double var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__I_NaK * var_sodium_potassium_pump__f_NaK * 1.0) / (1.0 + pow(var_sodium_potassium_pump__K_mNai / var_sodium_potassium_pump__Nai, 2.0))) * var_sodium_potassium_pump__Ko) / (var_sodium_potassium_pump__Ko + var_sodium_potassium_pump__K_mKo);
        double var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
        double var_non_specific_calcium_activated_current__F = var_membrane__F;
        double var_non_specific_calcium_activated_current__R = var_membrane__R;
        double var_non_specific_calcium_activated_current__T = var_membrane__T;
        double var_non_specific_calcium_activated_current__V = var_membrane__V;
        double var_non_specific_calcium_activated_current__P_ns_Ca = 0.0 * 1.75e-07;
        double var_non_specific_calcium_activated_current__Nao = var_ionic_concentrations__Nao;
        double var_non_specific_calcium_activated_current__gamma_Nai = var_L_type_Ca_channel__gamma_Nai;
        double var_non_specific_calcium_activated_current__gamma_Nao = var_L_type_Ca_channel__gamma_Nao;
        double var_non_specific_calcium_activated_current__Nai = var_ionic_concentrations__Nai;
        double var_non_specific_calcium_activated_current__I_ns_Na = (((var_non_specific_calcium_activated_current__P_ns_Ca * pow(1.0, 2.0) * var_non_specific_calcium_activated_current__V * pow(var_non_specific_calcium_activated_current__F, 2.0)) / (var_non_specific_calcium_activated_current__R * var_non_specific_calcium_activated_current__T)) * ((var_non_specific_calcium_activated_current__gamma_Nai * var_non_specific_calcium_activated_current__Nai * exp((1.0 * var_non_specific_calcium_activated_current__V * var_non_specific_calcium_activated_current__F) / (var_non_specific_calcium_activated_current__R * var_non_specific_calcium_activated_current__T))) - (var_non_specific_calcium_activated_current__gamma_Nao * var_non_specific_calcium_activated_current__Nao))) / (exp((1.0 * var_non_specific_calcium_activated_current__V * var_non_specific_calcium_activated_current__F) / (var_non_specific_calcium_activated_current__R * var_non_specific_calcium_activated_current__T)) - 1.0);
        const double var_non_specific_calcium_activated_current__K_m_ns_Ca = 0.0012;
        double var_non_specific_calcium_activated_current__Cai = var_calcium_dynamics__Cai;
        double var_non_specific_calcium_activated_current__i_ns_Na = (var_non_specific_calcium_activated_current__I_ns_Na * 1.0) / (1.0 + pow(var_non_specific_calcium_activated_current__K_m_ns_Ca / var_non_specific_calcium_activated_current__Cai, 3.0));
        double var_non_specific_calcium_activated_current__Ko = var_ionic_concentrations__Ko;
        double var_non_specific_calcium_activated_current__Ki = var_ionic_concentrations__Ki;
        double var_non_specific_calcium_activated_current__gamma_Ko = var_L_type_Ca_channel__gamma_Ko;
        double var_non_specific_calcium_activated_current__gamma_Ki = var_L_type_Ca_channel__gamma_Ki;
        double var_non_specific_calcium_activated_current__I_ns_K = (((var_non_specific_calcium_activated_current__P_ns_Ca * pow(1.0, 2.0) * var_non_specific_calcium_activated_current__V * pow(var_non_specific_calcium_activated_current__F, 2.0)) / (var_non_specific_calcium_activated_current__R * var_non_specific_calcium_activated_current__T)) * ((var_non_specific_calcium_activated_current__gamma_Ki * var_non_specific_calcium_activated_current__Ki * exp((1.0 * var_non_specific_calcium_activated_current__V * var_non_specific_calcium_activated_current__F) / (var_non_specific_calcium_activated_current__R * var_non_specific_calcium_activated_current__T))) - (var_non_specific_calcium_activated_current__gamma_Ko * var_non_specific_calcium_activated_current__Ko))) / (exp((1.0 * var_non_specific_calcium_activated_current__V * var_non_specific_calcium_activated_current__F) / (var_non_specific_calcium_activated_current__R * var_non_specific_calcium_activated_current__T)) - 1.0);
        double var_non_specific_calcium_activated_current__i_ns_K = (var_non_specific_calcium_activated_current__I_ns_K * 1.0) / (1.0 + pow(var_non_specific_calcium_activated_current__K_m_ns_Ca / var_non_specific_calcium_activated_current__Cai, 3.0));
        double var_non_specific_calcium_activated_current__i_ns_Ca = var_non_specific_calcium_activated_current__i_ns_Na + var_non_specific_calcium_activated_current__i_ns_K;
        double var_membrane__i_ns_Ca = var_non_specific_calcium_activated_current__i_ns_Ca;

        return var_membrane__i_Na+var_membrane__i_Ca_L+var_membrane__i_Ca_T+var_membrane__i_Kr+var_membrane__i_Ks+var_membrane__i_K_Na+var_membrane__i_K_ATP+var_membrane__i_to+var_membrane__i_NaCa+var_membrane__i_K1+var_membrane__i_Kp+var_membrane__i_p_Ca+var_membrane__i_Na_b+var_membrane__i_Ca_b+var_membrane__i_NaK+var_membrane__i_ns_Ca;
    }

    void EvaluateYDerivatives (
            double var_environment__time,
            const std::vector<double> &rY,
            std::vector<double> &rDY)
    {
        // Inputs:
        // Time units: second, but input in ms
        var_environment__time *= 1e-3;
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -90
        double var_fast_sodium_current_m_gate__m = rY[1];
        // Units: dimensionless; Initial value: 0.0008
        double var_fast_sodium_current_h_gate__h = rY[2];
        // Units: dimensionless; Initial value: 0.993771
        double var_fast_sodium_current_j_gate__j = rY[3];
        // Units: dimensionless; Initial value: 0.995727
        double var_L_type_Ca_channel_d_gate__d = rY[4];
        // Units: dimensionless; Initial value: 3.210618e-6
        double var_L_type_Ca_channel_f_gate__f = rY[5];
        // Units: dimensionless; Initial value: 0.999837
        double var_T_type_Ca_channel_b_gate__b = rY[6];
        // Units: dimensionless; Initial value: 0.000970231
        double var_T_type_Ca_channel_g_gate__g = rY[7];
        // Units: dimensionless; Initial value: 0.994305
        double var_rapid_delayed_rectifier_potassium_current_xr_gate__xr = rY[8];
        // Units: dimensionless; Initial value: 0.000124042
        double var_slow_delayed_rectifier_potassium_current_xs1_gate__xs1 = rY[9];
        // Units: dimensionless; Initial value: 0.00445683
        double var_slow_delayed_rectifier_potassium_current_xs2_gate__xs2 = rY[10];
        // Units: dimensionless; Initial value: 0.00445683
        double var_transient_outward_current_zdv_gate__zdv = rY[11];
        // Units: dimensionless; Initial value: 0.5
        double var_transient_outward_current_ydv_gate__ydv = rY[12];
        // Units: dimensionless; Initial value: 0.5
        double var_calcium_dynamics__Cai = rY[13];
        // Units: millimolar; Initial value: 6e-5
        double var_calcium_dynamics__Ca_JSR = rY[14];
        // Units: millimolar; Initial value: 1.8
        double var_calcium_dynamics__Ca_NSR = rY[15];
        // Units: millimolar; Initial value: 1.8
        double var_calcium_dynamics__APtrack = rY[16];
        // Units: dimensionless; Initial value: 0
        double var_calcium_dynamics__APtrack2 = rY[17];
        // Units: dimensionless; Initial value: 0
        double var_calcium_dynamics__APtrack3 = rY[18];
        // Units: dimensionless; Initial value: 0
        double var_calcium_dynamics__Cainfluxtrack = rY[19];
        // Units: dimensionless; Initial value: 0
        double var_calcium_dynamics__OVRLDtrack = rY[20];
        // Units: dimensionless; Initial value: 0
        double var_calcium_dynamics__OVRLDtrack2 = rY[21];
        // Units: dimensionless; Initial value: 0
        double var_calcium_dynamics__OVRLDtrack3 = rY[22];
        // Units: dimensionless; Initial value: 0
        double var_ionic_concentrations__Nai = rY[23];
        // Units: millimolar; Initial value: 9
        double var_ionic_concentrations__Ki = rY[24];
        // Units: millimolar; Initial value: 141.2


        // Mathematics
        const double var_membrane__R = 8314.0;
        const double var_membrane__T = 310.0;
        const double var_membrane__F = 96485.0;
        const double var_membrane__Cm = 0.001;
//        const double var_membrane__stim_end = 100000.0;
//        const double var_membrane__stim_amplitude =  -25.5;
//        double var_membrane__time = var_environment__time;
//        const double var_membrane__stim_duration = 0.002;
//        const double var_membrane__stim_period = 1.0;
//        const double var_membrane__stim_start = 0.1;
//        double var_membrane__I_st = ((var_membrane__time >= var_membrane__stim_start) && (var_membrane__time <= var_membrane__stim_end) && (((var_membrane__time - var_membrane__stim_start) - (floor((var_membrane__time - var_membrane__stim_start) / var_membrane__stim_period) * var_membrane__stim_period)) <= var_membrane__stim_duration)) ? var_membrane__stim_amplitude : 0.0;
        double var_membrane__I_st = GetStimulus(var_environment__time*1000);
        double var_fast_sodium_current__j = var_fast_sodium_current_j_gate__j;
        double var_fast_sodium_current__h = var_fast_sodium_current_h_gate__h;
        const double var_fast_sodium_current__g_Na = 16.0;
        double var_fast_sodium_current__m = var_fast_sodium_current_m_gate__m;
        double var_fast_sodium_current__V = var_membrane__V;
        double var_fast_sodium_current__R = var_membrane__R;
        double var_fast_sodium_current__F = var_membrane__F;
        const double var_ionic_concentrations__Nao = 132.0;
        double var_fast_sodium_current__Nao = var_ionic_concentrations__Nao;
        double var_fast_sodium_current__Nai = var_ionic_concentrations__Nai;
        double var_fast_sodium_current__T = var_membrane__T;
        double var_fast_sodium_current__E_Na = ((var_fast_sodium_current__R * var_fast_sodium_current__T) / var_fast_sodium_current__F) * log(var_fast_sodium_current__Nao / var_fast_sodium_current__Nai);
        double var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na);
        double var_membrane__i_Na = var_fast_sodium_current__i_Na;
        double var_L_type_Ca_channel__d = var_L_type_Ca_channel_d_gate__d;
        double var_L_type_Ca_channel__f = var_L_type_Ca_channel_f_gate__f;
        const double var_L_type_Ca_channel_f_Ca_gate__Km_Ca = 0.0006;
        double var_L_type_Ca_channel__Cai = var_calcium_dynamics__Cai;
        double var_L_type_Ca_channel_f_Ca_gate__Cai = var_L_type_Ca_channel__Cai;
        double var_L_type_Ca_channel_f_Ca_gate__f_Ca = 1.0 / (1.0 + (var_L_type_Ca_channel_f_Ca_gate__Cai / var_L_type_Ca_channel_f_Ca_gate__Km_Ca));
        double var_L_type_Ca_channel__f_Ca = var_L_type_Ca_channel_f_Ca_gate__f_Ca;
        const double var_L_type_Ca_channel__gamma_Cai = 1.0;
        const double var_L_type_Ca_channel__gamma_Cao = 0.341;
        const double var_calcium_dynamics__Cao = 1.8;
        double var_L_type_Ca_channel__Cao = var_calcium_dynamics__Cao;
        double var_L_type_Ca_channel__F = var_membrane__F;
        const double var_L_type_Ca_channel__P_Ca = 0.00054;
        double var_L_type_Ca_channel__T = var_membrane__T;
        double var_L_type_Ca_channel__V = var_membrane__V;
        double var_L_type_Ca_channel__R = var_membrane__R;
        double var_L_type_Ca_channel__I_CaCa = (((var_L_type_Ca_channel__P_Ca * pow(2.0, 2.0) * var_L_type_Ca_channel__V * pow(var_L_type_Ca_channel__F, 2.0)) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) * ((var_L_type_Ca_channel__gamma_Cai * var_L_type_Ca_channel__Cai * exp((2.0 * var_L_type_Ca_channel__V * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__gamma_Cao * var_L_type_Ca_channel__Cao))) / (exp((2.0 * var_L_type_Ca_channel__V * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) - 1.0);
        double var_L_type_Ca_channel__i_CaCa = var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f_Ca * var_L_type_Ca_channel__I_CaCa;
        const double var_L_type_Ca_channel__P_Na = 6.75e-07;
        const double var_L_type_Ca_channel__gamma_Nao = 0.75;
        double var_L_type_Ca_channel__Nao = var_ionic_concentrations__Nao;
        double var_L_type_Ca_channel__Nai = var_ionic_concentrations__Nai;
        const double var_L_type_Ca_channel__gamma_Nai = 0.75;
        double var_L_type_Ca_channel__I_CaNa = (((var_L_type_Ca_channel__P_Na * pow(1.0, 2.0) * var_L_type_Ca_channel__V * pow(var_L_type_Ca_channel__F, 2.0)) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) * ((var_L_type_Ca_channel__gamma_Nai * var_L_type_Ca_channel__Nai * exp((1.0 * var_L_type_Ca_channel__V * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__gamma_Nao * var_L_type_Ca_channel__Nao))) / (exp((1.0 * var_L_type_Ca_channel__V * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) - 1.0);
        double var_L_type_Ca_channel__i_CaNa = var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f_Ca * var_L_type_Ca_channel__I_CaNa;
        const double var_ionic_concentrations__Ko = 4.5;
        double var_L_type_Ca_channel__Ko = var_ionic_concentrations__Ko;
        double var_L_type_Ca_channel__Ki = var_ionic_concentrations__Ki;
        const double var_L_type_Ca_channel__P_K = 1.93e-07;
        const double var_L_type_Ca_channel__gamma_Ki = 0.75;
        const double var_L_type_Ca_channel__gamma_Ko = 0.75;
        double var_L_type_Ca_channel__I_CaK = (((var_L_type_Ca_channel__P_K * pow(1.0, 2.0) * var_L_type_Ca_channel__V * pow(var_L_type_Ca_channel__F, 2.0)) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) * ((var_L_type_Ca_channel__gamma_Ki * var_L_type_Ca_channel__Ki * exp((1.0 * var_L_type_Ca_channel__V * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__gamma_Ko * var_L_type_Ca_channel__Ko))) / (exp((1.0 * var_L_type_Ca_channel__V * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) - 1.0);
        double var_L_type_Ca_channel__i_CaK = var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f_Ca * var_L_type_Ca_channel__I_CaK;
        double var_L_type_Ca_channel__i_Ca_L = var_L_type_Ca_channel__i_CaCa + var_L_type_Ca_channel__i_CaK + var_L_type_Ca_channel__i_CaNa;
        double var_membrane__i_Ca_L = var_L_type_Ca_channel__i_Ca_L;
        double var_T_type_Ca_channel__V = var_membrane__V;
        double var_T_type_Ca_channel__g = var_T_type_Ca_channel_g_gate__g;
        const double var_T_type_Ca_channel__g_CaT = 0.05;
        double var_T_type_Ca_channel__b = var_T_type_Ca_channel_b_gate__b;
        double var_calcium_background_current__R = var_membrane__R;
        double var_calcium_background_current__Cai = var_calcium_dynamics__Cai;
        double var_calcium_background_current__F = var_membrane__F;
        double var_calcium_background_current__T = var_membrane__T;
        double var_calcium_background_current__Cao = var_calcium_dynamics__Cao;
        double var_calcium_background_current__E_Ca = ((var_calcium_background_current__R * var_calcium_background_current__T) / (2.0 * var_calcium_background_current__F)) * log(var_calcium_background_current__Cao / var_calcium_background_current__Cai);
        double var_T_type_Ca_channel__E_Ca = var_calcium_background_current__E_Ca;
        double var_T_type_Ca_channel__i_Ca_T = var_T_type_Ca_channel__g_CaT * var_T_type_Ca_channel__b * var_T_type_Ca_channel__b * var_T_type_Ca_channel__g * (var_T_type_Ca_channel__V - var_T_type_Ca_channel__E_Ca);
        double var_membrane__i_Ca_T = var_T_type_Ca_channel__i_Ca_T;
        double var_rapid_delayed_rectifier_potassium_current__Ko = var_ionic_concentrations__Ko;
        double var_rapid_delayed_rectifier_potassium_current__g_Kr = 0.02614 * sqrt(var_rapid_delayed_rectifier_potassium_current__Ko / 5.4);
        double var_rapid_delayed_rectifier_potassium_current__V = var_membrane__V;
        double var_rapid_delayed_rectifier_potassium_current__xr = var_rapid_delayed_rectifier_potassium_current_xr_gate__xr;
        double var_rapid_delayed_rectifier_potassium_current__Rect = 1.0 / (1.0 + exp((var_rapid_delayed_rectifier_potassium_current__V + 9.0) / 22.4));
        double var_time_independent_potassium_current__Ki = var_ionic_concentrations__Ki;
        double var_time_independent_potassium_current__R = var_membrane__R;
        double var_time_independent_potassium_current__F = var_membrane__F;
        double var_time_independent_potassium_current__Ko = var_ionic_concentrations__Ko;
        double var_time_independent_potassium_current__T = var_membrane__T;
        double var_time_independent_potassium_current__E_K = ((var_time_independent_potassium_current__R * var_time_independent_potassium_current__T) / var_time_independent_potassium_current__F) * log(var_time_independent_potassium_current__Ko / var_time_independent_potassium_current__Ki);
        double var_rapid_delayed_rectifier_potassium_current__E_K = var_time_independent_potassium_current__E_K;
        double var_rapid_delayed_rectifier_potassium_current__i_Kr = var_rapid_delayed_rectifier_potassium_current__g_Kr * var_rapid_delayed_rectifier_potassium_current__xr * var_rapid_delayed_rectifier_potassium_current__Rect * (var_rapid_delayed_rectifier_potassium_current__V - var_rapid_delayed_rectifier_potassium_current__E_K);
        double var_membrane__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
        double var_slow_delayed_rectifier_potassium_current__xs2 = var_slow_delayed_rectifier_potassium_current_xs2_gate__xs2;
        double var_slow_delayed_rectifier_potassium_current__Nai = var_ionic_concentrations__Nai;
        double var_slow_delayed_rectifier_potassium_current__Nao = var_ionic_concentrations__Nao;
        double var_slow_delayed_rectifier_potassium_current__F = var_membrane__F;
        const double var_slow_delayed_rectifier_potassium_current__PNaK = 0.01833;
        double var_slow_delayed_rectifier_potassium_current__Ko = var_ionic_concentrations__Ko;
        double var_slow_delayed_rectifier_potassium_current__R = var_membrane__R;
        double var_slow_delayed_rectifier_potassium_current__T = var_membrane__T;
        double var_slow_delayed_rectifier_potassium_current__Ki = var_ionic_concentrations__Ki;
        double var_slow_delayed_rectifier_potassium_current__E_Ks = ((var_slow_delayed_rectifier_potassium_current__R * var_slow_delayed_rectifier_potassium_current__T) / var_slow_delayed_rectifier_potassium_current__F) * log((var_slow_delayed_rectifier_potassium_current__Ko + (var_slow_delayed_rectifier_potassium_current__PNaK * var_slow_delayed_rectifier_potassium_current__Nao)) / (var_slow_delayed_rectifier_potassium_current__Ki + (var_slow_delayed_rectifier_potassium_current__PNaK * var_slow_delayed_rectifier_potassium_current__Nai)));
        double var_slow_delayed_rectifier_potassium_current__xs1 = var_slow_delayed_rectifier_potassium_current_xs1_gate__xs1;
        double var_slow_delayed_rectifier_potassium_current__Cai = var_calcium_dynamics__Cai;
        double var_slow_delayed_rectifier_potassium_current__g_Ks = 0.433 * (1.0 + (0.6 / (1.0 + pow(3.8e-05 / var_slow_delayed_rectifier_potassium_current__Cai, 1.4))));
        double var_slow_delayed_rectifier_potassium_current__V = var_membrane__V;
        double var_slow_delayed_rectifier_potassium_current__i_Ks = mScaleFactorGks*var_slow_delayed_rectifier_potassium_current__g_Ks * var_slow_delayed_rectifier_potassium_current__xs1 * var_slow_delayed_rectifier_potassium_current__xs2 * (var_slow_delayed_rectifier_potassium_current__V - var_slow_delayed_rectifier_potassium_current__E_Ks);
        double var_membrane__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
        double var_sodium_activated_potassium_current__V = var_membrane__V;
        double var_sodium_activated_potassium_current__pov = 0.8 - (0.65 / (1.0 + exp((var_sodium_activated_potassium_current__V + 125.0) / 15.0)));
        double var_sodium_activated_potassium_current__g_K_Na = 0.0 * 0.12848;
        const double var_sodium_activated_potassium_current__nKNa = 2.8;
        double var_sodium_activated_potassium_current__Nai = var_ionic_concentrations__Nai;
        const double var_sodium_activated_potassium_current__kdKNa = 66.0;
        double var_sodium_activated_potassium_current__pona = 0.85 / (1.0 + pow(var_sodium_activated_potassium_current__kdKNa / var_sodium_activated_potassium_current__Nai, var_sodium_activated_potassium_current__nKNa));
        double var_sodium_activated_potassium_current__E_K = var_time_independent_potassium_current__E_K;
        double var_sodium_activated_potassium_current__i_K_Na = var_sodium_activated_potassium_current__g_K_Na * var_sodium_activated_potassium_current__pona * var_sodium_activated_potassium_current__pov * (var_sodium_activated_potassium_current__V - var_sodium_activated_potassium_current__E_K);
        double var_membrane__i_K_Na = var_sodium_activated_potassium_current__i_K_Na;
        double var_ATP_sensitive_potassium_current__V = var_membrane__V;
        const double var_ATP_sensitive_potassium_current__nATP = 0.24;
        const double var_ATP_sensitive_potassium_current__hATP = 2.0;
        const double var_ATP_sensitive_potassium_current__ATPi = 3.0;
        const double var_ATP_sensitive_potassium_current__kATP = 0.00025;
        double var_ATP_sensitive_potassium_current__pATP = 1.0 / (1.0 + pow(var_ATP_sensitive_potassium_current__ATPi / var_ATP_sensitive_potassium_current__kATP, var_ATP_sensitive_potassium_current__hATP));
        const double var_ATP_sensitive_potassium_current__i_K_ATP_on = 1.0;
        const double var_ATP_sensitive_potassium_current__nicholsarea = 5e-05;
        double var_ATP_sensitive_potassium_current__g_K_ATP = (var_ATP_sensitive_potassium_current__i_K_ATP_on * 0.000193) / var_ATP_sensitive_potassium_current__nicholsarea;
        double var_ATP_sensitive_potassium_current__Ko = var_ionic_concentrations__Ko;
        double var_ATP_sensitive_potassium_current__GKbaraATP = var_ATP_sensitive_potassium_current__g_K_ATP * var_ATP_sensitive_potassium_current__pATP * pow(var_ATP_sensitive_potassium_current__Ko / 4.0, var_ATP_sensitive_potassium_current__nATP);
        double var_ATP_sensitive_potassium_current__E_K = var_time_independent_potassium_current__E_K;
        double var_ATP_sensitive_potassium_current__i_K_ATP = var_ATP_sensitive_potassium_current__GKbaraATP * (var_ATP_sensitive_potassium_current__V - var_ATP_sensitive_potassium_current__E_K);
        double var_membrane__i_K_ATP = var_ATP_sensitive_potassium_current__i_K_ATP;
        double var_transient_outward_current__g_to = mScaleFactorIto* 0.5;
        double var_transient_outward_current__V = var_membrane__V;
        double var_transient_outward_current__E_K = var_time_independent_potassium_current__E_K;
        double var_transient_outward_current__zdv = var_transient_outward_current_zdv_gate__zdv;
        double var_transient_outward_current__ydv = var_transient_outward_current_ydv_gate__ydv;
        double var_transient_outward_current__rvdv = exp(var_transient_outward_current__V / 100.0);
        double var_transient_outward_current__i_to = var_transient_outward_current__g_to * pow(var_transient_outward_current__zdv, 3.0) * var_transient_outward_current__ydv * var_transient_outward_current__rvdv * (var_transient_outward_current__V - var_transient_outward_current__E_K);
        double var_membrane__i_to = var_transient_outward_current__i_to;
        double var_Na_Ca_exchanger__Nao = var_ionic_concentrations__Nao;
        double var_Na_Ca_exchanger__F = var_membrane__F;
        const double var_Na_Ca_exchanger__gamma = 0.15;
        double var_Na_Ca_exchanger__Nai = var_ionic_concentrations__Nai;
        const double var_Na_Ca_exchanger__c1 = 0.00025;
        const double var_Na_Ca_exchanger__c2 = 0.0001;
        double var_Na_Ca_exchanger__T = var_membrane__T;
        double var_Na_Ca_exchanger__Cao = var_calcium_dynamics__Cao;
        double var_Na_Ca_exchanger__V = var_membrane__V;
        double var_Na_Ca_exchanger__R = var_membrane__R;
        double var_Na_Ca_exchanger__Cai = var_calcium_dynamics__Cai;
        double var_Na_Ca_exchanger__i_NaCa = (var_Na_Ca_exchanger__c1 * exp(((var_Na_Ca_exchanger__gamma - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * ((exp((var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Nai, 3.0) * var_Na_Ca_exchanger__Cao) - (pow(var_Na_Ca_exchanger__Nao, 3.0) * var_Na_Ca_exchanger__Cai))) / (1.0 + (var_Na_Ca_exchanger__c2 * exp(((var_Na_Ca_exchanger__gamma - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * ((exp((var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Nai, 3.0) * var_Na_Ca_exchanger__Cao) + (pow(var_Na_Ca_exchanger__Nao, 3.0) * var_Na_Ca_exchanger__Cai))));
        double var_membrane__i_NaCa = var_Na_Ca_exchanger__i_NaCa;
        double var_time_independent_potassium_current__V = var_membrane__V;
        double var_time_independent_potassium_current_K1_gate__V = var_time_independent_potassium_current__V;
        double var_time_independent_potassium_current_K1_gate__E_K = var_time_independent_potassium_current__E_K;
        double var_time_independent_potassium_current_K1_gate__beta_K1 = (1000.0 * ((0.49124 * exp(0.08032 * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K) + 5.476))) + exp(0.06175 * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K) - 594.31)))) / (1.0 + exp((-0.5143) * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K) + 4.753)));
        double var_time_independent_potassium_current_K1_gate__alpha_K1 = 1020.0 / (1.0 + exp(0.2385 * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K) - 59.215)));
        double var_time_independent_potassium_current_K1_gate__K1_infinity = var_time_independent_potassium_current_K1_gate__alpha_K1 / (var_time_independent_potassium_current_K1_gate__alpha_K1 + var_time_independent_potassium_current_K1_gate__beta_K1);
        double var_time_independent_potassium_current__K1_infinity = var_time_independent_potassium_current_K1_gate__K1_infinity;
        double var_time_independent_potassium_current__g_K1 = 0.75 * sqrt(var_time_independent_potassium_current__Ko / 5.4);
        double var_time_independent_potassium_current__i_K1 = var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K1_infinity * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
        const double var_plateau_potassium_current__g_Kp = 0.00552;
        double var_plateau_potassium_current__V = var_membrane__V;
        double var_plateau_potassium_current__Kp = 1.0 / (1.0 + exp((7.488 - var_plateau_potassium_current__V) / 5.98));
        double var_plateau_potassium_current__E_K = var_time_independent_potassium_current__E_K;
        double var_plateau_potassium_current__i_Kp = var_plateau_potassium_current__g_Kp * var_plateau_potassium_current__Kp * (var_plateau_potassium_current__V - var_plateau_potassium_current__E_K);
        double var_membrane__i_Kp = var_plateau_potassium_current__i_Kp;
        const double var_sarcolemmal_calcium_pump__I_pCa = 1.15;
        double var_sarcolemmal_calcium_pump__Cai = var_calcium_dynamics__Cai;
        const double var_sarcolemmal_calcium_pump__K_mpCa = 0.0005;
        double var_sarcolemmal_calcium_pump__i_p_Ca = (var_sarcolemmal_calcium_pump__I_pCa * var_sarcolemmal_calcium_pump__Cai) / (var_sarcolemmal_calcium_pump__K_mpCa + var_sarcolemmal_calcium_pump__Cai);
        double var_membrane__i_p_Ca = var_sarcolemmal_calcium_pump__i_p_Ca;
        const double var_sodium_background_current__g_Nab = 0.004;
        double var_sodium_background_current__V = var_membrane__V;
        double var_sodium_background_current__E_Na = var_fast_sodium_current__E_Na;
        double var_sodium_background_current__i_Na_b = var_sodium_background_current__g_Nab * (var_sodium_background_current__V - var_sodium_background_current__E_Na);
        double var_membrane__i_Na_b = var_sodium_background_current__i_Na_b;
        const double var_calcium_background_current__g_Cab = 0.003016;
        double var_calcium_background_current__V = var_membrane__V;
        double var_calcium_background_current__i_Ca_b = var_calcium_background_current__g_Cab * (var_calcium_background_current__V - var_calcium_background_current__E_Ca);
        double var_membrane__i_Ca_b = var_calcium_background_current__i_Ca_b;
        double var_sodium_potassium_pump__Ko = var_ionic_concentrations__Ko;
        const double var_sodium_potassium_pump__K_mNai = 10.0;
        double var_sodium_potassium_pump__Nai = var_ionic_concentrations__Nai;
        double var_sodium_potassium_pump__V = var_membrane__V;
        double var_sodium_potassium_pump__F = var_membrane__F;
        double var_sodium_potassium_pump__T = var_membrane__T;
        double var_sodium_potassium_pump__Nao = var_ionic_concentrations__Nao;
        double var_sodium_potassium_pump__sigma = (1.0 / 7.0) * (exp(var_sodium_potassium_pump__Nao / 67.3) - 1.0);
        double var_sodium_potassium_pump__R = var_membrane__R;
        double var_sodium_potassium_pump__f_NaK = 1.0 / (1.0 + (0.1245 * exp(((-0.1) * var_sodium_potassium_pump__V * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))) + (0.0365 * var_sodium_potassium_pump__sigma * exp(((-var_sodium_potassium_pump__V) * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))));
        const double var_sodium_potassium_pump__I_NaK = 2.25;
        const double var_sodium_potassium_pump__K_mKo = 1.5;
        double var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__I_NaK * var_sodium_potassium_pump__f_NaK * 1.0) / (1.0 + pow(var_sodium_potassium_pump__K_mNai / var_sodium_potassium_pump__Nai, 2.0))) * var_sodium_potassium_pump__Ko) / (var_sodium_potassium_pump__Ko + var_sodium_potassium_pump__K_mKo);
        double var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
        double var_non_specific_calcium_activated_current__F = var_membrane__F;
        double var_non_specific_calcium_activated_current__R = var_membrane__R;
        double var_non_specific_calcium_activated_current__T = var_membrane__T;
        double var_non_specific_calcium_activated_current__V = var_membrane__V;
        double var_non_specific_calcium_activated_current__P_ns_Ca = 0.0 * 1.75e-07;
        double var_non_specific_calcium_activated_current__Nao = var_ionic_concentrations__Nao;
        double var_non_specific_calcium_activated_current__gamma_Nai = var_L_type_Ca_channel__gamma_Nai;
        double var_non_specific_calcium_activated_current__gamma_Nao = var_L_type_Ca_channel__gamma_Nao;
        double var_non_specific_calcium_activated_current__Nai = var_ionic_concentrations__Nai;
        double var_non_specific_calcium_activated_current__I_ns_Na = (((var_non_specific_calcium_activated_current__P_ns_Ca * pow(1.0, 2.0) * var_non_specific_calcium_activated_current__V * pow(var_non_specific_calcium_activated_current__F, 2.0)) / (var_non_specific_calcium_activated_current__R * var_non_specific_calcium_activated_current__T)) * ((var_non_specific_calcium_activated_current__gamma_Nai * var_non_specific_calcium_activated_current__Nai * exp((1.0 * var_non_specific_calcium_activated_current__V * var_non_specific_calcium_activated_current__F) / (var_non_specific_calcium_activated_current__R * var_non_specific_calcium_activated_current__T))) - (var_non_specific_calcium_activated_current__gamma_Nao * var_non_specific_calcium_activated_current__Nao))) / (exp((1.0 * var_non_specific_calcium_activated_current__V * var_non_specific_calcium_activated_current__F) / (var_non_specific_calcium_activated_current__R * var_non_specific_calcium_activated_current__T)) - 1.0);
        const double var_non_specific_calcium_activated_current__K_m_ns_Ca = 0.0012;
        double var_non_specific_calcium_activated_current__Cai = var_calcium_dynamics__Cai;
        double var_non_specific_calcium_activated_current__i_ns_Na = (var_non_specific_calcium_activated_current__I_ns_Na * 1.0) / (1.0 + pow(var_non_specific_calcium_activated_current__K_m_ns_Ca / var_non_specific_calcium_activated_current__Cai, 3.0));
        double var_non_specific_calcium_activated_current__Ko = var_ionic_concentrations__Ko;
        double var_non_specific_calcium_activated_current__Ki = var_ionic_concentrations__Ki;
        double var_non_specific_calcium_activated_current__gamma_Ko = var_L_type_Ca_channel__gamma_Ko;
        double var_non_specific_calcium_activated_current__gamma_Ki = var_L_type_Ca_channel__gamma_Ki;
        double var_non_specific_calcium_activated_current__I_ns_K = (((var_non_specific_calcium_activated_current__P_ns_Ca * pow(1.0, 2.0) * var_non_specific_calcium_activated_current__V * pow(var_non_specific_calcium_activated_current__F, 2.0)) / (var_non_specific_calcium_activated_current__R * var_non_specific_calcium_activated_current__T)) * ((var_non_specific_calcium_activated_current__gamma_Ki * var_non_specific_calcium_activated_current__Ki * exp((1.0 * var_non_specific_calcium_activated_current__V * var_non_specific_calcium_activated_current__F) / (var_non_specific_calcium_activated_current__R * var_non_specific_calcium_activated_current__T))) - (var_non_specific_calcium_activated_current__gamma_Ko * var_non_specific_calcium_activated_current__Ko))) / (exp((1.0 * var_non_specific_calcium_activated_current__V * var_non_specific_calcium_activated_current__F) / (var_non_specific_calcium_activated_current__R * var_non_specific_calcium_activated_current__T)) - 1.0);
        double var_non_specific_calcium_activated_current__i_ns_K = (var_non_specific_calcium_activated_current__I_ns_K * 1.0) / (1.0 + pow(var_non_specific_calcium_activated_current__K_m_ns_Ca / var_non_specific_calcium_activated_current__Cai, 3.0));
        double var_non_specific_calcium_activated_current__i_ns_Ca = var_non_specific_calcium_activated_current__i_ns_Na + var_non_specific_calcium_activated_current__i_ns_K;
        double var_membrane__i_ns_Ca = var_non_specific_calcium_activated_current__i_ns_Ca;
        double var_membrane__dVdt = ((-1.0) / var_membrane__Cm) * (var_membrane__i_Na + var_membrane__i_Ca_L + var_membrane__i_Ca_T + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_K_Na + var_membrane__i_K_ATP + var_membrane__i_to + var_membrane__i_K1 + var_membrane__i_Kp + var_membrane__i_NaCa + var_membrane__i_p_Ca + var_membrane__i_Na_b + var_membrane__i_Ca_b + var_membrane__i_NaK + var_membrane__i_ns_Ca + var_membrane__I_st);
        double var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V;
        double var_fast_sodium_current_m_gate__E0_m = var_fast_sodium_current_m_gate__V + 47.13;
        const double var_fast_sodium_current_m_gate__delta_m = 1e-05;
        double var_fast_sodium_current_m_gate__alpha_m = (fabs(var_fast_sodium_current_m_gate__E0_m) >= var_fast_sodium_current_m_gate__delta_m) ? ((320.0 * var_fast_sodium_current_m_gate__E0_m) / (1.0 - exp((-0.1) * var_fast_sodium_current_m_gate__E0_m))) : 3200.0;
        double var_fast_sodium_current_m_gate__beta_m = 80.0 * exp((-var_fast_sodium_current_m_gate__V) / 11.0);
        double var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V;
        double var_fast_sodium_current_h_gate__alpha_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? (135.0 * exp((80.0 + var_fast_sodium_current_h_gate__V) / (-6.8))) : 0.0;
        double var_fast_sodium_current_h_gate__beta_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? ((3560.0 * exp(0.079 * var_fast_sodium_current_h_gate__V)) + (310000000.0 * exp(0.35 * var_fast_sodium_current_h_gate__V))) : (1000.0 / (0.13 * (1.0 + exp((var_fast_sodium_current_h_gate__V + 10.66) / (-11.1)))));
        double var_fast_sodium_current_j_gate__V = var_fast_sodium_current__V;
        double var_fast_sodium_current_j_gate__alpha_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((1000.0 * (-((127140.0 * exp(0.2444 * var_fast_sodium_current_j_gate__V)) + (3.474e-05 * exp((-0.04391) * var_fast_sodium_current_j_gate__V)))) * (var_fast_sodium_current_j_gate__V + 37.78)) / (1.0 + exp(0.311 * (var_fast_sodium_current_j_gate__V + 79.23)))) : 0.0;
        double var_fast_sodium_current_j_gate__beta_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((121.2 * exp((-0.01052) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1378) * (var_fast_sodium_current_j_gate__V + 40.14)))) : ((300.0 * exp((-2.535e-07) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1) * (var_fast_sodium_current_j_gate__V + 32.0))));
        double var_L_type_Ca_channel_d_gate__V = var_L_type_Ca_channel__V;
        double var_L_type_Ca_channel_d_gate__E0_d = var_L_type_Ca_channel_d_gate__V + 10.0;
        double var_L_type_Ca_channel_d_gate__d_infinity = 1.0 / (1.0 + exp((-var_L_type_Ca_channel_d_gate__E0_d) / 6.24));
        double var_L_type_Ca_channel_d_gate__tau_d = (fabs(var_L_type_Ca_channel_d_gate__E0_d) < 1e-05) ? (0.001 / (0.035 * 6.24)) : ((0.001 * var_L_type_Ca_channel_d_gate__d_infinity * (1.0 - exp((-var_L_type_Ca_channel_d_gate__E0_d) / 6.24))) / (0.035 * var_L_type_Ca_channel_d_gate__E0_d));
        double var_L_type_Ca_channel_d_gate__alpha_d = var_L_type_Ca_channel_d_gate__d_infinity / var_L_type_Ca_channel_d_gate__tau_d;
        double var_L_type_Ca_channel_d_gate__beta_d = (1.0 - var_L_type_Ca_channel_d_gate__d_infinity) / var_L_type_Ca_channel_d_gate__tau_d;
        double var_L_type_Ca_channel_f_gate__V = var_L_type_Ca_channel__V;
        double var_L_type_Ca_channel_f_gate__f_infinity = (1.0 / (1.0 + exp((var_L_type_Ca_channel_f_gate__V + 32.0) / 8.0))) + (0.6 / (1.0 + exp((50.0 - var_L_type_Ca_channel_f_gate__V) / 20.0)));
        double var_L_type_Ca_channel_f_gate__tau_f = 0.001 / ((0.0197 * exp(-pow(0.0337 * (var_L_type_Ca_channel_f_gate__V + 10.0), 2.0))) + 0.02);
        double var_L_type_Ca_channel_f_gate__alpha_f = var_L_type_Ca_channel_f_gate__f_infinity / var_L_type_Ca_channel_f_gate__tau_f;
        double var_L_type_Ca_channel_f_gate__beta_f = (1.0 - var_L_type_Ca_channel_f_gate__f_infinity) / var_L_type_Ca_channel_f_gate__tau_f;
        double var_T_type_Ca_channel_b_gate__V = var_T_type_Ca_channel__V;
        double var_T_type_Ca_channel_b_gate__b_inf = 1.0 / (1.0 + exp((-(var_T_type_Ca_channel_b_gate__V + 14.0)) / 10.8));
        double var_T_type_Ca_channel_b_gate__tau_b = 0.0037 + (0.0061 / (1.0 + exp((var_T_type_Ca_channel_b_gate__V + 25.0) / 4.5)));
        double var_T_type_Ca_channel_g_gate__V = var_T_type_Ca_channel__V;
        double var_T_type_Ca_channel_g_gate__g_inf = 1.0 / (1.0 + exp((var_T_type_Ca_channel_g_gate__V + 60.0) / 5.6));
        double var_T_type_Ca_channel_g_gate__tau_g = (var_T_type_Ca_channel_g_gate__V <= 0.0) ? (((-0.000875) * var_T_type_Ca_channel_g_gate__V) + 0.012) : 0.012;
        double var_rapid_delayed_rectifier_potassium_current_xr_gate__V = var_rapid_delayed_rectifier_potassium_current__V;
        double var_rapid_delayed_rectifier_potassium_current_xr_gate__xr_infinity = 1.0 / (1.0 + exp((-(var_rapid_delayed_rectifier_potassium_current_xr_gate__V + 21.5)) / 7.5));
        double var_rapid_delayed_rectifier_potassium_current_xr_gate__tau_xr = 0.001 / (((0.00138 * (var_rapid_delayed_rectifier_potassium_current_xr_gate__V + 14.2)) / (1.0 - exp((-0.123) * (var_rapid_delayed_rectifier_potassium_current_xr_gate__V + 14.2)))) + ((0.00061 * (var_rapid_delayed_rectifier_potassium_current_xr_gate__V + 38.9)) / (exp(0.145 * (var_rapid_delayed_rectifier_potassium_current_xr_gate__V + 38.9)) - 1.0)));
        double var_slow_delayed_rectifier_potassium_current_xs1_gate__V = var_slow_delayed_rectifier_potassium_current__V;
        double var_slow_delayed_rectifier_potassium_current_xs1_gate__xs1_infinity = 1.0 / (1.0 + exp((-(var_slow_delayed_rectifier_potassium_current_xs1_gate__V - 1.5)) / 16.7));
        double var_slow_delayed_rectifier_potassium_current_xs1_gate__tau_xs1 = 0.001 / (((7.19e-05 * (var_slow_delayed_rectifier_potassium_current_xs1_gate__V + 30.0)) / (1.0 - exp((-0.148) * (var_slow_delayed_rectifier_potassium_current_xs1_gate__V + 30.0)))) + ((0.000131 * (var_slow_delayed_rectifier_potassium_current_xs1_gate__V + 30.0)) / (exp(0.0687 * (var_slow_delayed_rectifier_potassium_current_xs1_gate__V + 30.0)) - 1.0)));
        double var_slow_delayed_rectifier_potassium_current_xs2_gate__V = var_slow_delayed_rectifier_potassium_current__V;
        double var_slow_delayed_rectifier_potassium_current_xs2_gate__xs2_infinity = 1.0 / (1.0 + exp((-(var_slow_delayed_rectifier_potassium_current_xs2_gate__V - 1.5)) / 16.7));
        double var_slow_delayed_rectifier_potassium_current_xs2_gate__tau_xs2 = (4.0 * 0.001) / (((7.19e-05 * (var_slow_delayed_rectifier_potassium_current_xs2_gate__V + 30.0)) / (1.0 - exp((-0.148) * (var_slow_delayed_rectifier_potassium_current_xs2_gate__V + 30.0)))) + ((0.000131 * (var_slow_delayed_rectifier_potassium_current_xs2_gate__V + 30.0)) / (exp(0.0687 * (var_slow_delayed_rectifier_potassium_current_xs2_gate__V + 30.0)) - 1.0)));
        double var_transient_outward_current_zdv_gate__V = var_transient_outward_current__V;
        double var_transient_outward_current_zdv_gate__alpha_zdv = (10000.0 * exp((var_transient_outward_current_zdv_gate__V - 40.0) / 25.0)) / (1.0 + exp((var_transient_outward_current_zdv_gate__V - 40.0) / 25.0));
        double var_transient_outward_current_zdv_gate__beta_zdv = (10000.0 * exp((-(var_transient_outward_current_zdv_gate__V + 90.0)) / 25.0)) / (1.0 + exp((-(var_transient_outward_current_zdv_gate__V + 90.0)) / 25.0));
        double var_transient_outward_current_zdv_gate__tau_zdv = 1.0 / (var_transient_outward_current_zdv_gate__alpha_zdv + var_transient_outward_current_zdv_gate__beta_zdv);
        double var_transient_outward_current_zdv_gate__zdv_ss = var_transient_outward_current_zdv_gate__alpha_zdv / (var_transient_outward_current_zdv_gate__alpha_zdv + var_transient_outward_current_zdv_gate__beta_zdv);
        double var_transient_outward_current_ydv_gate__V = var_transient_outward_current__V;
        double var_transient_outward_current_ydv_gate__alpha_ydv = 15.0 / (1.0 + exp((var_transient_outward_current_ydv_gate__V + 60.0) / 5.0));
        double var_transient_outward_current_ydv_gate__beta_ydv = (100.0 * exp((var_transient_outward_current_ydv_gate__V + 25.0) / 5.0)) / (1.0 + exp((var_transient_outward_current_ydv_gate__V + 25.0) / 5.0));
        double var_transient_outward_current_ydv_gate__tau_ydv = 1.0 / (var_transient_outward_current_ydv_gate__alpha_ydv + var_transient_outward_current_ydv_gate__beta_ydv);
        double var_transient_outward_current_ydv_gate__ydv_ss = var_transient_outward_current_ydv_gate__alpha_ydv / (var_transient_outward_current_ydv_gate__alpha_ydv + var_transient_outward_current_ydv_gate__beta_ydv);
        const double var_calcium_dynamics__delta_Ca_ith = 0.00018;
        const double var_calcium_dynamics__K_mrel = 0.0008;
        const double var_calcium_dynamics__G_rel_max = 60000.0;
        const double var_calcium_dynamics__G_rel_overload = 4000.0;
        double var_calcium_dynamics__G_rel = (var_calcium_dynamics__Cainfluxtrack > var_calcium_dynamics__delta_Ca_ith) ? (((var_calcium_dynamics__G_rel_max * (var_calcium_dynamics__Cainfluxtrack - var_calcium_dynamics__delta_Ca_ith)) / ((var_calcium_dynamics__K_mrel + var_calcium_dynamics__Cainfluxtrack) - var_calcium_dynamics__delta_Ca_ith)) * (1.0 - var_calcium_dynamics__APtrack2) * var_calcium_dynamics__APtrack2) : ((var_calcium_dynamics__Cainfluxtrack <= var_calcium_dynamics__delta_Ca_ith) && (var_calcium_dynamics__OVRLDtrack2 > 0.0)) ? (var_calcium_dynamics__G_rel_overload * (1.0 - var_calcium_dynamics__OVRLDtrack2) * var_calcium_dynamics__OVRLDtrack2) : 0.0;
        double var_calcium_dynamics__i_rel = var_calcium_dynamics__G_rel * (var_calcium_dynamics__Ca_JSR - var_calcium_dynamics__Cai);
        const double var_calcium_dynamics__K_mup = 0.00092;
        const double var_calcium_dynamics__I_up = 8.75;
        double var_calcium_dynamics__i_up = (var_calcium_dynamics__I_up * var_calcium_dynamics__Cai) / (var_calcium_dynamics__Cai + var_calcium_dynamics__K_mup);
        const double var_calcium_dynamics__Ca_NSR_max = 15.0;
        double var_calcium_dynamics__K_leak = var_calcium_dynamics__I_up / var_calcium_dynamics__Ca_NSR_max;
        double var_calcium_dynamics__i_leak = var_calcium_dynamics__K_leak * var_calcium_dynamics__Ca_NSR;
        const double var_calcium_dynamics__tau_tr = 0.18;
        double var_calcium_dynamics__i_tr = (var_calcium_dynamics__Ca_NSR - var_calcium_dynamics__Ca_JSR) / var_calcium_dynamics__tau_tr;
        const double var_calcium_dynamics__CSQN_max = 10.0;
        const double var_calcium_dynamics__K_mCSQN = 0.8;
        double var_calcium_dynamics__F = var_membrane__F;
        const double var_ionic_concentrations__preplength = 0.001;
        const double var_ionic_concentrations__radius = 0.00011;
        double var_ionic_concentrations__volume = M_PI * var_ionic_concentrations__preplength * pow(var_ionic_concentrations__radius, 2.0);
        double var_ionic_concentrations__V_myo = 0.68 * var_ionic_concentrations__volume;
        double var_calcium_dynamics__V_myo = var_ionic_concentrations__V_myo;
        const double var_ionic_concentrations__A_cap = 1.434e-07;
        double var_calcium_dynamics__A_cap = var_ionic_concentrations__A_cap;
        double var_calcium_dynamics__V_JSR = (0.0048 / 0.68) * var_calcium_dynamics__V_myo;
        double var_calcium_dynamics__V_NSR = (0.0552 / 0.68) * var_calcium_dynamics__V_myo;
        double var_calcium_dynamics__i_NaCa = var_Na_Ca_exchanger__i_NaCa;
        double var_calcium_dynamics__i_CaCa = var_L_type_Ca_channel__i_CaCa;
        double var_calcium_dynamics__i_p_Ca = var_sarcolemmal_calcium_pump__i_p_Ca;
        double var_calcium_dynamics__i_Ca_b = var_calcium_background_current__i_Ca_b;
        double var_calcium_dynamics__i_Ca_T = var_T_type_Ca_channel__i_Ca_T;
        const double var_calcium_dynamics__K_mTn = 0.0005;
        const double var_calcium_dynamics__K_mCMDN = 0.00238;
        const double var_calcium_dynamics__Tn_max = 0.07;
        const double var_calcium_dynamics__CMDN_max = 0.05;
        double var_calcium_dynamics__dVdt = var_membrane__dVdt;
        const double var_calcium_dynamics__CSQNthresh = 0.7;
        const double var_calcium_dynamics__Logicthresh = 0.98;
        double var_ionic_concentrations__F = var_membrane__F;
        double var_ionic_concentrations__i_Na = var_fast_sodium_current__i_Na;
        double var_ionic_concentrations__i_CaNa = var_L_type_Ca_channel__i_CaNa;
        double var_ionic_concentrations__i_Na_b = var_sodium_background_current__i_Na_b;
        double var_ionic_concentrations__i_ns_Na = var_non_specific_calcium_activated_current__i_ns_Na;
        double var_ionic_concentrations__i_NaCa = var_Na_Ca_exchanger__i_NaCa;
        double var_ionic_concentrations__i_NaK = var_sodium_potassium_pump__i_NaK;
        double var_ionic_concentrations__i_CaK = var_L_type_Ca_channel__i_CaK;
        double var_ionic_concentrations__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
        double var_ionic_concentrations__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
        double var_ionic_concentrations__i_K1 = var_time_independent_potassium_current__i_K1;
        double var_ionic_concentrations__i_Kp = var_plateau_potassium_current__i_Kp;
        double var_ionic_concentrations__i_K_Na = var_sodium_activated_potassium_current__i_K_Na;
        double var_ionic_concentrations__i_K_ATP = var_ATP_sensitive_potassium_current__i_K_ATP;
        double var_ionic_concentrations__i_ns_K = var_non_specific_calcium_activated_current__i_ns_K;
        double var_ionic_concentrations__i_to = var_transient_outward_current__i_to;
        double d_dt_membrane__V = ((-1.0) / var_membrane__Cm) * (var_membrane__i_Na + var_membrane__i_Ca_L + var_membrane__i_Ca_T + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_K_Na + var_membrane__i_K_ATP + var_membrane__i_to + var_membrane__i_K1 + var_membrane__i_Kp + var_membrane__i_NaCa + var_membrane__i_p_Ca + var_membrane__i_Na_b + var_membrane__i_Ca_b + var_membrane__i_NaK + var_membrane__i_ns_Ca + var_membrane__I_st);
        double d_dt_fast_sodium_current_m_gate__m = (var_fast_sodium_current_m_gate__alpha_m * (1.0 - var_fast_sodium_current_m_gate__m)) - (var_fast_sodium_current_m_gate__beta_m * var_fast_sodium_current_m_gate__m);
        double d_dt_fast_sodium_current_h_gate__h = (var_fast_sodium_current_h_gate__alpha_h * (1.0 - var_fast_sodium_current_h_gate__h)) - (var_fast_sodium_current_h_gate__beta_h * var_fast_sodium_current_h_gate__h);
        double d_dt_fast_sodium_current_j_gate__j = (var_fast_sodium_current_j_gate__alpha_j * (1.0 - var_fast_sodium_current_j_gate__j)) - (var_fast_sodium_current_j_gate__beta_j * var_fast_sodium_current_j_gate__j);
        double d_dt_L_type_Ca_channel_d_gate__d = (var_L_type_Ca_channel_d_gate__alpha_d * (1.0 - var_L_type_Ca_channel_d_gate__d)) - (var_L_type_Ca_channel_d_gate__beta_d * var_L_type_Ca_channel_d_gate__d);
        double d_dt_L_type_Ca_channel_f_gate__f = (var_L_type_Ca_channel_f_gate__alpha_f * (1.0 - var_L_type_Ca_channel_f_gate__f)) - (var_L_type_Ca_channel_f_gate__beta_f * var_L_type_Ca_channel_f_gate__f);
        double d_dt_T_type_Ca_channel_b_gate__b = (var_T_type_Ca_channel_b_gate__b_inf - var_T_type_Ca_channel_b_gate__b) / var_T_type_Ca_channel_b_gate__tau_b;
        double d_dt_T_type_Ca_channel_g_gate__g = (var_T_type_Ca_channel_g_gate__g_inf - var_T_type_Ca_channel_g_gate__g) / var_T_type_Ca_channel_g_gate__tau_g;
        double d_dt_rapid_delayed_rectifier_potassium_current_xr_gate__xr = (var_rapid_delayed_rectifier_potassium_current_xr_gate__xr_infinity - var_rapid_delayed_rectifier_potassium_current_xr_gate__xr) / var_rapid_delayed_rectifier_potassium_current_xr_gate__tau_xr;
        double d_dt_slow_delayed_rectifier_potassium_current_xs1_gate__xs1 = (var_slow_delayed_rectifier_potassium_current_xs1_gate__xs1_infinity - var_slow_delayed_rectifier_potassium_current_xs1_gate__xs1) / var_slow_delayed_rectifier_potassium_current_xs1_gate__tau_xs1;
        double d_dt_slow_delayed_rectifier_potassium_current_xs2_gate__xs2 = (var_slow_delayed_rectifier_potassium_current_xs2_gate__xs2_infinity - var_slow_delayed_rectifier_potassium_current_xs2_gate__xs2) / var_slow_delayed_rectifier_potassium_current_xs2_gate__tau_xs2;
        double d_dt_transient_outward_current_zdv_gate__zdv = (var_transient_outward_current_zdv_gate__zdv_ss - var_transient_outward_current_zdv_gate__zdv) / var_transient_outward_current_zdv_gate__tau_zdv;
        double d_dt_transient_outward_current_ydv_gate__ydv = (var_transient_outward_current_ydv_gate__ydv_ss - var_transient_outward_current_ydv_gate__ydv) / var_transient_outward_current_ydv_gate__tau_ydv;
        double d_dt_calcium_dynamics__APtrack = (var_calcium_dynamics__dVdt > 150000.0) ? ((100000.0 * (1.0 - var_calcium_dynamics__APtrack)) - (500.0 * var_calcium_dynamics__APtrack)) : ((-500.0) * var_calcium_dynamics__APtrack);
        double d_dt_calcium_dynamics__APtrack2 = ((var_calcium_dynamics__APtrack < 0.2) && (var_calcium_dynamics__APtrack > 0.18)) ? ((100000.0 * (1.0 - var_calcium_dynamics__APtrack2)) - (500.0 * var_calcium_dynamics__APtrack2)) : ((-500.0) * var_calcium_dynamics__APtrack2);
        double d_dt_calcium_dynamics__APtrack3 = ((var_calcium_dynamics__APtrack < 0.2) && (var_calcium_dynamics__APtrack > 0.18)) ? ((100000.0 * (1.0 - var_calcium_dynamics__APtrack3)) - (500.0 * var_calcium_dynamics__APtrack3)) : ((-10.0) * var_calcium_dynamics__APtrack3);
        double d_dt_calcium_dynamics__Cainfluxtrack = (var_calcium_dynamics__APtrack > 0.2) ? (((-var_calcium_dynamics__A_cap) * (((var_calcium_dynamics__i_CaCa + var_calcium_dynamics__i_Ca_T) - var_calcium_dynamics__i_NaCa) + var_calcium_dynamics__i_p_Ca + var_calcium_dynamics__i_Ca_b)) / (2.0 * var_calcium_dynamics__V_myo * var_calcium_dynamics__F)) : ((var_calcium_dynamics__APtrack2 > 0.01) && (var_calcium_dynamics__APtrack <= 0.2)) ? 0.0 : ((-500.0) * var_calcium_dynamics__Cainfluxtrack);
        double d_dt_calcium_dynamics__OVRLDtrack = (((1.0 / (1.0 + (var_calcium_dynamics__K_mCSQN / var_calcium_dynamics__Ca_JSR))) > var_calcium_dynamics__CSQNthresh) && (var_calcium_dynamics__OVRLDtrack3 < 0.37) && (var_calcium_dynamics__APtrack3 < 0.37)) ? (50000.0 * (1.0 - var_calcium_dynamics__OVRLDtrack)) : ((-500.0) * var_calcium_dynamics__OVRLDtrack);
        double d_dt_calcium_dynamics__OVRLDtrack2 = ((var_calcium_dynamics__OVRLDtrack > var_calcium_dynamics__Logicthresh) && (var_calcium_dynamics__OVRLDtrack2 < var_calcium_dynamics__Logicthresh)) ? (50000.0 * (1.0 - var_calcium_dynamics__OVRLDtrack2)) : ((-500.0) * var_calcium_dynamics__OVRLDtrack2);
        double d_dt_calcium_dynamics__OVRLDtrack3 = ((var_calcium_dynamics__OVRLDtrack > var_calcium_dynamics__Logicthresh) && (var_calcium_dynamics__OVRLDtrack3 < var_calcium_dynamics__Logicthresh)) ? (50000.0 * (1.0 - var_calcium_dynamics__OVRLDtrack3)) : ((-10.0) * var_calcium_dynamics__OVRLDtrack3);
        double d_dt_calcium_dynamics__Ca_JSR = (1.0 / (1.0 + ((var_calcium_dynamics__CSQN_max * var_calcium_dynamics__K_mCSQN) / pow(var_calcium_dynamics__K_mCSQN + var_calcium_dynamics__Ca_JSR, 2.0)))) * (var_calcium_dynamics__i_tr - var_calcium_dynamics__i_rel);
        double d_dt_calcium_dynamics__Ca_NSR = ((((-var_calcium_dynamics__i_tr) * var_calcium_dynamics__V_JSR) / var_calcium_dynamics__V_NSR) - var_calcium_dynamics__i_leak) + var_calcium_dynamics__i_up;
        double d_dt_calcium_dynamics__Cai = (1.0 / (1.0 + ((var_calcium_dynamics__CMDN_max * var_calcium_dynamics__K_mCMDN) / pow(var_calcium_dynamics__K_mCMDN + var_calcium_dynamics__Cai, 2.0)) + ((var_calcium_dynamics__Tn_max * var_calcium_dynamics__K_mTn) / pow(var_calcium_dynamics__K_mTn + var_calcium_dynamics__Cai, 2.0)))) * ((((-var_calcium_dynamics__A_cap) * (((var_calcium_dynamics__i_CaCa + var_calcium_dynamics__i_Ca_T) - (2.0 * var_calcium_dynamics__i_NaCa)) + var_calcium_dynamics__i_p_Ca + var_calcium_dynamics__i_Ca_b)) / (2.0 * var_calcium_dynamics__V_myo * var_calcium_dynamics__F)) + ((var_calcium_dynamics__i_rel * var_calcium_dynamics__V_JSR) / var_calcium_dynamics__V_myo) + (((var_calcium_dynamics__i_leak - var_calcium_dynamics__i_up) * var_calcium_dynamics__V_NSR) / var_calcium_dynamics__V_myo));
        double d_dt_ionic_concentrations__Nai = ((-(var_ionic_concentrations__i_Na + var_ionic_concentrations__i_CaNa + var_ionic_concentrations__i_Na_b + var_ionic_concentrations__i_ns_Na + (var_ionic_concentrations__i_NaCa * 3.0) + (var_ionic_concentrations__i_NaK * 3.0))) * var_ionic_concentrations__A_cap) / (var_ionic_concentrations__V_myo * var_ionic_concentrations__F);
        double d_dt_ionic_concentrations__Ki = ((-(var_ionic_concentrations__i_CaK + var_ionic_concentrations__i_Kr + var_ionic_concentrations__i_Ks + var_ionic_concentrations__i_K1 + var_ionic_concentrations__i_Kp + var_ionic_concentrations__i_K_Na + var_ionic_concentrations__i_K_ATP + var_ionic_concentrations__i_to + var_ionic_concentrations__i_ns_K + ((-var_ionic_concentrations__i_NaK) * 2.0))) * var_ionic_concentrations__A_cap) / (var_ionic_concentrations__V_myo * var_ionic_concentrations__F);

        rDY[0] = d_dt_membrane__V*1e-3;
        rDY[1] = d_dt_fast_sodium_current_m_gate__m*1e-3;
        rDY[2] = d_dt_fast_sodium_current_h_gate__h*1e-3;
        rDY[3] = d_dt_fast_sodium_current_j_gate__j*1e-3;
        rDY[4] = d_dt_L_type_Ca_channel_d_gate__d*1e-3;
        rDY[5] = d_dt_L_type_Ca_channel_f_gate__f*1e-3;
        rDY[6] = d_dt_T_type_Ca_channel_b_gate__b*1e-3;
        rDY[7] = d_dt_T_type_Ca_channel_g_gate__g*1e-3;
        rDY[8] = d_dt_rapid_delayed_rectifier_potassium_current_xr_gate__xr*1e-3;
        rDY[9] = d_dt_slow_delayed_rectifier_potassium_current_xs1_gate__xs1*1e-3;
        rDY[10] = d_dt_slow_delayed_rectifier_potassium_current_xs2_gate__xs2*1e-3;
        rDY[11] = d_dt_transient_outward_current_zdv_gate__zdv*1e-3;
        rDY[12] = d_dt_transient_outward_current_ydv_gate__ydv*1e-3;
        rDY[13] = d_dt_calcium_dynamics__Cai*1e-3;
        rDY[14] = d_dt_calcium_dynamics__Ca_JSR*1e-3;
        rDY[15] = d_dt_calcium_dynamics__Ca_NSR*1e-3;
        rDY[16] = d_dt_calcium_dynamics__APtrack*1e-3;
        rDY[17] = d_dt_calcium_dynamics__APtrack2*1e-3;
        rDY[18] = d_dt_calcium_dynamics__APtrack3*1e-3;
        rDY[19] = d_dt_calcium_dynamics__Cainfluxtrack*1e-3;
        rDY[20] = d_dt_calcium_dynamics__OVRLDtrack*1e-3;
        rDY[21] = d_dt_calcium_dynamics__OVRLDtrack2*1e-3;
        rDY[22] = d_dt_calcium_dynamics__OVRLDtrack3*1e-3;
        rDY[23] = d_dt_ionic_concentrations__Nai*1e-3;
        rDY[24] = d_dt_ionic_concentrations__Ki*1e-3;
    }

};



template<>
void OdeSystemInformation<FaberRudy2000Version3>::Initialise(void)
{
    // Time units: second
    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("millivolt");
    this->mInitialConditions.push_back(-90);

    this->mVariableNames.push_back("m");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0008);

    this->mVariableNames.push_back("h");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.993771);

    this->mVariableNames.push_back("j");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.995727);

    this->mVariableNames.push_back("d");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(3.210618e-6);

    this->mVariableNames.push_back("f");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.999837);

    this->mVariableNames.push_back("b");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.000970231);

    this->mVariableNames.push_back("g");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.994305);

    this->mVariableNames.push_back("xr");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.000124042);

    this->mVariableNames.push_back("xs1");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.00445683);

    this->mVariableNames.push_back("xs2");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.00445683);

    this->mVariableNames.push_back("zdv");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.5);

    this->mVariableNames.push_back("ydv");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.5);

    this->mVariableNames.push_back("CaI");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(6e-5);

    this->mVariableNames.push_back("Ca_JSR");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(1.8);

    this->mVariableNames.push_back("Ca_NSR");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(1.8);

    this->mVariableNames.push_back("APtrack");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("APtrack2");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("APtrack3");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("Cainfluxtrack");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("OVRLDtrack");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("OVRLDtrack2");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("OVRLDtrack3");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("Nai");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(9);

    this->mVariableNames.push_back("Ki");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(141.2);

    this->mInitialised = true;
}



#endif
