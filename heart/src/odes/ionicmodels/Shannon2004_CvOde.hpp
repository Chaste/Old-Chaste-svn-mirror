#ifdef CHASTE_CVODE

#ifndef _CVODECELLSHANNON2004FROMCELLML_
#define _CVODECELLSHANNON2004FROMCELLML_

// Model: shannon_wang_puglisi_weber_bers_2004_model_updated
// Processed by pycml - CellML Tools in Python
//     (translate: 6039, pycml: 5150)
// on Thu Jan 21 12:50:32 2010

#include <cmath>
#include <cassert>
#include "AbstractStimulusFunction.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractCvodeCell.hpp"

class CvOdeCellShannon2004FromCellML : public AbstractCvodeCell
{
public:
    CvOdeCellShannon2004FromCellML(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
                                   boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
        : AbstractCvodeCell(pOdeSolver, 45, 0, pIntracellularStimulus)
    {
        // Time units: millisecond
        //
        mpSystemInfo = OdeSystemInformation<CvOdeCellShannon2004FromCellML>::Instance();
        Init();
    }

    ~CvOdeCellShannon2004FromCellML(void)
    {
    }

    void EvaluateRhs(realtype var_environment__time, N_Vector rY, N_Vector ydot)
    {
        // Inputs:
        // Time units: millisecond
        double var_cell__V = NV_Ith_S(rY,0);
        // Units: millivolt; Initial value: -85.719687955637
        double var_INa_h_gate__h = NV_Ith_S(rY, 1);
        // Units: dimensionless; Initial value: 0.987140350343
        double var_INa_j_gate__j = NV_Ith_S(rY, 2);
        // Units: dimensionless; Initial value: 0.991822731369
        double var_INa_m_gate__m = NV_Ith_S(rY, 3);
        // Units: dimensionless; Initial value: 0.001370685156
        double var_IKr_Xr_gate__Xr = NV_Ith_S(rY, 4);
        // Units: dimensionless; Initial value: 0.008471550841
        double var_IKs_Xs_gate__Xs = NV_Ith_S(rY, 5);
        // Units: dimensionless; Initial value: 0.00687399199
        double var_Itos_X_gate__X_tos = NV_Ith_S(rY, 6);
        // Units: dimensionless; Initial value: 0.004011272375
        double var_Itos_Y_gate__Y_tos = NV_Ith_S(rY, 7);
        // Units: dimensionless; Initial value: 0.293519921626
        double var_Itos_R_gate__R_tos = NV_Ith_S(rY, 8);
        // Units: dimensionless; Initial value: 0.383430556383
        double var_Itof_X_gate__X_tof = NV_Ith_S(rY, 9);
        // Units: dimensionless; Initial value: 0.00401120993
        double var_Itof_Y_gate__Y_tof = NV_Ith_S(rY, 10);
        // Units: dimensionless; Initial value: 0.9946314893
        double var_ICaL_d_gate__d = NV_Ith_S(rY, 11);
        // Units: dimensionless; Initial value: 0.000006997531
        double var_ICaL_f_gate__f = NV_Ith_S(rY, 12);
        // Units: dimensionless; Initial value: 1.000675515962
        double var_ICaL_fCa_gate__fCaB_SL = NV_Ith_S(rY, 13);
        // Units: dimensionless; Initial value: 0.015352888928
        double var_ICaL_fCa_gate__fCaB_jct = NV_Ith_S(rY, 14);
        // Units: dimensionless; Initial value: 0.024609183734
        double var_Jrel_SR__R = NV_Ith_S(rY, 15);
        // Units: dimensionless; Initial value: 0.884673513138
        double var_Jrel_SR__I = NV_Ith_S(rY, 16);
        // Units: dimensionless; Initial value: 0.00000009272
        double var_Jrel_SR__O = NV_Ith_S(rY, 17);
        // Units: dimensionless; Initial value: 0.000000711264
        double var_Na_buffer__Na_SL = NV_Ith_S(rY, 18);
        // Units: millimolar; Initial value: 8.874077316753
        double var_Na_buffer__Na_jct = NV_Ith_S(rY, 19);
        // Units: millimolar; Initial value: 8.872823559072
        double var_Na_buffer__Na_SL_buf = NV_Ith_S(rY, 20);
        // Units: millimolar; Initial value: 0.776121392467
        double var_Na_buffer__Na_jct_buf = NV_Ith_S(rY, 21);
        // Units: millimolar; Initial value: 3.557055389701
        double var_Na_buffer__Nai = NV_Ith_S(rY, 22);
        // Units: millimolar; Initial value: 8.874461106492
        double var_Ca_buffer__Ca_SR = NV_Ith_S(rY, 23);
        // Units: millimolar; Initial value: 0.545611267699
        double var_Ca_buffer__Ca_SL = NV_Ith_S(rY, 24);
        // Units: millimolar; Initial value: 0.000106395937
        double var_Ca_buffer__Ca_jct = NV_Ith_S(rY, 25);
        // Units: millimolar; Initial value: 0.000174843061
        double var_Ca_buffer__Cai = NV_Ith_S(rY, 26);
        // Units: millimolar; Initial value: 0.000087350002
        double var_Ca_buffer__Ca_SLB_SL = NV_Ith_S(rY, 27);
        // Units: millimolar; Initial value: 0.009868629147
        double var_Ca_buffer__Ca_SLB_jct = NV_Ith_S(rY, 28);
        // Units: millimolar; Initial value: 0.007780801995
        double var_Ca_buffer__Ca_SLHigh_SL = NV_Ith_S(rY, 29);
        // Units: millimolar; Initial value: 0.114438990328
        double var_Ca_buffer__Ca_SLHigh_jct = NV_Ith_S(rY, 30);
        // Units: millimolar; Initial value: 0.077503874257
        double var_Ca_buffer__Ca_Calsequestrin = NV_Ith_S(rY, 31);
        // Units: millimolar; Initial value: 1.186496899338
        double var_cytosolic_Ca_buffer__Ca_TroponinC = NV_Ith_S(rY, 32);
        // Units: millimolar; Initial value: 0.008963736337
        double var_cytosolic_Ca_buffer__Ca_TroponinC_Ca_Mg = NV_Ith_S(rY, 33);
        // Units: millimolar; Initial value: 0.117995194438
        double var_cytosolic_Ca_buffer__Mg_TroponinC_Ca_Mg = NV_Ith_S(rY, 34);
        // Units: millimolar; Initial value: 0.010337654274
        double var_cytosolic_Ca_buffer__Ca_Calmodulin = NV_Ith_S(rY, 35);
        // Units: millimolar; Initial value: 0.000295961245
        double var_cytosolic_Ca_buffer__Ca_Myosin = NV_Ith_S(rY, 36);
        // Units: millimolar; Initial value: 0.001984672275
        double var_cytosolic_Ca_buffer__Mg_Myosin = NV_Ith_S(rY, 37);
        // Units: millimolar; Initial value: 0.137497736234
        double var_cytosolic_Ca_buffer__Ca_SRB = NV_Ith_S(rY, 38);
        // Units: millimolar; Initial value: 0.002177112381
        double var_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_Cytosol = NV_Ith_S(rY, 39);
        // Units: millimolar; Initial value: 0
        double var_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_SL = NV_Ith_S(rY, 40);
        // Units: millimolar; Initial value: 0
        double var_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_jct = NV_Ith_S(rY, 41);
        // Units: millimolar; Initial value: 0
        double var_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_Cytosol = NV_Ith_S(rY, 42);
        // Units: millimolar; Initial value: 0
        double var_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_SL = NV_Ith_S(rY, 43);
        // Units: millimolar; Initial value: 0
        double var_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_jct = NV_Ith_S(rY, 44);
        // Units: millimolar; Initial value: 0

        // Mathematics
        const double var_INa__V = var_cell__V;
        const double var_reversal_potentials__Na_jct = var_Na_buffer__Na_jct;
        const double var_model_parameters__R = 8314.3;
        const double var_reversal_potentials__R = var_model_parameters__R;
        const double var_model_parameters__Nao = 140.0;
        const double var_reversal_potentials__Nao = var_model_parameters__Nao;
        const double var_model_parameters__F = 96486.7;
        const double var_reversal_potentials__F = var_model_parameters__F;
        const double var_model_parameters__T = 310.0;
        const double var_reversal_potentials__T = var_model_parameters__T;
        const double var_reversal_potentials__E_Na_jct = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Nao / var_reversal_potentials__Na_jct);
        const double var_INa__E_Na_jct = var_reversal_potentials__E_Na_jct;
        double var_INa__G_INa = 16.0;
        const double var_INa__m = var_INa_m_gate__m;
        const double var_INa__j = var_INa_j_gate__j;
        const double var_INa__h = var_INa_h_gate__h;
        const double var_INa__openProb = pow(var_INa__m, 3.0) * var_INa__h * var_INa__j;
        const double var_INa__Fx_Na_jct = 0.11;
        const double var_INa__i_Na_jct = var_INa__Fx_Na_jct * var_INa__G_INa * var_INa__openProb * (var_INa__V - var_INa__E_Na_jct);
        const double var_reversal_potentials__Na_SL = var_Na_buffer__Na_SL;
        const double var_reversal_potentials__E_Na_SL = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Nao / var_reversal_potentials__Na_SL);
        const double var_INa__E_Na_SL = var_reversal_potentials__E_Na_SL;
        const double var_INa__Fx_Na_SL = 0.89;
        const double var_INa__i_Na_SL = var_INa__Fx_Na_SL * var_INa__G_INa * var_INa__openProb * (var_INa__V - var_INa__E_Na_SL);
        const double var_INa__i_Na = var_INa__i_Na_jct + var_INa__i_Na_SL;
        const double var_cell__i_Na = var_INa__i_Na;
        const double var_INab__Fx_NaBk_jct = 0.11;
        const double var_INab__E_Na_jct = var_reversal_potentials__E_Na_jct;
        const double var_INab__G_NaBk = 0.000297;
        const double var_INab__V = var_cell__V;
        const double var_INab__i_Nab_jct = var_INab__Fx_NaBk_jct * var_INab__G_NaBk * (var_INab__V - var_INab__E_Na_jct);
        const double var_INab__E_Na_SL = var_reversal_potentials__E_Na_SL;
        const double var_INab__Fx_NaBk_SL = 0.89;
        const double var_INab__i_Nab_SL = var_INab__Fx_NaBk_SL * var_INab__G_NaBk * (var_INab__V - var_INab__E_Na_SL);
        const double var_INab__i_Nab = var_INab__i_Nab_jct + var_INab__i_Nab_SL;
        const double var_cell__i_Nab = var_INab__i_Nab;
        const double var_INaK__V = var_cell__V;
        const double var_INaK__F = var_model_parameters__F;
        const double var_INaK__T = var_model_parameters__T;
        const double var_INaK__Nao = var_model_parameters__Nao;
        const double var_INaK__sigma = (exp(var_INaK__Nao / 67.3) - 1.0) / 7.0;
        const double var_INaK__R = var_model_parameters__R;
        const double var_INaK__f_NaK = 1.0 / (1.0 + (0.1245 * exp(((-0.1) * var_INaK__V * var_INaK__F) / (var_INaK__R * var_INaK__T))) + (0.0365 * var_INaK__sigma * exp(((-var_INaK__V) * var_INaK__F) / (var_INaK__R * var_INaK__T))));
        const double var_INaK__H_NaK = 4.0;
        const double var_INaK__Q10_NaK = 1.63;
        const double var_INaK__Q_NaK = pow(var_INaK__Q10_NaK, (var_INaK__T - 310.0) / 10.0);
        const double var_model_parameters__Ko = 5.4;
        const double var_INaK__Ko = var_model_parameters__Ko;
        const double var_INaK__Km_Nai = 11.0;
        const double var_INaK__Fx_NaK_SL = 0.89;
        const double var_INaK__Na_SL = var_Na_buffer__Na_SL;
        const double var_INaK__I_NaK_max = 1.91;
        const double var_INaK__Q10_Km_Nai = 1.49;
        const double var_INaK__Q_Km_Nai = pow(var_INaK__Q10_Km_Nai, (var_INaK__T - 310.0) / 10.0);
        const double var_INaK__Km_Ko = 1.5;
        const double var_INaK__i_NaK_SL = (((var_INaK__Fx_NaK_SL * var_INaK__Q_NaK * var_INaK__I_NaK_max * var_INaK__f_NaK) / (1.0 + pow((var_INaK__Q_Km_Nai * var_INaK__Km_Nai) / var_INaK__Na_SL, var_INaK__H_NaK))) * var_INaK__Ko) / (var_INaK__Ko + var_INaK__Km_Ko);
        const double var_INaK__Fx_NaK_jct = 0.11;
        const double var_INaK__Na_jct = var_Na_buffer__Na_jct;
        const double var_INaK__i_NaK_jct = (((var_INaK__Fx_NaK_jct * var_INaK__Q_NaK * var_INaK__I_NaK_max * var_INaK__f_NaK) / (1.0 + pow((var_INaK__Q_Km_Nai * var_INaK__Km_Nai) / var_INaK__Na_jct, var_INaK__H_NaK))) * var_INaK__Ko) / (var_INaK__Ko + var_INaK__Km_Ko);
        const double var_INaK__i_NaK = var_INaK__i_NaK_jct + var_INaK__i_NaK_SL;
        const double var_cell__i_NaK = var_INaK__i_NaK;
        const double var_IKr__Xr = var_IKr_Xr_gate__Xr;
        const double var_IKr__V = var_cell__V;
        const double var_IKr_Rr_gate__V = var_IKr__V;
        const double var_IKr_Rr_gate__Rr = 1.0 / (1.0 + exp((33.0 + var_IKr_Rr_gate__V) / 22.4));
        const double var_IKr__Rr = var_IKr_Rr_gate__Rr;
        double var_IKr__G_IKr_const = 0.03;
        const double var_IKr__Ko = var_model_parameters__Ko;
        const double var_IKr__G_IKr = var_IKr__G_IKr_const * sqrt(var_IKr__Ko / 5.4);
        const double var_model_parameters__Ki = 135.0;
        const double var_reversal_potentials__Ki = var_model_parameters__Ki;
        const double var_reversal_potentials__Ko = var_model_parameters__Ko;
        const double var_reversal_potentials__E_K = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Ko / var_reversal_potentials__Ki);
        const double var_IKr__E_K = var_reversal_potentials__E_K;
        const double var_IKr__i_Kr = var_IKr__G_IKr * var_IKr__Xr * var_IKr__Rr * (var_IKr__V - var_IKr__E_K);
        const double var_cell__i_Kr = var_IKr__i_Kr;
        const double var_IKs__Fx_Ks_jct = 0.11;
        const double var_IKs__Xs = var_IKs_Xs_gate__Xs;
        const double var_IKs__Ca_jct = var_Ca_buffer__Ca_jct;
        const double var_IKs__pCa_jct = (-log10(var_IKs__Ca_jct / 1.0)) + 3.0;
        const double var_IKs__G_Ks_jct = 0.07 * (0.057 + (0.19 / (1.0 + exp(((-7.2) + var_IKs__pCa_jct) / 0.6))));
        const double var_IKs__R = var_model_parameters__R;
        const double var_IKs__T = var_model_parameters__T;
        const double var_IKs__pKNa = 0.01833;
        const double var_IKs__F = var_model_parameters__F;
        const double var_IKs__Nao = var_model_parameters__Nao;
        const double var_IKs__Ki = var_model_parameters__Ki;
        const double var_IKs__Ko = var_model_parameters__Ko;
        const double var_IKs__Nai = var_Na_buffer__Nai;
        const double var_IKs__E_Ks = ((var_IKs__R * var_IKs__T) / var_IKs__F) * log((var_IKs__Ko + (var_IKs__pKNa * var_IKs__Nao)) / (var_IKs__Ki + (var_IKs__pKNa * var_IKs__Nai)));
        const double var_IKs__V = var_cell__V;
        const double var_IKs__i_Ks_jct = var_IKs__Fx_Ks_jct * var_IKs__G_Ks_jct * pow(var_IKs__Xs, 2.0) * (var_IKs__V - var_IKs__E_Ks);
        const double var_IKs__Fx_Ks_SL = 0.89;
        const double var_IKs__Ca_SL = var_Ca_buffer__Ca_SL;
        const double var_IKs__pCa_SL = (-log10(var_IKs__Ca_SL / 1.0)) + 3.0;
        const double var_IKs__G_Ks_SL = 0.07 * (0.057 + (0.19 / (1.0 + exp(((-7.2) + var_IKs__pCa_SL) / 0.6))));
        const double var_IKs__i_Ks_SL = var_IKs__Fx_Ks_SL * var_IKs__G_Ks_SL * pow(var_IKs__Xs, 2.0) * (var_IKs__V - var_IKs__E_Ks);
        const double var_IKs__i_Ks = var_IKs__i_Ks_jct + var_IKs__i_Ks_SL;
        const double var_cell__i_Ks = var_IKs__i_Ks;
        const double var_Itos__Y_tos = var_Itos_Y_gate__Y_tos;
        const double var_Itos__X_tos = var_Itos_X_gate__X_tos;
        const double var_Itos__R_tos = var_Itos_R_gate__R_tos;
        const double var_Itos__G_tos = 0.06;
        const double var_Itos__E_K = var_reversal_potentials__E_K;
        const double var_Itos__V = var_cell__V;
        const double var_Itos__i_tos = var_Itos__G_tos * var_Itos__X_tos * (var_Itos__Y_tos + (0.5 * var_Itos__R_tos)) * (var_Itos__V - var_Itos__E_K);
        const double var_cell__i_tos = var_Itos__i_tos;
        const double var_Itof__G_tof = 0.02;
        const double var_Itof__E_K = var_reversal_potentials__E_K;
        const double var_Itof__X_tof = var_Itof_X_gate__X_tof;
        const double var_Itof__Y_tof = var_Itof_Y_gate__Y_tof;
        const double var_Itof__V = var_cell__V;
        const double var_Itof__i_tof = var_Itof__G_tof * var_Itof__X_tof * var_Itof__Y_tof * (var_Itof__V - var_Itof__E_K);
        const double var_cell__i_tof = var_Itof__i_tof;
        const double var_IK1__V = var_cell__V;
        const double var_IK1_K1_gate__V = var_IK1__V;
        const double var_IK1__E_K = var_reversal_potentials__E_K;
        const double var_IK1_K1_gate__E_K = var_IK1__E_K;
        const double var_IK1_K1_gate__beta_K1 = ((0.49124 * exp(0.08032 * ((var_IK1_K1_gate__V - var_IK1_K1_gate__E_K) + 5.476))) + (1.0 * exp(0.06175 * (var_IK1_K1_gate__V - (var_IK1_K1_gate__E_K + 594.31))))) / (1.0 + exp((-0.5143) * ((var_IK1_K1_gate__V - var_IK1_K1_gate__E_K) + 4.753)));
        const double var_IK1_K1_gate__alpha_K1 = 1.02 / (1.0 + exp(0.2385 * (var_IK1_K1_gate__V - (var_IK1_K1_gate__E_K + 59.215))));
        const double var_IK1_K1_gate__K1_infinity = var_IK1_K1_gate__alpha_K1 / (var_IK1_K1_gate__alpha_K1 + var_IK1_K1_gate__beta_K1);
        const double var_IK1__K1_infinity = var_IK1_K1_gate__K1_infinity;
        const double var_IK1__Ko = var_model_parameters__Ko;
        const double var_IK1__G_K1 = 0.9 * sqrt(var_IK1__Ko / 5.4);
        const double var_IK1__i_K1 = var_IK1__G_K1 * var_IK1__K1_infinity * (var_IK1__V - var_IK1__E_K);
        const double var_cell__i_K1 = var_IK1__i_K1;
        const double var_INaCa__K_mNai = 12.29;
        const double var_INaCa__HNa = 3.0;
        const double var_INaCa__K_mNao = 87.5;
        const double var_INaCa__T = var_model_parameters__T;
        const double var_model_parameters__Cao = 1.8;
        const double var_INaCa__Cao = var_model_parameters__Cao;
        const double var_INaCa__V = var_cell__V;
        const double var_INaCa__R = var_model_parameters__R;
        const double var_INaCa__Nao = var_model_parameters__Nao;
        const double var_INaCa__eta = 0.35;
        const double var_INaCa__Ca_jct = var_Ca_buffer__Ca_jct;
        const double var_INaCa__F = var_model_parameters__F;
        const double var_INaCa__Na_jct = var_Na_buffer__Na_jct;
        const double var_INaCa__ksat = 0.27;
        const double var_INaCa__temp_jct = ((exp((var_INaCa__eta * var_INaCa__V * var_INaCa__F) / (var_INaCa__R * var_INaCa__T)) * pow(var_INaCa__Na_jct, var_INaCa__HNa) * var_INaCa__Cao) - (exp(((var_INaCa__eta - 1.0) * var_INaCa__V * var_INaCa__F) / (var_INaCa__R * var_INaCa__T)) * pow(var_INaCa__Nao, var_INaCa__HNa) * var_INaCa__Ca_jct)) / (1.0 + (var_INaCa__ksat * exp(((var_INaCa__eta - 1.0) * var_INaCa__V * var_INaCa__F) / (var_INaCa__R * var_INaCa__T))));
        const double var_INaCa__K_mCao = 1.3;
        const double var_INaCa__V_max = 9.0;
        const double var_INaCa__K_mCai = 0.00359;
        const double var_INaCa__Fx_NCX_jct = 0.11;
        const double var_INaCa__Kd_act = 0.000256;
        const double var_INaCa__Ka_jct = 1.0 / (1.0 + pow(var_INaCa__Kd_act / var_INaCa__Ca_jct, 3.0));
        const double var_INaCa__Q10_NCX = 1.57;
        const double var_INaCa__Q_NCX = pow(var_INaCa__Q10_NCX, (var_INaCa__T - 310.0) / 10.0);
        const double var_INaCa__i_NaCa_jct = (var_INaCa__Fx_NCX_jct * var_INaCa__V_max * var_INaCa__Ka_jct * var_INaCa__Q_NCX * var_INaCa__temp_jct) / ((var_INaCa__K_mCai * pow(var_INaCa__Nao, var_INaCa__HNa) * (1.0 + pow(var_INaCa__Na_jct / var_INaCa__K_mNai, var_INaCa__HNa))) + (pow(var_INaCa__K_mNao, var_INaCa__HNa) * var_INaCa__Ca_jct * (1.0 + (var_INaCa__Ca_jct / var_INaCa__K_mCai))) + (var_INaCa__K_mCao * pow(var_INaCa__Na_jct, var_INaCa__HNa)) + (pow(var_INaCa__Na_jct, var_INaCa__HNa) * var_INaCa__Cao) + (pow(var_INaCa__Nao, var_INaCa__HNa) * var_INaCa__Ca_jct));
        const double var_INaCa__Fx_NCX_SL = 0.89;
        const double var_INaCa__Na_SL = var_Na_buffer__Na_SL;
        const double var_INaCa__Ca_SL = var_Ca_buffer__Ca_SL;
        const double var_INaCa__temp_SL = ((exp((var_INaCa__eta * var_INaCa__V * var_INaCa__F) / (var_INaCa__R * var_INaCa__T)) * pow(var_INaCa__Na_SL, var_INaCa__HNa) * var_INaCa__Cao) - (exp(((var_INaCa__eta - 1.0) * var_INaCa__V * var_INaCa__F) / (var_INaCa__R * var_INaCa__T)) * pow(var_INaCa__Nao, var_INaCa__HNa) * var_INaCa__Ca_SL)) / (1.0 + (var_INaCa__ksat * exp(((var_INaCa__eta - 1.0) * var_INaCa__V * var_INaCa__F) / (var_INaCa__R * var_INaCa__T))));
        const double var_INaCa__Ka_SL = 1.0 / (1.0 + pow(var_INaCa__Kd_act / var_INaCa__Ca_SL, 3.0));
        const double var_INaCa__i_NaCa_SL = (var_INaCa__Fx_NCX_SL * var_INaCa__V_max * var_INaCa__Ka_SL * var_INaCa__Q_NCX * var_INaCa__temp_SL) / ((var_INaCa__K_mCai * pow(var_INaCa__Nao, var_INaCa__HNa) * (1.0 + pow(var_INaCa__Na_SL / var_INaCa__K_mNai, var_INaCa__HNa))) + (pow(var_INaCa__K_mNao, var_INaCa__HNa) * var_INaCa__Ca_SL * (1.0 + (var_INaCa__Ca_SL / var_INaCa__K_mCai))) + (var_INaCa__K_mCao * pow(var_INaCa__Na_SL, var_INaCa__HNa)) + (pow(var_INaCa__Na_SL, var_INaCa__HNa) * var_INaCa__Cao) + (pow(var_INaCa__Nao, var_INaCa__HNa) * var_INaCa__Ca_SL));
        const double var_INaCa__i_NaCa = var_INaCa__i_NaCa_jct + var_INaCa__i_NaCa_SL;
        const double var_cell__i_NaCa = var_INaCa__i_NaCa;
        const double var_ICl_Ca__V = var_cell__V;
        const double var_ICl_Ca__G_Cl = 0.109625;
        const double var_ICl_Ca__Ca_jct = var_Ca_buffer__Ca_jct;
        const double var_ICl_Ca__Kd_ClCa = 0.1;
        const double var_ICl_Ca__Fx_Cl_jct = 0.11;
        const double var_model_parameters__Cli = 15.0;
        const double var_reversal_potentials__Cli = var_model_parameters__Cli;
        const double var_model_parameters__Clo = 150.0;
        const double var_reversal_potentials__Clo = var_model_parameters__Clo;
        const double var_reversal_potentials__E_Cl = (((-var_reversal_potentials__R) * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Clo / var_reversal_potentials__Cli);
        const double var_ICl_Ca__E_Cl = var_reversal_potentials__E_Cl;
        const double var_ICl_Ca__Ca_SL = var_Ca_buffer__Ca_SL;
        const double var_ICl_Ca__Fx_Cl_SL = 0.89;
        const double var_ICl_Ca__i_Cl_Ca = var_ICl_Ca__G_Cl * (var_ICl_Ca__V - var_ICl_Ca__E_Cl) * ((var_ICl_Ca__Fx_Cl_jct / (1.0 + (var_ICl_Ca__Kd_ClCa / var_ICl_Ca__Ca_jct))) + (var_ICl_Ca__Fx_Cl_SL / (1.0 + (var_ICl_Ca__Kd_ClCa / var_ICl_Ca__Ca_SL))));
        const double var_cell__i_Cl_Ca = var_ICl_Ca__i_Cl_Ca;
        const double var_IClb__V = var_cell__V;
        const double var_IClb__G_ClBk = 0.009;
        const double var_IClb__E_Cl = var_reversal_potentials__E_Cl;
        const double var_IClb__i_Clb = var_IClb__G_ClBk * (var_IClb__V - var_IClb__E_Cl);
        const double var_cell__i_Clb = var_IClb__i_Clb;
        const double var_ICaL__T = var_model_parameters__T;
        const double var_ICaL__V = var_cell__V;
        const double var_ICaL__Nao = var_model_parameters__Nao;
        const double var_ICaL__Fx_ICaL_SL = 0.1;
        const double var_ICaL__R = var_model_parameters__R;
        const double var_ICaL__Na_SL = var_Na_buffer__Na_SL;
        const double var_ICaL__d = var_ICaL_d_gate__d;
        const double var_ICaL__f = var_ICaL_f_gate__f;
        const double var_ICaL__Q10_CaL = 1.8;
        const double var_ICaL__Q_CaL = pow(var_ICaL__Q10_CaL, (var_ICaL__T - 310.0) / 10.0);
        const double var_ICaL__F = var_model_parameters__F;
        const double var_ICaL__temp = (0.45 * var_ICaL__d * var_ICaL__f * var_ICaL__Q_CaL * var_ICaL__V * pow(var_ICaL__F, 2.0)) / (var_ICaL__R * var_ICaL__T);
        const double var_ICaL__gamma_Nao = 0.75;
        const double var_ICaL__gamma_Nai = 0.75;
        const double var_ICaL_fCa_gate__fCa_SL = 1.0 - var_ICaL_fCa_gate__fCaB_SL;
        const double var_ICaL__fCa_SL = var_ICaL_fCa_gate__fCa_SL;
        const double var_ICaL__PNa = 1.5e-08;
        const double var_ICaL__i_CaL_Na_SL = (var_ICaL__temp * var_ICaL__fCa_SL * var_ICaL__Fx_ICaL_SL * var_ICaL__PNa * ((var_ICaL__gamma_Nai * var_ICaL__Na_SL * exp((var_ICaL__V * var_ICaL__F) / (var_ICaL__R * var_ICaL__T))) - (var_ICaL__gamma_Nao * var_ICaL__Nao))) / (exp((var_ICaL__V * var_ICaL__F) / (var_ICaL__R * var_ICaL__T)) - 1.0);
        const double var_ICaL__Ki = var_model_parameters__Ki;
        const double var_ICaL__PK = 2.7e-07;
        const double var_ICaL_fCa_gate__fCa_jct = 1.0 - var_ICaL_fCa_gate__fCaB_jct;
        const double var_ICaL__fCa_jct = var_ICaL_fCa_gate__fCa_jct;
        const double var_ICaL__Ko = var_model_parameters__Ko;
        const double var_ICaL__Fx_ICaL_jct = 0.9;
        const double var_ICaL__gamma_Ko = 0.75;
        const double var_ICaL__gamma_Ki = 0.75;
        const double var_ICaL__i_CaL_K = (var_ICaL__temp * ((var_ICaL__fCa_SL * var_ICaL__Fx_ICaL_SL) + (var_ICaL__fCa_jct * var_ICaL__Fx_ICaL_jct)) * var_ICaL__PK * ((var_ICaL__gamma_Ki * var_ICaL__Ki * exp((var_ICaL__V * var_ICaL__F) / (var_ICaL__R * var_ICaL__T))) - (var_ICaL__gamma_Ko * var_ICaL__Ko))) / (exp((var_ICaL__V * var_ICaL__F) / (var_ICaL__R * var_ICaL__T)) - 1.0);
        const double var_ICaL__Ca_jct = var_Ca_buffer__Ca_jct;
        const double var_ICaL__gamma_Cao = 0.341;
        const double var_ICaL__gamma_Cai = 0.341;
        const double var_ICaL__Cao = var_model_parameters__Cao;
        const double var_ICaL__PCa = 0.00054;
        const double var_ICaL__i_CaL_Ca_jct = (var_ICaL__temp * var_ICaL__fCa_jct * var_ICaL__Fx_ICaL_jct * var_ICaL__PCa * 4.0 * ((var_ICaL__gamma_Cai * var_ICaL__Ca_jct * exp((2.0 * var_ICaL__V * var_ICaL__F) / (var_ICaL__R * var_ICaL__T))) - (var_ICaL__gamma_Cao * var_ICaL__Cao))) / (exp((2.0 * var_ICaL__V * var_ICaL__F) / (var_ICaL__R * var_ICaL__T)) - 1.0);
        const double var_ICaL__Na_jct = var_Na_buffer__Na_jct;
        const double var_ICaL__i_CaL_Na_jct = (var_ICaL__temp * var_ICaL__fCa_jct * var_ICaL__Fx_ICaL_jct * var_ICaL__PNa * ((var_ICaL__gamma_Nai * var_ICaL__Na_jct * exp((var_ICaL__V * var_ICaL__F) / (var_ICaL__R * var_ICaL__T))) - (var_ICaL__gamma_Nao * var_ICaL__Nao))) / (exp((var_ICaL__V * var_ICaL__F) / (var_ICaL__R * var_ICaL__T)) - 1.0);
        double var_ICaL__G_CaL_mult = 1.0;
        const double var_ICaL__Ca_SL = var_Ca_buffer__Ca_SL;
        const double var_ICaL__i_CaL_Ca_SL = (var_ICaL__temp * var_ICaL__fCa_SL * var_ICaL__Fx_ICaL_SL * var_ICaL__PCa * 4.0 * ((var_ICaL__gamma_Cai * var_ICaL__Ca_SL * exp((2.0 * var_ICaL__V * var_ICaL__F) / (var_ICaL__R * var_ICaL__T))) - (var_ICaL__gamma_Cao * var_ICaL__Cao))) / (exp((2.0 * var_ICaL__V * var_ICaL__F) / (var_ICaL__R * var_ICaL__T)) - 1.0);
        const double var_ICaL__i_CaL = var_ICaL__G_CaL_mult * (var_ICaL__i_CaL_Ca_SL + var_ICaL__i_CaL_Ca_jct + var_ICaL__i_CaL_Na_SL + var_ICaL__i_CaL_Na_jct + var_ICaL__i_CaL_K);
        const double var_cell__i_CaL = var_ICaL__i_CaL;
        const double var_ICab__G_CaBk = 0.0002513;
        const double var_ICab__V = var_cell__V;
        const double var_ICab__Fx_CaBk_SL = 0.89;
        const double var_reversal_potentials__Ca_SL = var_Ca_buffer__Ca_SL;
        const double var_reversal_potentials__Cao = var_model_parameters__Cao;
        const double var_reversal_potentials__E_Ca_SL = ((var_reversal_potentials__R * var_reversal_potentials__T) / (2.0 * var_reversal_potentials__F)) * log(var_reversal_potentials__Cao / var_reversal_potentials__Ca_SL);
        const double var_ICab__E_Ca_SL = var_reversal_potentials__E_Ca_SL;
        const double var_ICab__i_Cab_SL = var_ICab__G_CaBk * var_ICab__Fx_CaBk_SL * (var_ICab__V - var_ICab__E_Ca_SL);
        const double var_ICab__Fx_CaBk_jct = 0.11;
        const double var_reversal_potentials__Ca_jct = var_Ca_buffer__Ca_jct;
        const double var_reversal_potentials__E_Ca_jct = ((var_reversal_potentials__R * var_reversal_potentials__T) / (2.0 * var_reversal_potentials__F)) * log(var_reversal_potentials__Cao / var_reversal_potentials__Ca_jct);
        const double var_ICab__E_Ca_jct = var_reversal_potentials__E_Ca_jct;
        const double var_ICab__i_Cab_jct = var_ICab__G_CaBk * var_ICab__Fx_CaBk_jct * (var_ICab__V - var_ICab__E_Ca_jct);
        const double var_ICab__i_Cab = var_ICab__i_Cab_SL + var_ICab__i_Cab_jct;
        const double var_cell__i_Cab = var_ICab__i_Cab;
        const double var_ICap__H = 1.6;
        const double var_ICap__Fx_SLCaP_SL = 0.89;
        const double var_ICap__Km = 0.0005;
        const double var_ICap__V_maxAF = 0.0673;
        const double var_ICap__Q10_SLCaP = 2.35;
        const double var_ICap__T = var_model_parameters__T;
        const double var_ICap__Q_SLCaP = pow(var_ICap__Q10_SLCaP, (var_ICap__T - 310.0) / 10.0);
        const double var_ICap__Ca_SL = var_Ca_buffer__Ca_SL;
        const double var_ICap__i_Cap_SL = (var_ICap__Q_SLCaP * var_ICap__V_maxAF * var_ICap__Fx_SLCaP_SL) / (1.0 + pow(var_ICap__Km / var_ICap__Ca_SL, var_ICap__H));
        const double var_ICap__Fx_SLCaP_jct = 0.11;
        const double var_ICap__Ca_jct = var_Ca_buffer__Ca_jct;
        const double var_ICap__i_Cap_jct = (var_ICap__Q_SLCaP * var_ICap__V_maxAF * var_ICap__Fx_SLCaP_jct) / (1.0 + pow(var_ICap__Km / var_ICap__Ca_jct, var_ICap__H));
        const double var_ICap__i_Cap = var_ICap__i_Cap_jct + var_ICap__i_Cap_SL;
        const double var_cell__i_Cap = var_ICap__i_Cap;
        double var_cell__i_Stim = GetStimulus((1.0/1)*var_environment__time);
        const double var_model_parameters__Mgi = 1.0;
        const double var_model_parameters__cell_length = 100.0;
        const double var_model_parameters__cell_radius = 10.25;
        const double var_model_parameters__Cm_per_area = 2e-06;
        const double var_model_parameters__Cm = (((var_model_parameters__Cm_per_area * 2.0 * var_model_parameters__cell_radius) / 10000.0) * M_PI * var_model_parameters__cell_length) / 10000.0;
        const double var_model_parameters__Vol_Cell = (3.141592654 * pow(var_model_parameters__cell_radius / 1000.0, 2.0) * var_model_parameters__cell_length) / pow(1000.0, 3.0);
        const double var_model_parameters__Vol_SR = 0.035 * var_model_parameters__Vol_Cell;
        const double var_model_parameters__Vol_SL = 0.02 * var_model_parameters__Vol_Cell;
        const double var_model_parameters__Vol_jct = 0.00051 * var_model_parameters__Vol_Cell;
        const double var_model_parameters__Vol_cytosol = 0.65 * var_model_parameters__Vol_Cell;
        const double var_INa_h_gate__V = var_INa__V;
        const double var_INa_h_gate__alpha_h = (var_INa_h_gate__V < (-40.0)) ? (0.135 * exp((80.0 + var_INa_h_gate__V) / (-6.8))) : 0.0;
        const double var_INa_h_gate__beta_h = (var_INa_h_gate__V < (-40.0)) ? ((3.56 * exp(0.079 * var_INa_h_gate__V)) + (310000.0 * exp(0.35 * var_INa_h_gate__V))) : (1.0 / (0.13 * (1.0 + exp((var_INa_h_gate__V + 10.66) / (-11.1)))));
        const double var_INa_h_gate__tau_h = 1.0 / (var_INa_h_gate__alpha_h + var_INa_h_gate__beta_h);
        const double var_INa_h_gate__h_infinity = var_INa_h_gate__alpha_h / (var_INa_h_gate__alpha_h + var_INa_h_gate__beta_h);
        const double var_INa_j_gate__V = var_INa__V;
        const double var_INa_j_gate__alpha_j = (var_INa_j_gate__V < (-40.0)) ? ((((((-127140.0) * exp(0.2444 * var_INa_j_gate__V)) - (3.474e-05 * exp((-0.04391) * var_INa_j_gate__V))) * (var_INa_j_gate__V + 37.78)) / 1.0) / (1.0 + exp(0.311 * (var_INa_j_gate__V + 79.23)))) : 0.0;
        const double var_INa_j_gate__beta_j = (var_INa_j_gate__V < (-40.0)) ? ((0.1212 * exp((-0.01052) * var_INa_j_gate__V)) / (1.0 + exp((-0.1378) * (var_INa_j_gate__V + 40.14)))) : ((0.3 * exp((-2.535e-07) * var_INa_j_gate__V)) / (1.0 + exp((-0.1) * (var_INa_j_gate__V + 32.0))));
        const double var_INa_j_gate__tau_j = 1.0 / (var_INa_j_gate__alpha_j + var_INa_j_gate__beta_j);
        const double var_INa_j_gate__j_infinity = var_INa_j_gate__alpha_j / (var_INa_j_gate__alpha_j + var_INa_j_gate__beta_j);
        const double var_INa_m_gate__V = var_INa__V;
        const double var_INa_m_gate__alpha_m = ((0.32 * (var_INa_m_gate__V + 47.13)) / 1.0) / (1.0 - exp((-0.1) * (var_INa_m_gate__V + 47.13)));
        const double var_INa_m_gate__beta_m = 0.08 * exp((-var_INa_m_gate__V) / 11.0);
        const double var_INa_m_gate__tau_m = 1.0 / (var_INa_m_gate__alpha_m + var_INa_m_gate__beta_m);
        const double var_INa_m_gate__m_infinity = var_INa_m_gate__alpha_m / (var_INa_m_gate__alpha_m + var_INa_m_gate__beta_m);
        const double var_IKr_Xr_gate__V = var_IKr__V;
        const double var_IKr_Xr_gate__Xr_infinity = 1.0 / (1.0 + exp((-(50.0 + var_IKr_Xr_gate__V)) / 7.5));
        const double var_IKr_Xr_gate__tau_Xr = 1.0 / (((0.00138 * (var_IKr_Xr_gate__V + 7.0)) / (1.0 - exp((-0.123) * (var_IKr_Xr_gate__V + 7.0)))) + ((0.00061 * (var_IKr_Xr_gate__V + 10.0)) / (exp(0.145 * (var_IKr_Xr_gate__V + 10.0)) - 1.0)));
        const double var_IKs_Xs_gate__V = var_IKs__V;
        const double var_IKs_Xs_gate__Xs_infinity = 1.0 / (1.0 + exp((-(var_IKs_Xs_gate__V - 1.5)) / 16.7));
        const double var_IKs_Xs_gate__tau_Xs = 1.0 / (((7.19e-05 * (var_IKs_Xs_gate__V + 30.0)) / (1.0 - exp((-0.148) * (var_IKs_Xs_gate__V + 30.0)))) + ((0.000131 * (var_IKs_Xs_gate__V + 30.0)) / ((-1.0) + exp(0.0687 * (var_IKs_Xs_gate__V + 30.0)))));
        const double var_Itos_X_gate__V = var_Itos__V;
        const double var_Itos_X_gate__X_tos_infinity = 1.0 / (1.0 + exp((-(var_Itos_X_gate__V + 3.0)) / 15.0));
        const double var_Itos_X_gate__tau_X_tos = (9.0 / (1.0 + exp((var_Itos_X_gate__V + 3.0) / 15.0))) + 0.5;
        const double var_Itos_Y_gate__V = var_Itos__V;
        const double var_Itos_Y_gate__Y_tos_infinity = 1.0 / (1.0 + exp((var_Itos_Y_gate__V + 33.5) / 10.0));
        const double var_Itos_Y_gate__tau_Y_tos = (3000.0 / (1.0 + exp((var_Itos_Y_gate__V + 60.0) / 10.0))) + 30.0;
        const double var_Itos_R_gate__V = var_Itos__V;
        const double var_Itos_R_gate__R_tos_infinity = 1.0 / (1.0 + exp((var_Itos_R_gate__V + 33.5) / 10.0));
        const double var_Itos_R_gate__tau_R_tos = (2800.0 / (1.0 + exp((var_Itos_R_gate__V + 60.0) / 10.0))) + 220.0;
        const double var_Itof_X_gate__V = var_Itof__V;
        const double var_Itof_X_gate__X_tof_infinity = 1.0 / (1.0 + exp((-(var_Itof_X_gate__V + 3.0)) / 15.0));
        const double var_Itof_X_gate__tau_X_tof = (3.5 * exp(-pow(var_Itof_X_gate__V / 30.0, 2.0))) + 1.5;
        const double var_Itof_Y_gate__V = var_Itof__V;
        const double var_Itof_Y_gate__Y_tof_infinity = 1.0 / (1.0 + exp((var_Itof_Y_gate__V + 33.5) / 10.0));
        const double var_Itof_Y_gate__tau_Y_tof = (20.0 / (1.0 + exp((var_Itof_Y_gate__V + 33.5) / 10.0))) + 20.0;
        const double var_ICaL_d_gate__V = var_ICaL__V;
        const double var_ICaL_d_gate__d_infinity = 1.0 / (1.0 + exp((-(var_ICaL_d_gate__V + 14.5)) / 6.0));
        const double var_ICaL_d_gate__tau_d = (1.0 * var_ICaL_d_gate__d_infinity * (1.0 - exp((-(var_ICaL_d_gate__V + 14.5)) / 6.0))) / (0.035 * (var_ICaL_d_gate__V + 14.5));
        const double var_ICaL_f_gate__V = var_ICaL__V;
        const double var_ICaL_f_gate__f_infinity = (1.0 / (1.0 + exp((var_ICaL_f_gate__V + 35.06) / 3.6))) + (0.6 / (1.0 + exp((50.0 - var_ICaL_f_gate__V) / 20.0)));
        const double var_ICaL_f_gate__tau_f = 1.0 / ((0.0197 * exp(-pow(0.0337 * (var_ICaL_f_gate__V + 14.5), 2.0))) + 0.02);
        const double var_ICaL_fCa_gate__Ca_SL = var_ICaL__Ca_SL;
        const double var_ICaL_fCa_gate__Ca_jct = var_ICaL__Ca_jct;
        const double var_Jrel_SR__ks = 25.0;
        const double var_Jrel_SR__Ca_jct = var_Ca_buffer__Ca_jct;
        const double var_Jrel_SR__Ca_SR = var_Ca_buffer__Ca_SR;
        const double var_Jrel_SR__j_rel_SR = var_Jrel_SR__ks * var_Jrel_SR__O * (var_Jrel_SR__Ca_SR - var_Jrel_SR__Ca_jct);
        const double var_Jrel_SR__Max_SR = 15.0;
        const double var_Jrel_SR__Min_SR = 1.0;
        const double var_Jrel_SR__EC50_SR = 0.45;
        const double var_Jrel_SR__RI = ((1.0 - var_Jrel_SR__R) - var_Jrel_SR__O) - var_Jrel_SR__I;
        const double var_Jrel_SR__koCa = 10.0;
        const double var_Jrel_SR__kom = 0.06;
        const double var_Jrel_SR__kiCa = 0.5;
        const double var_Jrel_SR__kim = 0.005;
        const double var_Jrel_SR__HSR = 2.5;
        const double var_Jrel_SR__kCaSR = var_Jrel_SR__Max_SR - ((var_Jrel_SR__Max_SR - var_Jrel_SR__Min_SR) / (1.0 + pow(var_Jrel_SR__EC50_SR / var_Jrel_SR__Ca_SR, var_Jrel_SR__HSR)));
        const double var_Jrel_SR__koSRCa = var_Jrel_SR__koCa / var_Jrel_SR__kCaSR;
        const double var_Jrel_SR__kiSRCa = var_Jrel_SR__kiCa * var_Jrel_SR__kCaSR;
        const double var_Jleak_SR__Ca_jct = var_Ca_buffer__Ca_jct;
        const double var_Jleak_SR__Ca_SR = var_Ca_buffer__Ca_SR;
        const double var_Jleak_SR__KSRleak = 5.348e-06;
        const double var_Jleak_SR__j_leak_SR = var_Jleak_SR__KSRleak * (var_Jleak_SR__Ca_SR - var_Jleak_SR__Ca_jct);
        const double var_Jpump_SR__Ca_SR = var_Ca_buffer__Ca_SR;
        const double var_Jpump_SR__T = var_model_parameters__T;
        const double var_Jpump_SR__Q10_SRCaP = 2.6;
        const double var_Jpump_SR__Q_SRCaP = pow(var_Jpump_SR__Q10_SRCaP, (var_Jpump_SR__T - 310.0) / 10.0);
        const double var_Jpump_SR__Kmf = 0.000246;
        const double var_Jpump_SR__Vol_cytosol = var_model_parameters__Vol_cytosol;
        const double var_Jpump_SR__Vol_SR = var_model_parameters__Vol_SR;
        const double var_Jpump_SR__H = 1.787;
        const double var_Jpump_SR__Cai = var_Ca_buffer__Cai;
        const double var_Jpump_SR__Kmr = 1.7;
        const double var_Jpump_SR__V_max = 0.000286;
        const double var_Jpump_SR__j_pump_SR = (((var_Jpump_SR__Q_SRCaP * var_Jpump_SR__V_max * var_Jpump_SR__Vol_cytosol) / var_Jpump_SR__Vol_SR) * (pow(var_Jpump_SR__Cai / var_Jpump_SR__Kmf, var_Jpump_SR__H) - pow(var_Jpump_SR__Ca_SR / var_Jpump_SR__Kmr, var_Jpump_SR__H))) / (1.0 + pow(var_Jpump_SR__Cai / var_Jpump_SR__Kmf, var_Jpump_SR__H) + pow(var_Jpump_SR__Ca_SR / var_Jpump_SR__Kmr, var_Jpump_SR__H));
        const double var_ion_diffusion__Na_jct = var_Na_buffer__Na_jct;
        const double var_ion_diffusion__Na_SL = var_Na_buffer__Na_SL;
        const double var_ion_diffusion__J_Na_jct_SL = (var_ion_diffusion__Na_jct - var_ion_diffusion__Na_SL) * 1.8313e-14;
        const double var_ion_diffusion__Nai = var_Na_buffer__Nai;
        const double var_ion_diffusion__J_Na_SL_cytosol = (var_ion_diffusion__Na_SL - var_ion_diffusion__Nai) * 1.6386e-12;
        const double var_ion_diffusion__Ca_SL = var_Ca_buffer__Ca_SL;
        const double var_ion_diffusion__Ca_jct = var_Ca_buffer__Ca_jct;
        const double var_ion_diffusion__J_Ca_jct_SL = (var_ion_diffusion__Ca_jct - var_ion_diffusion__Ca_SL) * 8.2413e-13;
        const double var_ion_diffusion__Cai = var_Ca_buffer__Cai;
        const double var_ion_diffusion__J_Ca_SL_cytosol = (var_ion_diffusion__Ca_SL - var_ion_diffusion__Cai) * 3.7243e-12;
        const double var_Na_buffer__Bmax_SL = 1.65;
        const double var_Na_buffer__Bmax_jct = 7.561;
        const double var_Na_buffer__kon = 0.0001;
        const double var_Na_buffer__koff = 0.001;
        const double var_Na_buffer__J_Na_jct_SL = var_ion_diffusion__J_Na_jct_SL;
        const double var_Na_buffer__J_Na_SL_cytosol = var_ion_diffusion__J_Na_SL_cytosol;
        const double var_Na_buffer__i_Na_jct = var_INa__i_Na_jct;
        const double var_Na_buffer__i_NaCa_jct = var_INaCa__i_NaCa_jct;
        const double var_Na_buffer__i_Nab_jct = var_INab__i_Nab_jct;
        const double var_Na_buffer__i_NaK_jct = var_INaK__i_NaK_jct;
        const double var_Na_buffer__i_CaL_Na_jct = var_ICaL__i_CaL_Na_jct;
        const double var_Na_buffer__i_Na_SL = var_INa__i_Na_SL;
        const double var_Na_buffer__i_NaCa_SL = var_INaCa__i_NaCa_SL;
        const double var_Na_buffer__i_Nab_SL = var_INab__i_Nab_SL;
        const double var_Na_buffer__i_NaK_SL = var_INaK__i_NaK_SL;
        const double var_Na_buffer__i_CaL_Na_SL = var_ICaL__i_CaL_Na_SL;
        const double var_Na_buffer__Vol_SL = var_model_parameters__Vol_SL;
        const double var_Na_buffer__Vol_jct = var_model_parameters__Vol_jct;
        const double var_Na_buffer__Vol_cytosol = var_model_parameters__Vol_cytosol;
        const double var_Na_buffer__F = var_model_parameters__F;
        const double var_Na_buffer__Cm = var_model_parameters__Cm;
        const double var_Na_buffer__dNa_jct_buf = (var_Na_buffer__kon * var_Na_buffer__Na_jct * (var_Na_buffer__Bmax_jct - var_Na_buffer__Na_jct_buf)) - (var_Na_buffer__koff * var_Na_buffer__Na_jct_buf);
        const double var_Na_buffer__dNa_SL_buf = (var_Na_buffer__kon * var_Na_buffer__Na_SL * (var_Na_buffer__Bmax_SL - var_Na_buffer__Na_SL_buf)) - (var_Na_buffer__koff * var_Na_buffer__Na_SL_buf);
        const double var_Ca_buffer__Mgi = var_model_parameters__Mgi;
        const double var_Ca_buffer__Bmax_SLB_SL = 0.0374;
        const double var_Ca_buffer__Bmax_SLB_jct = 0.0046;
        const double var_Ca_buffer__Bmax_SLHigh_SL = 0.0134;
        const double var_Ca_buffer__Bmax_SLHigh_jct = 0.00165;
        const double var_Ca_buffer__Bmax_Calsequestrin = 0.14;
        const double var_Ca_buffer__kon_SL = 100.0;
        const double var_Ca_buffer__kon_Calsequestrin = 100.0;
        const double var_Ca_buffer__koff_SLB = 1.3;
        const double var_Ca_buffer__koff_SLHigh = 0.03;
        const double var_Ca_buffer__koff_Calsequestrin = 65.0;
        const double var_Ca_buffer__i_CaL_Ca_jct = var_ICaL__i_CaL_Ca_jct;
        const double var_Ca_buffer__i_NaCa_jct = var_INaCa__i_NaCa_jct;
        const double var_Ca_buffer__i_Cab_jct = var_ICab__i_Cab_jct;
        const double var_Ca_buffer__i_Cap_jct = var_ICap__i_Cap_jct;
        const double var_Ca_buffer__i_CaL_Ca_SL = var_ICaL__i_CaL_Ca_SL;
        const double var_Ca_buffer__i_NaCa_SL = var_INaCa__i_NaCa_SL;
        const double var_Ca_buffer__i_Cab_SL = var_ICab__i_Cab_SL;
        const double var_Ca_buffer__i_Cap_SL = var_ICap__i_Cap_SL;
        const double var_Ca_buffer__j_pump_SR = var_Jpump_SR__j_pump_SR;
        const double var_Ca_buffer__j_rel_SR = var_Jrel_SR__j_rel_SR;
        const double var_Ca_buffer__j_leak_SR = var_Jleak_SR__j_leak_SR;
        const double var_Ca_buffer__J_Ca_jct_SL = var_ion_diffusion__J_Ca_jct_SL;
        const double var_Ca_buffer__J_Ca_SL_cytosol = var_ion_diffusion__J_Ca_SL_cytosol;
        const double var_Ca_buffer__Vol_SR = var_model_parameters__Vol_SR;
        const double var_Ca_buffer__Vol_SL = var_model_parameters__Vol_SL;
        const double var_Ca_buffer__Vol_jct = var_model_parameters__Vol_jct;
        const double var_Ca_buffer__Vol_cytosol = var_model_parameters__Vol_cytosol;
        const double var_Ca_buffer__F = var_model_parameters__F;
        const double var_Ca_buffer__Cm = var_model_parameters__Cm;
        const double var_Ca_buffer__dCalsequestrin = (var_Ca_buffer__kon_Calsequestrin * var_Ca_buffer__Ca_SR * (((var_Ca_buffer__Bmax_Calsequestrin * var_Ca_buffer__Vol_cytosol) / var_Ca_buffer__Vol_SR) - var_Ca_buffer__Ca_Calsequestrin)) - (var_Ca_buffer__koff_Calsequestrin * var_Ca_buffer__Ca_Calsequestrin);
        const double var_cytosolic_Ca_buffer__Bmax_Myosin_Mg = 0.14;
        const double var_cytosolic_Ca_buffer__koff_Myosin_Mg = 5.7e-05;
        const double var_cytosolic_Ca_buffer__Mgi = var_Ca_buffer__Mgi;
        const double var_cytosolic_Ca_buffer__kon_Myosin_Mg = 0.0157;
        const double var_cytosolic_Ca_buffer__dMg_Myosin = (var_cytosolic_Ca_buffer__kon_Myosin_Mg * var_cytosolic_Ca_buffer__Mgi * (var_cytosolic_Ca_buffer__Bmax_Myosin_Mg - (var_cytosolic_Ca_buffer__Ca_Myosin + var_cytosolic_Ca_buffer__Mg_Myosin))) - (var_cytosolic_Ca_buffer__koff_Myosin_Mg * var_cytosolic_Ca_buffer__Mg_Myosin);
        const double var_cytosolic_Ca_buffer__Bmax_TroponinC_Ca_Mg_Mg = 0.14;
        const double var_cytosolic_Ca_buffer__kon_TroponinC_Ca_Mg_Mg = 0.003;
        const double var_cytosolic_Ca_buffer__koff_TroponinC_Ca_Mg_Mg = 0.00333;
        const double var_cytosolic_Ca_buffer__dMg_TroponinC_Ca_Mg = (var_cytosolic_Ca_buffer__kon_TroponinC_Ca_Mg_Mg * var_cytosolic_Ca_buffer__Mgi * (var_cytosolic_Ca_buffer__Bmax_TroponinC_Ca_Mg_Mg - (var_cytosolic_Ca_buffer__Ca_TroponinC_Ca_Mg + var_cytosolic_Ca_buffer__Mg_TroponinC_Ca_Mg))) - (var_cytosolic_Ca_buffer__koff_TroponinC_Ca_Mg_Mg * var_cytosolic_Ca_buffer__Mg_TroponinC_Ca_Mg);
        const double var_cytosolic_Ca_buffer__Cai = var_Ca_buffer__Cai;
        const double var_cytosolic_Ca_buffer__kon_TroponinC_Ca_Mg_Ca = 2.37;
        const double var_cytosolic_Ca_buffer__koff_TroponinC_Ca_Mg_Ca = 3.2e-05;
        const double var_cytosolic_Ca_buffer__Bmax_TroponinC_Ca_Mg_Ca = 0.14;
        const double var_cytosolic_Ca_buffer__dCa_TroponinC_Ca_Mg = (var_cytosolic_Ca_buffer__kon_TroponinC_Ca_Mg_Ca * var_cytosolic_Ca_buffer__Cai * (var_cytosolic_Ca_buffer__Bmax_TroponinC_Ca_Mg_Ca - (var_cytosolic_Ca_buffer__Ca_TroponinC_Ca_Mg + var_cytosolic_Ca_buffer__Mg_TroponinC_Ca_Mg))) - (var_cytosolic_Ca_buffer__koff_TroponinC_Ca_Mg_Ca * var_cytosolic_Ca_buffer__Ca_TroponinC_Ca_Mg);
        const double var_cytosolic_Ca_buffer__kon_Myosin_Ca = 13.8;
        const double var_cytosolic_Ca_buffer__Bmax_Myosin_Ca = 0.14;
        const double var_cytosolic_Ca_buffer__koff_Myosin_Ca = 0.00046;
        const double var_cytosolic_Ca_buffer__dCa_Myosin = (var_cytosolic_Ca_buffer__kon_Myosin_Ca * var_cytosolic_Ca_buffer__Cai * (var_cytosolic_Ca_buffer__Bmax_Myosin_Ca - (var_cytosolic_Ca_buffer__Ca_Myosin + var_cytosolic_Ca_buffer__Mg_Myosin))) - (var_cytosolic_Ca_buffer__koff_Myosin_Ca * var_cytosolic_Ca_buffer__Ca_Myosin);
        const double var_cytosolic_Ca_buffer__koff_Calmodulin = 0.238;
        const double var_cytosolic_Ca_buffer__kon_Calmodulin = 34.0;
        const double var_cytosolic_Ca_buffer__Bmax_Calmodulin = 0.024;
        const double var_cytosolic_Ca_buffer__dCa_Calmodulin = (var_cytosolic_Ca_buffer__kon_Calmodulin * var_cytosolic_Ca_buffer__Cai * (var_cytosolic_Ca_buffer__Bmax_Calmodulin - var_cytosolic_Ca_buffer__Ca_Calmodulin)) - (var_cytosolic_Ca_buffer__koff_Calmodulin * var_cytosolic_Ca_buffer__Ca_Calmodulin);
        const double var_cytosolic_Ca_buffer__Bmax_TroponinC = 0.07;
        const double var_cytosolic_Ca_buffer__kon_TroponinC = 32.7;
        const double var_cytosolic_Ca_buffer__koff_TroponinC = 0.0196;
        const double var_cytosolic_Ca_buffer__dCa_TroponinC = (var_cytosolic_Ca_buffer__kon_TroponinC * var_cytosolic_Ca_buffer__Cai * (var_cytosolic_Ca_buffer__Bmax_TroponinC - var_cytosolic_Ca_buffer__Ca_TroponinC)) - (var_cytosolic_Ca_buffer__koff_TroponinC * var_cytosolic_Ca_buffer__Ca_TroponinC);
        const double var_cytosolic_Ca_buffer__koff_SRB = 0.06;
        const double var_cytosolic_Ca_buffer__Bmax_SRB = 0.0171;
        const double var_cytosolic_Ca_buffer__kon_SRB = 100.0;
        const double var_cytosolic_Ca_buffer__dCa_SRB = (var_cytosolic_Ca_buffer__kon_SRB * var_cytosolic_Ca_buffer__Cai * (var_cytosolic_Ca_buffer__Bmax_SRB - var_cytosolic_Ca_buffer__Ca_SRB)) - (var_cytosolic_Ca_buffer__koff_SRB * var_cytosolic_Ca_buffer__Ca_SRB);
        const double var_cytosolic_Ca_buffer__dCa_cytosol_tot_bound = var_cytosolic_Ca_buffer__dCa_TroponinC + var_cytosolic_Ca_buffer__dCa_TroponinC_Ca_Mg + var_cytosolic_Ca_buffer__dMg_TroponinC_Ca_Mg + var_cytosolic_Ca_buffer__dCa_Calmodulin + var_cytosolic_Ca_buffer__dCa_Myosin + var_cytosolic_Ca_buffer__dMg_Myosin + var_cytosolic_Ca_buffer__dCa_SRB;
        const double var_Ca_buffer__dCa_cytosol_tot_bound = var_cytosolic_Ca_buffer__dCa_cytosol_tot_bound;
        const double var_Ca_buffer__dCa_SLB_SL = (var_Ca_buffer__kon_SL * var_Ca_buffer__Ca_SL * (((var_Ca_buffer__Bmax_SLB_SL * var_Ca_buffer__Vol_cytosol) / var_Ca_buffer__Vol_SL) - var_Ca_buffer__Ca_SLB_SL)) - (var_Ca_buffer__koff_SLB * var_Ca_buffer__Ca_SLB_SL);
        const double var_Ca_buffer__dCa_SLB_jct = (var_Ca_buffer__kon_SL * var_Ca_buffer__Ca_jct * (((var_Ca_buffer__Bmax_SLB_jct * 0.1 * var_Ca_buffer__Vol_cytosol) / var_Ca_buffer__Vol_jct) - var_Ca_buffer__Ca_SLB_jct)) - (var_Ca_buffer__koff_SLB * var_Ca_buffer__Ca_SLB_jct);
        const double var_Ca_buffer__dCa_SLHigh_SL = (var_Ca_buffer__kon_SL * var_Ca_buffer__Ca_SL * (((var_Ca_buffer__Bmax_SLHigh_SL * var_Ca_buffer__Vol_cytosol) / var_Ca_buffer__Vol_SL) - var_Ca_buffer__Ca_SLHigh_SL)) - (var_Ca_buffer__koff_SLHigh * var_Ca_buffer__Ca_SLHigh_SL);
        const double var_Ca_buffer__dCa_SLHigh_jct = (var_Ca_buffer__kon_SL * var_Ca_buffer__Ca_jct * (((var_Ca_buffer__Bmax_SLHigh_jct * 0.1 * var_Ca_buffer__Vol_cytosol) / var_Ca_buffer__Vol_jct) - var_Ca_buffer__Ca_SLHigh_jct)) - (var_Ca_buffer__koff_SLHigh * var_Ca_buffer__Ca_SLHigh_jct);
        const double var_Ca_buffer__dCa_jct_tot_bound = var_Ca_buffer__dCa_SLB_jct + var_Ca_buffer__dCa_SLHigh_jct;
        const double var_Ca_buffer__dCa_SL_tot_bound = var_Ca_buffer__dCa_SLB_SL + var_Ca_buffer__dCa_SLHigh_SL;
        const double var_Ca_buffer__i_Ca_jct_tot = (var_Ca_buffer__i_CaL_Ca_jct - (2.0 * var_Ca_buffer__i_NaCa_jct)) + var_Ca_buffer__i_Cab_jct + var_Ca_buffer__i_Cap_jct;
        const double var_Ca_buffer__i_Ca_SL_tot = (var_Ca_buffer__i_CaL_Ca_SL - (2.0 * var_Ca_buffer__i_NaCa_SL)) + var_Ca_buffer__i_Cab_SL + var_Ca_buffer__i_Cap_SL;
        const double var_indo_fluo_Ca_buffer_not_connected__Cai = var_Ca_buffer__Cai;
        const double var_indo_fluo_Ca_buffer_not_connected__Ca_SL = var_Ca_buffer__Ca_SL;
        const double var_indo_fluo_Ca_buffer_not_connected__Ca_jct = var_Ca_buffer__Ca_jct;
        const double var_indo_fluo_Ca_buffer_not_connected__Vol_SL = var_Ca_buffer__Vol_SL;
        const double var_indo_fluo_Ca_buffer_not_connected__Vol_jct = var_Ca_buffer__Vol_jct;
        const double var_indo_fluo_Ca_buffer_not_connected__Vol_cytosol = var_Ca_buffer__Vol_cytosol;
        const double var_indo_fluo_Ca_buffer_not_connected__Indo1 = 0.0;
        const double var_indo_fluo_Ca_buffer_not_connected__Fluo3 = 0.0;
        const double var_indo_fluo_Ca_buffer_not_connected__Bmax_Indo1_Cytosol = 0.025;
        const double var_indo_fluo_Ca_buffer_not_connected__Bmax_Indo1_SL = 0.00077;
        const double var_indo_fluo_Ca_buffer_not_connected__Bmax_Indo1_jct = 2e-05;
        const double var_indo_fluo_Ca_buffer_not_connected__Bmax_Fluo3_Cytosol = 0.025;
        const double var_indo_fluo_Ca_buffer_not_connected__Bmax_Fluo3_SL = 0.00077;
        const double var_indo_fluo_Ca_buffer_not_connected__Bmax_Fluo3_jct = 2e-05;
        const double var_indo_fluo_Ca_buffer_not_connected__kon_Indo1 = 100.0;
        const double var_indo_fluo_Ca_buffer_not_connected__kon_Fluo3 = 100.0;
        const double var_indo_fluo_Ca_buffer_not_connected__koff_Indo1 = 0.06;
        const double var_indo_fluo_Ca_buffer_not_connected__koff_Fluo3 = 0.11;
        const double var_indo_fluo_Ca_buffer_not_connected__Indo1Bound = var_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_Cytosol + ((var_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_jct * var_indo_fluo_Ca_buffer_not_connected__Vol_jct) / var_indo_fluo_Ca_buffer_not_connected__Vol_cytosol) + ((var_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_SL * var_indo_fluo_Ca_buffer_not_connected__Vol_SL) / var_indo_fluo_Ca_buffer_not_connected__Vol_cytosol);
        const double var_indo_fluo_Ca_buffer_not_connected__Fluo3Bound = var_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_Cytosol + ((var_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_jct * var_indo_fluo_Ca_buffer_not_connected__Vol_jct) / var_indo_fluo_Ca_buffer_not_connected__Vol_cytosol) + ((var_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_SL * var_indo_fluo_Ca_buffer_not_connected__Vol_SL) / var_indo_fluo_Ca_buffer_not_connected__Vol_cytosol);
        const double var_indo_fluo_Ca_buffer_not_connected__dCa_Indo1_Cytosol = (var_indo_fluo_Ca_buffer_not_connected__kon_Indo1 * (var_indo_fluo_Ca_buffer_not_connected__Indo1 - var_indo_fluo_Ca_buffer_not_connected__Indo1Bound) * var_indo_fluo_Ca_buffer_not_connected__Cai * (var_indo_fluo_Ca_buffer_not_connected__Bmax_Indo1_Cytosol - var_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_Cytosol)) - (var_indo_fluo_Ca_buffer_not_connected__koff_Indo1 * var_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_Cytosol);
        const double var_indo_fluo_Ca_buffer_not_connected__dCa_Indo1_jct = (var_indo_fluo_Ca_buffer_not_connected__kon_Indo1 * (var_indo_fluo_Ca_buffer_not_connected__Indo1 - var_indo_fluo_Ca_buffer_not_connected__Indo1Bound) * var_indo_fluo_Ca_buffer_not_connected__Ca_jct * (((var_indo_fluo_Ca_buffer_not_connected__Bmax_Indo1_jct * var_indo_fluo_Ca_buffer_not_connected__Vol_cytosol) / var_indo_fluo_Ca_buffer_not_connected__Vol_jct) - var_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_jct)) - (var_indo_fluo_Ca_buffer_not_connected__koff_Indo1 * var_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_jct);
        const double var_indo_fluo_Ca_buffer_not_connected__dCa_Indo1_SL = (var_indo_fluo_Ca_buffer_not_connected__kon_Indo1 * (var_indo_fluo_Ca_buffer_not_connected__Indo1 - var_indo_fluo_Ca_buffer_not_connected__Indo1Bound) * var_indo_fluo_Ca_buffer_not_connected__Ca_SL * (((var_indo_fluo_Ca_buffer_not_connected__Bmax_Indo1_SL * var_indo_fluo_Ca_buffer_not_connected__Vol_cytosol) / var_indo_fluo_Ca_buffer_not_connected__Vol_SL) - var_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_SL)) - (var_indo_fluo_Ca_buffer_not_connected__koff_Indo1 * var_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_SL);
        const double var_indo_fluo_Ca_buffer_not_connected__dCa_Fluo3_Cytosol = (var_indo_fluo_Ca_buffer_not_connected__kon_Fluo3 * (var_indo_fluo_Ca_buffer_not_connected__Fluo3 - var_indo_fluo_Ca_buffer_not_connected__Fluo3Bound) * var_indo_fluo_Ca_buffer_not_connected__Cai * (var_indo_fluo_Ca_buffer_not_connected__Bmax_Fluo3_Cytosol - var_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_Cytosol)) - (var_indo_fluo_Ca_buffer_not_connected__koff_Fluo3 * var_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_Cytosol);
        const double var_indo_fluo_Ca_buffer_not_connected__dCa_Fluo3_jct = (var_indo_fluo_Ca_buffer_not_connected__kon_Fluo3 * (var_indo_fluo_Ca_buffer_not_connected__Fluo3 - var_indo_fluo_Ca_buffer_not_connected__Fluo3Bound) * var_indo_fluo_Ca_buffer_not_connected__Ca_jct * (((var_indo_fluo_Ca_buffer_not_connected__Bmax_Fluo3_jct * var_indo_fluo_Ca_buffer_not_connected__Vol_cytosol) / var_indo_fluo_Ca_buffer_not_connected__Vol_jct) - var_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_jct)) - (var_indo_fluo_Ca_buffer_not_connected__koff_Fluo3 * var_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_jct);
        const double var_indo_fluo_Ca_buffer_not_connected__dCa_Fluo3_SL = (var_indo_fluo_Ca_buffer_not_connected__kon_Fluo3 * (var_indo_fluo_Ca_buffer_not_connected__Fluo3 - var_indo_fluo_Ca_buffer_not_connected__Fluo3Bound) * var_indo_fluo_Ca_buffer_not_connected__Ca_SL * (((var_indo_fluo_Ca_buffer_not_connected__Bmax_Fluo3_SL * var_indo_fluo_Ca_buffer_not_connected__Vol_cytosol) / var_indo_fluo_Ca_buffer_not_connected__Vol_SL) - var_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_SL)) - (var_indo_fluo_Ca_buffer_not_connected__koff_Fluo3 * var_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_SL);

        double d_dt_cell__V;
        if (mSetVoltageDerivativeToZero)
        {
            d_dt_cell__V = 0.0;
        }
        else
        {
            d_dt_cell__V = -(var_cell__i_Na + var_cell__i_Nab + var_cell__i_NaK + var_cell__i_Kr + var_cell__i_Ks + var_cell__i_tos + var_cell__i_tof + var_cell__i_K1 + var_cell__i_NaCa + var_cell__i_Cl_Ca + var_cell__i_Clb + var_cell__i_CaL + var_cell__i_Cab + var_cell__i_Cap + var_cell__i_Stim);
        }

        const double d_dt_INa_h_gate__h = (var_INa_h_gate__h_infinity - var_INa_h_gate__h) / var_INa_h_gate__tau_h;
        const double d_dt_INa_j_gate__j = (var_INa_j_gate__j_infinity - var_INa_j_gate__j) / var_INa_j_gate__tau_j;
        const double d_dt_INa_m_gate__m = (var_INa_m_gate__m_infinity - var_INa_m_gate__m) / var_INa_m_gate__tau_m;
        const double d_dt_IKr_Xr_gate__Xr = (var_IKr_Xr_gate__Xr_infinity - var_IKr_Xr_gate__Xr) / var_IKr_Xr_gate__tau_Xr;
        const double d_dt_IKs_Xs_gate__Xs = (var_IKs_Xs_gate__Xs_infinity - var_IKs_Xs_gate__Xs) / var_IKs_Xs_gate__tau_Xs;
        const double d_dt_Itos_X_gate__X_tos = (var_Itos_X_gate__X_tos_infinity - var_Itos_X_gate__X_tos) / var_Itos_X_gate__tau_X_tos;
        const double d_dt_Itos_Y_gate__Y_tos = (var_Itos_Y_gate__Y_tos_infinity - var_Itos_Y_gate__Y_tos) / var_Itos_Y_gate__tau_Y_tos;
        const double d_dt_Itos_R_gate__R_tos = (var_Itos_R_gate__R_tos_infinity - var_Itos_R_gate__R_tos) / var_Itos_R_gate__tau_R_tos;
        const double d_dt_Itof_X_gate__X_tof = (var_Itof_X_gate__X_tof_infinity - var_Itof_X_gate__X_tof) / var_Itof_X_gate__tau_X_tof;
        const double d_dt_Itof_Y_gate__Y_tof = (var_Itof_Y_gate__Y_tof_infinity - var_Itof_Y_gate__Y_tof) / var_Itof_Y_gate__tau_Y_tof;
        const double d_dt_ICaL_d_gate__d = (var_ICaL_d_gate__d_infinity - var_ICaL_d_gate__d) / var_ICaL_d_gate__tau_d;
        const double d_dt_ICaL_f_gate__f = (var_ICaL_f_gate__f_infinity - var_ICaL_f_gate__f) / var_ICaL_f_gate__tau_f;
        const double d_dt_ICaL_fCa_gate__fCaB_SL = (1.7 * var_ICaL_fCa_gate__Ca_SL * (1.0 - var_ICaL_fCa_gate__fCaB_SL)) - (0.0119 * var_ICaL_fCa_gate__fCaB_SL);
        const double d_dt_ICaL_fCa_gate__fCaB_jct = (1.7 * var_ICaL_fCa_gate__Ca_jct * (1.0 - var_ICaL_fCa_gate__fCaB_jct)) - (0.0119 * var_ICaL_fCa_gate__fCaB_jct);
        const double d_dt_Jrel_SR__R = ((var_Jrel_SR__kim * var_Jrel_SR__RI) - (var_Jrel_SR__kiSRCa * var_Jrel_SR__Ca_jct * var_Jrel_SR__R)) - ((var_Jrel_SR__koSRCa * pow(var_Jrel_SR__Ca_jct, 2.0) * var_Jrel_SR__R) - (var_Jrel_SR__kom * var_Jrel_SR__O));
        const double d_dt_Jrel_SR__O = ((var_Jrel_SR__koSRCa * pow(var_Jrel_SR__Ca_jct, 2.0) * var_Jrel_SR__R) - (var_Jrel_SR__kom * var_Jrel_SR__O)) - ((var_Jrel_SR__kiSRCa * var_Jrel_SR__Ca_jct * var_Jrel_SR__O) - (var_Jrel_SR__kim * var_Jrel_SR__I));
        const double d_dt_Jrel_SR__I = ((var_Jrel_SR__kiSRCa * var_Jrel_SR__Ca_jct * var_Jrel_SR__O) - (var_Jrel_SR__kim * var_Jrel_SR__I)) - ((var_Jrel_SR__kom * var_Jrel_SR__I) - (var_Jrel_SR__koSRCa * pow(var_Jrel_SR__Ca_jct, 2.0) * var_Jrel_SR__RI));
        const double d_dt_Na_buffer__Na_jct_buf = var_Na_buffer__dNa_jct_buf;
        const double d_dt_Na_buffer__Na_SL_buf = var_Na_buffer__dNa_SL_buf;
        const double d_dt_Na_buffer__Na_jct = ((((-var_Na_buffer__Cm) * (var_Na_buffer__i_Na_jct + (3.0 * var_Na_buffer__i_NaCa_jct) + var_Na_buffer__i_Nab_jct + (3.0 * var_Na_buffer__i_NaK_jct) + var_Na_buffer__i_CaL_Na_jct)) / (var_Na_buffer__Vol_jct * var_Na_buffer__F)) - (var_Na_buffer__J_Na_jct_SL / var_Na_buffer__Vol_jct)) - var_Na_buffer__dNa_jct_buf;
        const double d_dt_Na_buffer__Na_SL = ((((-var_Na_buffer__Cm) * (var_Na_buffer__i_Na_SL + (3.0 * var_Na_buffer__i_NaCa_SL) + var_Na_buffer__i_Nab_SL + (3.0 * var_Na_buffer__i_NaK_SL) + var_Na_buffer__i_CaL_Na_SL)) / (var_Na_buffer__Vol_SL * var_Na_buffer__F)) + ((var_Na_buffer__J_Na_jct_SL - var_Na_buffer__J_Na_SL_cytosol) / var_Na_buffer__Vol_SL)) - var_Na_buffer__dNa_SL_buf;
        const double d_dt_Na_buffer__Nai = var_Na_buffer__J_Na_SL_cytosol / var_Na_buffer__Vol_cytosol;
        const double d_dt_Ca_buffer__Ca_Calsequestrin = var_Ca_buffer__dCalsequestrin;
        const double d_dt_Ca_buffer__Ca_SLB_SL = var_Ca_buffer__dCa_SLB_SL;
        const double d_dt_Ca_buffer__Ca_SLB_jct = var_Ca_buffer__dCa_SLB_jct;
        const double d_dt_Ca_buffer__Ca_SLHigh_SL = var_Ca_buffer__dCa_SLHigh_SL;
        const double d_dt_Ca_buffer__Ca_SLHigh_jct = var_Ca_buffer__dCa_SLHigh_jct;
        const double d_dt_Ca_buffer__Ca_SR = (var_Ca_buffer__j_pump_SR - (((var_Ca_buffer__j_leak_SR * var_Ca_buffer__Vol_cytosol) / var_Ca_buffer__Vol_SR) + var_Ca_buffer__j_rel_SR)) - var_Ca_buffer__dCalsequestrin;
        const double d_dt_Ca_buffer__Ca_jct = (((((-var_Ca_buffer__i_Ca_jct_tot) * var_Ca_buffer__Cm) / (var_Ca_buffer__Vol_jct * 2.0 * var_Ca_buffer__F)) - (var_Ca_buffer__J_Ca_jct_SL / var_Ca_buffer__Vol_jct)) + ((var_Ca_buffer__j_rel_SR * var_Ca_buffer__Vol_SR) / var_Ca_buffer__Vol_jct) + ((var_Ca_buffer__j_leak_SR * var_Ca_buffer__Vol_cytosol) / var_Ca_buffer__Vol_jct)) - (1.0 * var_Ca_buffer__dCa_jct_tot_bound);
        const double d_dt_Ca_buffer__Ca_SL = ((((-var_Ca_buffer__i_Ca_SL_tot) * var_Ca_buffer__Cm) / (var_Ca_buffer__Vol_SL * 2.0 * var_Ca_buffer__F)) + ((var_Ca_buffer__J_Ca_jct_SL - var_Ca_buffer__J_Ca_SL_cytosol) / var_Ca_buffer__Vol_SL)) - (1.0 * var_Ca_buffer__dCa_SL_tot_bound);
        const double d_dt_Ca_buffer__Cai = ((((-var_Ca_buffer__j_pump_SR) * var_Ca_buffer__Vol_SR) / var_Ca_buffer__Vol_cytosol) + (var_Ca_buffer__J_Ca_SL_cytosol / var_Ca_buffer__Vol_cytosol)) - (1.0 * var_Ca_buffer__dCa_cytosol_tot_bound);
        const double d_dt_cytosolic_Ca_buffer__Ca_TroponinC = var_cytosolic_Ca_buffer__dCa_TroponinC;
        const double d_dt_cytosolic_Ca_buffer__Ca_TroponinC_Ca_Mg = var_cytosolic_Ca_buffer__dCa_TroponinC_Ca_Mg;
        const double d_dt_cytosolic_Ca_buffer__Mg_TroponinC_Ca_Mg = var_cytosolic_Ca_buffer__dMg_TroponinC_Ca_Mg;
        const double d_dt_cytosolic_Ca_buffer__Ca_Calmodulin = var_cytosolic_Ca_buffer__dCa_Calmodulin;
        const double d_dt_cytosolic_Ca_buffer__Ca_Myosin = var_cytosolic_Ca_buffer__dCa_Myosin;
        const double d_dt_cytosolic_Ca_buffer__Mg_Myosin = var_cytosolic_Ca_buffer__dMg_Myosin;
        const double d_dt_cytosolic_Ca_buffer__Ca_SRB = var_cytosolic_Ca_buffer__dCa_SRB;
        const double d_dt_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_Cytosol = var_indo_fluo_Ca_buffer_not_connected__dCa_Indo1_Cytosol;
        const double d_dt_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_SL = var_indo_fluo_Ca_buffer_not_connected__dCa_Indo1_SL;
        const double d_dt_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_jct = var_indo_fluo_Ca_buffer_not_connected__dCa_Indo1_jct;
        const double d_dt_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_Cytosol = var_indo_fluo_Ca_buffer_not_connected__dCa_Fluo3_Cytosol;
        const double d_dt_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_SL = var_indo_fluo_Ca_buffer_not_connected__dCa_Fluo3_SL;
        const double d_dt_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_jct = var_indo_fluo_Ca_buffer_not_connected__dCa_Fluo3_jct;

        NV_Ith_S(ydot, 0) = d_dt_cell__V;
        NV_Ith_S(ydot, 1) = d_dt_INa_h_gate__h;
        NV_Ith_S(ydot, 2) = d_dt_INa_j_gate__j;
        NV_Ith_S(ydot, 3) = d_dt_INa_m_gate__m;
        NV_Ith_S(ydot, 4) = d_dt_IKr_Xr_gate__Xr;
        NV_Ith_S(ydot, 5) = d_dt_IKs_Xs_gate__Xs;
        NV_Ith_S(ydot, 6) = d_dt_Itos_X_gate__X_tos;
        NV_Ith_S(ydot, 7) = d_dt_Itos_Y_gate__Y_tos;
        NV_Ith_S(ydot, 8) = d_dt_Itos_R_gate__R_tos;
        NV_Ith_S(ydot, 9) = d_dt_Itof_X_gate__X_tof;
        NV_Ith_S(ydot, 10) = d_dt_Itof_Y_gate__Y_tof;
        NV_Ith_S(ydot, 11) = d_dt_ICaL_d_gate__d;
        NV_Ith_S(ydot, 12) = d_dt_ICaL_f_gate__f;
        NV_Ith_S(ydot, 13) = d_dt_ICaL_fCa_gate__fCaB_SL;
        NV_Ith_S(ydot, 14) = d_dt_ICaL_fCa_gate__fCaB_jct;
        NV_Ith_S(ydot, 15) = d_dt_Jrel_SR__R;
        NV_Ith_S(ydot, 16) = d_dt_Jrel_SR__I;
        NV_Ith_S(ydot, 17) = d_dt_Jrel_SR__O;
        NV_Ith_S(ydot, 18) = d_dt_Na_buffer__Na_SL;
        NV_Ith_S(ydot, 19) = d_dt_Na_buffer__Na_jct;
        NV_Ith_S(ydot, 20) = d_dt_Na_buffer__Na_SL_buf;
        NV_Ith_S(ydot, 21) = d_dt_Na_buffer__Na_jct_buf;
        NV_Ith_S(ydot, 22) = d_dt_Na_buffer__Nai;
        NV_Ith_S(ydot, 23) = d_dt_Ca_buffer__Ca_SR;
        NV_Ith_S(ydot, 24) = d_dt_Ca_buffer__Ca_SL;
        NV_Ith_S(ydot, 25) = d_dt_Ca_buffer__Ca_jct;
        NV_Ith_S(ydot, 26) = d_dt_Ca_buffer__Cai;
        NV_Ith_S(ydot, 27) = d_dt_Ca_buffer__Ca_SLB_SL;
        NV_Ith_S(ydot, 28) = d_dt_Ca_buffer__Ca_SLB_jct;
        NV_Ith_S(ydot, 29) = d_dt_Ca_buffer__Ca_SLHigh_SL;
        NV_Ith_S(ydot, 30) = d_dt_Ca_buffer__Ca_SLHigh_jct;
        NV_Ith_S(ydot, 31) = d_dt_Ca_buffer__Ca_Calsequestrin;
        NV_Ith_S(ydot, 32) = d_dt_cytosolic_Ca_buffer__Ca_TroponinC;
        NV_Ith_S(ydot, 33) = d_dt_cytosolic_Ca_buffer__Ca_TroponinC_Ca_Mg;
        NV_Ith_S(ydot, 34) = d_dt_cytosolic_Ca_buffer__Mg_TroponinC_Ca_Mg;
        NV_Ith_S(ydot, 35) = d_dt_cytosolic_Ca_buffer__Ca_Calmodulin;
        NV_Ith_S(ydot, 36) = d_dt_cytosolic_Ca_buffer__Ca_Myosin;
        NV_Ith_S(ydot, 37) = d_dt_cytosolic_Ca_buffer__Mg_Myosin;
        NV_Ith_S(ydot, 38) = d_dt_cytosolic_Ca_buffer__Ca_SRB;
        NV_Ith_S(ydot, 39) = d_dt_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_Cytosol;
        NV_Ith_S(ydot, 40) = d_dt_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_SL;
        NV_Ith_S(ydot, 41) = d_dt_indo_fluo_Ca_buffer_not_connected__Ca_Indo1_jct;
        NV_Ith_S(ydot, 42) = d_dt_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_Cytosol;
        NV_Ith_S(ydot, 43) = d_dt_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_SL;
        NV_Ith_S(ydot, 44) = d_dt_indo_fluo_Ca_buffer_not_connected__Ca_Fluo3_jct;
    }

};


template<>
void OdeSystemInformation<CvOdeCellShannon2004FromCellML>::Initialise(void)
{
    // Time units: millisecond
    //
    this->mVariableNames.push_back("membrane_voltage");
    this->mVariableUnits.push_back("millivolt");
    this->mInitialConditions.push_back(-85.719687955637);

    this->mVariableNames.push_back("h");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.987140350343);

    this->mVariableNames.push_back("j");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.991822731369);

    this->mVariableNames.push_back("m");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.001370685156);

    this->mVariableNames.push_back("Xr");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.008471550841);

    this->mVariableNames.push_back("Xs");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.00687399199);

    this->mVariableNames.push_back("X_tos");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.004011272375);

    this->mVariableNames.push_back("Y_tos");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.293519921626);

    this->mVariableNames.push_back("R_tos");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.383430556383);

    this->mVariableNames.push_back("X_tof");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.00401120993);

    this->mVariableNames.push_back("Y_tof");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.9946314893);

    this->mVariableNames.push_back("d");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.000006997531);

    this->mVariableNames.push_back("f");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.000675515962);

    this->mVariableNames.push_back("fCaB_SL");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.015352888928);

    this->mVariableNames.push_back("fCaB_jct");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.024609183734);

    this->mVariableNames.push_back("R");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.884673513138);

    this->mVariableNames.push_back("I");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.00000009272);

    this->mVariableNames.push_back("O");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.000000711264);

    this->mVariableNames.push_back("Na_SL");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(8.874077316753);

    this->mVariableNames.push_back("Na_jct");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(8.872823559072);

    this->mVariableNames.push_back("Na_SL_buf");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.776121392467);

    this->mVariableNames.push_back("Na_jct_buf");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(3.557055389701);

    this->mVariableNames.push_back("Nai");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(8.874461106492);

    this->mVariableNames.push_back("Ca_SR");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.545611267699);

    this->mVariableNames.push_back("Ca_SL");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.000106395937);

    this->mVariableNames.push_back("Ca_jct");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.000174843061);

    this->mVariableNames.push_back("Cai");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.000087350002);

    this->mVariableNames.push_back("Ca_SLB_SL");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.009868629147);

    this->mVariableNames.push_back("Ca_SLB_jct");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.007780801995);

    this->mVariableNames.push_back("Ca_SLHigh_SL");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.114438990328);

    this->mVariableNames.push_back("Ca_SLHigh_jct");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.077503874257);

    this->mVariableNames.push_back("Ca_Calsequestrin");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(1.186496899338);

    this->mVariableNames.push_back("Ca_TroponinC");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.008963736337);

    this->mVariableNames.push_back("Ca_TroponinC_Ca_Mg");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.117995194438);

    this->mVariableNames.push_back("Mg_TroponinC_Ca_Mg");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.010337654274);

    this->mVariableNames.push_back("Ca_Calmodulin");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.000295961245);

    this->mVariableNames.push_back("Ca_Myosin");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.001984672275);

    this->mVariableNames.push_back("Mg_Myosin");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.137497736234);

    this->mVariableNames.push_back("Ca_SRB");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.002177112381);

    this->mVariableNames.push_back("Ca_Indo1_Cytosol");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("Ca_Indo1_SL");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("Ca_Indo1_jct");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("Ca_Fluo3_Cytosol");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("Ca_Fluo3_SL");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("Ca_Fluo3_jct");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0);

    this->mInitialised = true;
}

#endif //_CVODECELLSHANNON2004FROMCELLML_

#endif // CHASTE_CVODE

