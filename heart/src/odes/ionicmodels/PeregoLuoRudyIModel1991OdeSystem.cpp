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
#include "PeregoLuoRudyIModel1991OdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include <cmath>
#include "Exception.hpp"

//
// Model-scope constant parameters
//
const double PeregoLuoRudyIModel1991OdeSystem::membrane_C = 1.0;
const double PeregoLuoRudyIModel1991OdeSystem::membrane_F = 96484.6;
const double PeregoLuoRudyIModel1991OdeSystem::membrane_R = 8314;
const double PeregoLuoRudyIModel1991OdeSystem::membrane_T = 310.0;
const double PeregoLuoRudyIModel1991OdeSystem::background_current_E_b = -59.87;
const double PeregoLuoRudyIModel1991OdeSystem::background_current_g_b = 0.03921;
const double PeregoLuoRudyIModel1991OdeSystem::fast_sodium_current_g_Na = 23.0;
const double PeregoLuoRudyIModel1991OdeSystem::ionic_concentrations_Ki = 145.0;
const double PeregoLuoRudyIModel1991OdeSystem::ionic_concentrations_Ko = 5.4;
const double PeregoLuoRudyIModel1991OdeSystem::ionic_concentrations_Nai = 18.0;
const double PeregoLuoRudyIModel1991OdeSystem::ionic_concentrations_Nao = 140.0;
const double PeregoLuoRudyIModel1991OdeSystem::plateau_potassium_current_g_Kp = 0.0183;
const double PeregoLuoRudyIModel1991OdeSystem::time_dependent_potassium_current_PR_NaK = 0.01833;



/**
 * Constructor
 */
PeregoLuoRudyIModel1991OdeSystem::PeregoLuoRudyIModel1991OdeSystem(
    boost::shared_ptr<AbstractIvpOdeSolver> pSolver, // unused.
    boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus,
    bool useAdaptTimestep)
        : AbstractPeregoCardiacCell(8, 4, pIntracellularStimulus, useAdaptTimestep)
{
    
    assert(mGatingVariableIndices.size() == 0);
    mGatingVariableIndices.push_back(0);
    mGatingVariableIndices.push_back(1);
    mGatingVariableIndices.push_back(2);
    mGatingVariableIndices.push_back(5);
    mGatingVariableIndices.push_back(6);
    mGatingVariableIndices.push_back(7);
    assert(mGatingVariableIndices.size() == 6);
    
    mpSystemInfo = OdeSystemInformation<PeregoLuoRudyIModel1991OdeSystem>::Instance();
    // set the final paramter
    fast_sodium_current_E_Na = ((membrane_R * membrane_T) / membrane_F) *
                               log(ionic_concentrations_Nao / ionic_concentrations_Nai);
                       
    
    ma_current.resize(8);
    mb_current.resize(8);
    ma_predicted.resize(8);
    mb_predicted.resize(8);
    ma_error.resize(8);
    mb_error.resize(8);
    mWeightedErrorTolerances.resize(8);
    
    //set the tolerance weight to the default value. Can be modified by set method.
    double tolerance_weight = 1e-2;
    
    mWeightedErrorTolerances[0] = 1 * tolerance_weight; // Gate h error
    mWeightedErrorTolerances[1] = 1 * tolerance_weight; // Gate j error
    mWeightedErrorTolerances[2] = 1 * tolerance_weight; // Gate m error
    mWeightedErrorTolerances[3] = 7e-3 * tolerance_weight; // Calcium error
    mWeightedErrorTolerances[4] = 84 * tolerance_weight; // Voltage error
    mWeightedErrorTolerances[5] = 1 * tolerance_weight; // Gate d error
    mWeightedErrorTolerances[6] = 1 * tolerance_weight; // Gate f error
    mWeightedErrorTolerances[7] = 1 * tolerance_weight; // Gate X error

    Init();
}

/**
 * Destructor
 */
PeregoLuoRudyIModel1991OdeSystem::~PeregoLuoRudyIModel1991OdeSystem(void)
{
}

void PeregoLuoRudyIModel1991OdeSystem::ComputeSystemParameters(const std::vector<double>& rY, double currentTime)
{
    double fast_sodium_current_h_gate_h = rY[0];
    double fast_sodium_current_j_gate_j = rY[1];
    double fast_sodium_current_m_gate_m = rY[2];
    double intracellular_calcium_concentration_Cai = rY[3];
    double membrane_V = rY[4];
    double slow_inward_current_d_gate_d = rY[5];
    double slow_inward_current_f_gate_f = rY[6];
    double time_dependent_potassium_current_X_gate_X = rY[7];

    //PRINT_VECTOR(rY);
    for(unsigned i=0;i<8;i++)
    {
        assert(!std::isnan(rY[i]));
    }
    
    if (mIsTheCorrectorStep == false && mIsTheErrorEvaluationStep == false && mIsThereTooMuchError == false)
    {
        ma_previous = ma_current;
        mb_previous = mb_current;
    }
    
    
    VerifyStateVariables();

    double background_current_i_b = background_current_g_b*(membrane_V-background_current_E_b);

    double fast_sodium_current_h_gate_alpha_h;
    
    if (membrane_V < -40.0)
    {
        fast_sodium_current_h_gate_alpha_h = 0.135*exp((80.0+membrane_V)/-6.8);
    }
    else
    {
        fast_sodium_current_h_gate_alpha_h = 0.0;
    }

    double fast_sodium_current_h_gate_beta_h;

    if (membrane_V < -40.0)
    {
        fast_sodium_current_h_gate_beta_h = 3.56*exp(0.079*membrane_V)+3.1e5*exp(0.35*membrane_V);
    }
    else
    {
        fast_sodium_current_h_gate_beta_h = 1.0/(0.13*(1.0+exp((membrane_V+10.66)/-11.1)));
    }

    double fast_sodium_current_j_gate_alpha_j;

    if (membrane_V < -40.0)
    {
        fast_sodium_current_j_gate_alpha_j = (-1.2714e5*exp(0.2444*membrane_V)-3.474e-5*exp(-0.04391*membrane_V))*(membrane_V+37.78)/(1.0+exp(0.311*(membrane_V+79.23)));
    }
    else
    {
        fast_sodium_current_j_gate_alpha_j = 0.0;
    }

    double fast_sodium_current_j_gate_beta_j;

    if (membrane_V < -40.0)
    {
        fast_sodium_current_j_gate_beta_j = 0.1212*exp(-0.01052*membrane_V)/(1.0+exp(-0.1378*(membrane_V+40.14)));
    }
    else
    {
        fast_sodium_current_j_gate_beta_j = 0.3*exp(-2.535e-7*membrane_V)/(1.0+exp(-0.1*(membrane_V+32.0)));
    }

    double fast_sodium_current_m_gate_alpha_m = 0.32*(membrane_V+47.13)/(1.0-exp(-0.1*(membrane_V+47.13)));
    double fast_sodium_current_m_gate_beta_m = 0.08*exp(-membrane_V/11.0);

    double fast_sodium_current_i_Na = fast_sodium_current_g_Na*pow(fast_sodium_current_m_gate_m, 3.0)*fast_sodium_current_h_gate_h*fast_sodium_current_j_gate_j*(membrane_V-fast_sodium_current_E_Na);

    double slow_inward_current_d_gate_alpha_d = 0.095*exp(-0.01*(membrane_V-5.0))/(1.0+exp(-0.072*(membrane_V-5.0)));
    double slow_inward_current_d_gate_beta_d = 0.07*exp(-0.017*(membrane_V+44.0))/(1.0+exp(0.05*(membrane_V+44.0)));


    double slow_inward_current_f_gate_alpha_f = 0.012*exp(-0.008*(membrane_V+28.0))/(1.0+exp(0.15*(membrane_V+28.0)));
    double slow_inward_current_f_gate_beta_f = 0.0065*exp(-0.02*(membrane_V+30.0))/(1.0+exp(-0.2*(membrane_V+30.0)));


    double slow_inward_current_E_si = 7.7-13.0287*log(intracellular_calcium_concentration_Cai);
    double slow_inward_current_i_si = 0.09*slow_inward_current_d_gate_d*slow_inward_current_f_gate_f*(membrane_V-slow_inward_current_E_si);
    double intracellular_calcium_concentration_Cai_prime = -1e-4*slow_inward_current_i_si+0.07*(1e-4-intracellular_calcium_concentration_Cai);
    double time_dependent_potassium_current_g_K = 0.282*sqrt(ionic_concentrations_Ko/5.4);

    double time_dependent_potassium_current_Xi_gate_Xi;

    if (membrane_V > -100.0)
    {
        time_dependent_potassium_current_Xi_gate_Xi = 2.837*(exp(0.04*(membrane_V+77.0))-1.0)/((membrane_V+77.0)*exp(0.04*(membrane_V+35.0)));
    }
    else
    {
        #define COVERAGE_IGNORE
        time_dependent_potassium_current_Xi_gate_Xi = 1.0;
        #undef COVERAGE_IGNORE
    }

    double time_dependent_potassium_current_X_gate_alpha_X = 0.0005*exp(0.083*(membrane_V+50.0))/(1.0+exp(0.057*(membrane_V+50.0)));
    double time_dependent_potassium_current_X_gate_beta_X = 0.0013*exp(-0.06*(membrane_V+20.0))/(1.0+exp(-0.04*(membrane_V+20.0)));


    double time_dependent_potassium_current_E_K = ((membrane_R*membrane_T)/membrane_F)*log((ionic_concentrations_Ko+time_dependent_potassium_current_PR_NaK*ionic_concentrations_Nao)/(ionic_concentrations_Ki+time_dependent_potassium_current_PR_NaK*ionic_concentrations_Nai));
    double time_dependent_potassium_current_i_K = time_dependent_potassium_current_g_K*time_dependent_potassium_current_X_gate_X*time_dependent_potassium_current_Xi_gate_Xi*(membrane_V-time_dependent_potassium_current_E_K);

    double time_independent_potassium_current_g_K1 = 0.6047*sqrt(ionic_concentrations_Ko/5.4);
    double time_independent_potassium_current_E_K1 = ((membrane_R*membrane_T)/membrane_F)*log(ionic_concentrations_Ko/ionic_concentrations_Ki);
    double time_independent_potassium_current_K1_gate_alpha_K1 = 1.02/(1.0+exp(0.2385*(membrane_V-time_independent_potassium_current_E_K1-59.215)));
    double time_independent_potassium_current_K1_gate_beta_K1 = (0.49124*exp(0.08032*(membrane_V+5.476-time_independent_potassium_current_E_K1))+exp(0.06175*(membrane_V-(time_independent_potassium_current_E_K1+594.31))))/(1.0+exp(-0.5143*(membrane_V-time_independent_potassium_current_E_K1+4.753)));
    double time_independent_potassium_current_K1_gate_K1_infinity = time_independent_potassium_current_K1_gate_alpha_K1/(time_independent_potassium_current_K1_gate_alpha_K1+time_independent_potassium_current_K1_gate_beta_K1);
    double time_independent_potassium_current_i_K1 = time_independent_potassium_current_g_K1*time_independent_potassium_current_K1_gate_K1_infinity*(membrane_V-time_independent_potassium_current_E_K1);
    double plateau_potassium_current_Kp = 1.0/(1.0+exp((7.488-membrane_V)/5.98));
    double plateau_potassium_current_E_Kp = time_independent_potassium_current_E_K1;
    double plateau_potassium_current_i_Kp = plateau_potassium_current_g_Kp*plateau_potassium_current_Kp*(membrane_V-plateau_potassium_current_E_Kp);
    

    double i_stim;
    if (mIsTheCorrectorStep == false)
    {
        i_stim = GetStimulus(currentTime);
    }
    else
    {
        i_stim = GetStimulus(currentTime+mLocalTimeStep);
    };
    

    //calculate dV
    assert(!std::isnan(i_stim));
    assert(!std::isnan(fast_sodium_current_i_Na));
    assert(!std::isnan(slow_inward_current_i_si));
    assert(!std::isnan(time_dependent_potassium_current_i_K));
    assert(!std::isnan(time_independent_potassium_current_i_K1));
    assert(!std::isnan(plateau_potassium_current_i_Kp));
    assert(!std::isnan(background_current_i_b));
    double membrane_V_prime = (-1.0/membrane_C)*(fast_sodium_current_i_Na+slow_inward_current_i_si+time_dependent_potassium_current_i_K+time_independent_potassium_current_i_K1+plateau_potassium_current_i_Kp+background_current_i_b + i_stim);
    //PRINT_2_VARIABLES(membrane_V_prime,membrane_V);
    assert(!std::isnan(membrane_V_prime));
    
    // do not update voltage if the mSetVoltageDerivativeToZero flag has been set
//    if (mSetVoltageDerivativeToZero)
//    {
//        //NEVER_REACHED;
//        // this hasn't been tested in tissue yet
//
//        ///\ TODO Remove coverage ignore when ComputeExceptVoltage method will be implemented in the abstract class  
//        //#define COVERAGE_IGNORE
//        membrane_V_prime = 0;
//        //#undef COVERAGE_IGNORE
//    }
    
    if (mIsTheCorrectorStep == false && mIsTheErrorEvaluationStep == false)
    {
        // Compute the parameters for the gating variable updates...
        ma_current[0]= - fast_sodium_current_h_gate_alpha_h - fast_sodium_current_h_gate_beta_h;
        ma_current[1]= - fast_sodium_current_j_gate_alpha_j - fast_sodium_current_j_gate_beta_j;
        ma_current[2]= - fast_sodium_current_m_gate_alpha_m - fast_sodium_current_m_gate_beta_m;
        ma_current[5]= - slow_inward_current_d_gate_alpha_d - slow_inward_current_d_gate_beta_d;
        ma_current[6]= - slow_inward_current_f_gate_alpha_f - slow_inward_current_f_gate_beta_f;
        ma_current[7]= - time_dependent_potassium_current_X_gate_alpha_X - time_dependent_potassium_current_X_gate_beta_X;
    
        mb_current[0] = fast_sodium_current_h_gate_alpha_h;
        mb_current[1] = fast_sodium_current_j_gate_alpha_j;
        mb_current[2] = fast_sodium_current_m_gate_alpha_m;
        mb_current[5] = slow_inward_current_d_gate_alpha_d;
        mb_current[6] = slow_inward_current_f_gate_alpha_f;
        mb_current[7] = time_dependent_potassium_current_X_gate_alpha_X;
        
        // ...and add to ma_current the derivatives of the voltage and the calcium concentration    
        ma_current[4] = membrane_V_prime;
        ma_current[3] = intracellular_calcium_concentration_Cai_prime;
        //PRINT_VECTOR(ma_current);
        //PRINT_VECTOR(mb_current);  
    }
    if (mIsTheCorrectorStep == true && mIsTheErrorEvaluationStep == false)
    {

        // Compute the parameters for the gating variable updates...
        ma_predicted[0]= - fast_sodium_current_h_gate_alpha_h - fast_sodium_current_h_gate_beta_h;
        ma_predicted[1]= - fast_sodium_current_j_gate_alpha_j - fast_sodium_current_j_gate_beta_j;
        ma_predicted[2]= - fast_sodium_current_m_gate_alpha_m - fast_sodium_current_m_gate_beta_m;
        ma_predicted[5]= - slow_inward_current_d_gate_alpha_d - slow_inward_current_d_gate_beta_d;
        ma_predicted[6]= - slow_inward_current_f_gate_alpha_f - slow_inward_current_f_gate_beta_f;
        ma_predicted[7]= - time_dependent_potassium_current_X_gate_alpha_X - time_dependent_potassium_current_X_gate_beta_X;
    
        mb_predicted[0] = fast_sodium_current_h_gate_alpha_h;
        mb_predicted[1] = fast_sodium_current_j_gate_alpha_j;
        mb_predicted[2] = fast_sodium_current_m_gate_alpha_m;
        mb_predicted[5] = slow_inward_current_d_gate_alpha_d;
        mb_predicted[6] = slow_inward_current_f_gate_alpha_f;
        mb_predicted[7] = time_dependent_potassium_current_X_gate_alpha_X;
        
        // ...and add to ma_predicted the derivatives of the voltage and the calcium concentration    
        ma_predicted[4] = membrane_V_prime;
        ma_predicted[3] = intracellular_calcium_concentration_Cai_prime;
        //PRINT_VECTOR(ma_predicted);
        //PRINT_VECTOR(mb_predicted);        
    }
    if (mIsTheErrorEvaluationStep==true)
    {
        
        // Compute the parameters for the gating variable updates...
        ma_error[0]= - fast_sodium_current_h_gate_alpha_h - fast_sodium_current_h_gate_beta_h;
        ma_error[1]= - fast_sodium_current_j_gate_alpha_j - fast_sodium_current_j_gate_beta_j;
        ma_error[2]= - fast_sodium_current_m_gate_alpha_m - fast_sodium_current_m_gate_beta_m;
        ma_error[5]= - slow_inward_current_d_gate_alpha_d - slow_inward_current_d_gate_beta_d;
        ma_error[6]= - slow_inward_current_f_gate_alpha_f - slow_inward_current_f_gate_beta_f;
        ma_error[7]= - time_dependent_potassium_current_X_gate_alpha_X - time_dependent_potassium_current_X_gate_beta_X;
    
        mb_error[0] = fast_sodium_current_h_gate_alpha_h;
        mb_error[1] = fast_sodium_current_j_gate_alpha_j;
        mb_error[2] = fast_sodium_current_m_gate_alpha_m;
        mb_error[5] = slow_inward_current_d_gate_alpha_d;
        mb_error[6] = slow_inward_current_f_gate_alpha_f;
        mb_error[7] = time_dependent_potassium_current_X_gate_alpha_X;
        
        // ...and add to ma_error the derivatives of the voltage and the calcium concentration    
        ma_error[4] = membrane_V_prime;
        ma_error[3] = intracellular_calcium_concentration_Cai_prime;
        
        //PRINT_VECTOR(ma_error);
        //PRINT_VECTOR(mb_error);
    }
    
    if (mIsTheFirstStep == true)
    {
        ma_previous = ma_current;
        mb_previous = mb_current;
        
        mIsTheFirstStep = false;
    }
    
}

void PeregoLuoRudyIModel1991OdeSystem::SetToleranceWeight (double tolerance_weight) 
{   
    mWeightedErrorTolerances[0] = 1 * tolerance_weight; // Gate h error
    mWeightedErrorTolerances[1] = 1 * tolerance_weight; // Gate j error
    mWeightedErrorTolerances[2] = 1 * tolerance_weight; // Gate m error
    mWeightedErrorTolerances[3] = 7e-3 * tolerance_weight; // Calcium error
    mWeightedErrorTolerances[4] = 84 * tolerance_weight; // Voltage error
    mWeightedErrorTolerances[5] = 1 * tolerance_weight; // Gate d error
    mWeightedErrorTolerances[6] = 1 * tolerance_weight; // Gate f error
    mWeightedErrorTolerances[7] = 1 * tolerance_weight; // Gate X error
}

double PeregoLuoRudyIModel1991OdeSystem::GetIIonic()
{
    double fast_sodium_current_h_gate_h = mStateVariables[0];
    double fast_sodium_current_j_gate_j = mStateVariables[1];
    double fast_sodium_current_m_gate_m = mStateVariables[2];
    double intracellular_calcium_concentration_Cai = mStateVariables[3];
    double membrane_V = mStateVariables[4];
    double slow_inward_current_d_gate_d = mStateVariables[5];
    double slow_inward_current_f_gate_f = mStateVariables[6];
    double time_dependent_potassium_current_X_gate_X = mStateVariables[7];

    /*
     * Compute the PeregoLuoRudyIModel1991OdeSystem model
     */
    double background_current_i_b = background_current_g_b*(membrane_V-background_current_E_b);

    double fast_sodium_current_i_Na = fast_sodium_current_g_Na*pow(fast_sodium_current_m_gate_m, 3.0)*fast_sodium_current_h_gate_h*fast_sodium_current_j_gate_j*(membrane_V-fast_sodium_current_E_Na);

    double slow_inward_current_E_si = 7.7-13.0287*log(intracellular_calcium_concentration_Cai);
    double slow_inward_current_i_si = 0.09*slow_inward_current_d_gate_d*slow_inward_current_f_gate_f*(membrane_V-slow_inward_current_E_si);

    double time_dependent_potassium_current_g_K = 0.282*sqrt(ionic_concentrations_Ko/5.4);
    double time_dependent_potassium_current_Xi_gate_Xi;

    // Although the equation below looks strange (particularly the arguments of the
    // exponentials, it is in fact correct.
    if (membrane_V > -100.0)
    {
        time_dependent_potassium_current_Xi_gate_Xi = 2.837*(exp(0.04*(membrane_V+77.0))-1.0)/((membrane_V+77.0)*exp(0.04*(membrane_V+35.0)));
    }
    else
    {
        #define COVERAGE_IGNORE
        time_dependent_potassium_current_Xi_gate_Xi = 1.0;
        #undef COVERAGE_IGNORE
    }

    double time_dependent_potassium_current_E_K = ((membrane_R*membrane_T)/membrane_F)*log((ionic_concentrations_Ko+time_dependent_potassium_current_PR_NaK*ionic_concentrations_Nao)/(ionic_concentrations_Ki+time_dependent_potassium_current_PR_NaK*ionic_concentrations_Nai));
    double time_dependent_potassium_current_i_K = time_dependent_potassium_current_g_K*time_dependent_potassium_current_X_gate_X*time_dependent_potassium_current_Xi_gate_Xi*(membrane_V-time_dependent_potassium_current_E_K);

    double time_independent_potassium_current_g_K1 = 0.6047*sqrt(ionic_concentrations_Ko/5.4);
    double time_independent_potassium_current_E_K1 =((membrane_R*membrane_T)/membrane_F)*log(ionic_concentrations_Ko/ionic_concentrations_Ki);
    double time_independent_potassium_current_K1_gate_alpha_K1 = 1.02/(1.0+exp(0.2385*(membrane_V-time_independent_potassium_current_E_K1-59.215)));
    double time_independent_potassium_current_K1_gate_beta_K1 = (0.49124*exp(0.08032*(membrane_V+5.476-time_independent_potassium_current_E_K1))+exp(0.06175*(membrane_V-(time_independent_potassium_current_E_K1+594.31))))/(1.0+exp(-0.5143*(membrane_V-time_independent_potassium_current_E_K1+4.753)));
    double time_independent_potassium_current_K1_gate_K1_infinity = time_independent_potassium_current_K1_gate_alpha_K1/(time_independent_potassium_current_K1_gate_alpha_K1+time_independent_potassium_current_K1_gate_beta_K1);
    double time_independent_potassium_current_i_K1 = time_independent_potassium_current_g_K1*time_independent_potassium_current_K1_gate_K1_infinity*(membrane_V-time_independent_potassium_current_E_K1);

    double plateau_potassium_current_Kp = 1.0/(1.0+exp((7.488-membrane_V)/5.98));
    double plateau_potassium_current_E_Kp = time_independent_potassium_current_E_K1;
    double plateau_potassium_current_i_Kp = plateau_potassium_current_g_Kp*plateau_potassium_current_Kp*(membrane_V-plateau_potassium_current_E_Kp);

    double i_ionic = fast_sodium_current_i_Na+slow_inward_current_i_si+time_dependent_potassium_current_i_K+time_independent_potassium_current_i_K1+plateau_potassium_current_i_Kp+background_current_i_b;

    assert(!std::isnan(i_ionic));
    return i_ionic;
}

void PeregoLuoRudyIModel1991OdeSystem::VerifyStateVariables()
{
    const std::vector<double>& rY = rGetStateVariables();

    const double fast_sodium_current_h_gate_h = rY[0];            // gating
    const double fast_sodium_current_j_gate_j = rY[1];            // gating
    const double fast_sodium_current_m_gate_m = rY[2];            // gating
    const double intracellular_calcium_concentration_Cai = rY[3]; // concentration
    const double slow_inward_current_d_gate_d = rY[5];            // gating
    const double slow_inward_current_f_gate_f = rY[6];            // gating
    const double time_dependent_potassium_current_X_gate_X = rY[7]; // gating

    #define COVERAGE_IGNORE
    if (!(0.0<=fast_sodium_current_h_gate_h && fast_sodium_current_h_gate_h<=1.0))
    {
        EXCEPTION(DumpState("h gate for fast sodium current has gone out of range. Check model parameters, for example spatial stepsize"));
    }

    if (!(0.0<=fast_sodium_current_j_gate_j && fast_sodium_current_j_gate_j<=1.0))
    {
        EXCEPTION(DumpState("j gate for fast sodium current has gone out of range. Check model parameters, for example spatial stepsize"));
    }

    if (!(0.0<=fast_sodium_current_m_gate_m && fast_sodium_current_m_gate_m<=1.0))
    {
        EXCEPTION(DumpState("m gate for fast sodium current has gone out of range. Check model parameters, for example spatial stepsize"));
    }

    if (!(0.0<intracellular_calcium_concentration_Cai))
    {
        EXCEPTION(DumpState("intracellular_calcium_concentration_Cai has become non-positive, ie gone out of range. Check model parameters, for example spatial stepsize"));
    }

    if (!(0.0<=slow_inward_current_d_gate_d && slow_inward_current_d_gate_d<=1.0))
    {
        EXCEPTION(DumpState("d gate for slow inward current has gone out of range. Check model parameters, for example spatial stepsize"));
    }

    if (!(0.0<=slow_inward_current_f_gate_f && slow_inward_current_f_gate_f<=1.0))
    {
        EXCEPTION(DumpState("f gate for slow inward current has gone out of range. Check model parameters, for example spatial stepsize"));
    }

    if (!(0.0<=time_dependent_potassium_current_X_gate_X && time_dependent_potassium_current_X_gate_X<=1.0))
    {
        EXCEPTION(DumpState("X gate for time dependent potassium current has gone out of range. Check model parameters, for example spatial stepsize"));
    }
    #undef COVERAGE_IGNORE
}




template<>
void OdeSystemInformation<PeregoLuoRudyIModel1991OdeSystem>::Initialise(void)
{
    // State variables
    this->mVariableNames.push_back("h");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.9804713);

    this->mVariableNames.push_back("j");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.98767124);

    this->mVariableNames.push_back("m");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.00187018);

    this->mVariableNames.push_back("CaI");
    this->mVariableUnits.push_back("mMol");
    this->mInitialConditions.push_back(0.0002);

    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("mV");
    this->mInitialConditions.push_back(-83.853);

    this->mVariableNames.push_back("d");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.00316354);

    this->mVariableNames.push_back("f");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.99427859);

    this->mVariableNames.push_back("x");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.16647703);

    this->mInitialised = true;
}
