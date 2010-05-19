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
#include "DiFrancescoNoble1985OdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include <cmath>
#include <cstdio>
//#include <iostream>
#include "Exception.hpp"

//
// Model-scope constant parameters
//
const double DiFrancescoNoble1985OdeSystem::membrane_C = 75.0;
const double DiFrancescoNoble1985OdeSystem::radius = 0.05;
const double DiFrancescoNoble1985OdeSystem::length = 2;
const double DiFrancescoNoble1985OdeSystem::V_e_ratio = 0.1;
const double DiFrancescoNoble1985OdeSystem::R = 8314.472;
const double DiFrancescoNoble1985OdeSystem::T = 310;
const double DiFrancescoNoble1985OdeSystem::F = 96485.3415;
// conductances and the likes
const double DiFrancescoNoble1985OdeSystem::g_fna = 3;
const double DiFrancescoNoble1985OdeSystem::g_fk = 3;
const double DiFrancescoNoble1985OdeSystem::g_na_b = 0.18;
const double DiFrancescoNoble1985OdeSystem::g_ca_b = 0.02;
const double DiFrancescoNoble1985OdeSystem::g_na = 750.0;
const double DiFrancescoNoble1985OdeSystem::g_k1 = 920;
const double DiFrancescoNoble1985OdeSystem::g_to = 0.28;
const double DiFrancescoNoble1985OdeSystem::I_P = 125;
const double DiFrancescoNoble1985OdeSystem::i_kmax = 180;
const double DiFrancescoNoble1985OdeSystem::k_naca = 0.02;
const double DiFrancescoNoble1985OdeSystem::P_si = 15.0;
// concentrations
const double DiFrancescoNoble1985OdeSystem::Nao = 140;
const double DiFrancescoNoble1985OdeSystem::Cao = 2;
const double DiFrancescoNoble1985OdeSystem::Kb = 4;
const double DiFrancescoNoble1985OdeSystem::Kmf = 45;
const double DiFrancescoNoble1985OdeSystem::Km1 = 210;
const double DiFrancescoNoble1985OdeSystem::Kmto = 10;
const double DiFrancescoNoble1985OdeSystem::KmK = 1;
const double DiFrancescoNoble1985OdeSystem::KmNa = 40;
const double DiFrancescoNoble1985OdeSystem::Kmf2 = 0.001;
const double DiFrancescoNoble1985OdeSystem::KmCa = 0.001;
const double DiFrancescoNoble1985OdeSystem::Km_Ca = 0.0005;
// others
const double DiFrancescoNoble1985OdeSystem::n_naca = 3;
const double DiFrancescoNoble1985OdeSystem::gamma = 0.5;
const double DiFrancescoNoble1985OdeSystem::d_naca = 0.001;
const double DiFrancescoNoble1985OdeSystem::rCa = 2.0;
const double DiFrancescoNoble1985OdeSystem::tau_up = 0.025;
const double DiFrancescoNoble1985OdeSystem::tau_rep = 2;
const double DiFrancescoNoble1985OdeSystem::tau_rel = 0.05;
const double DiFrancescoNoble1985OdeSystem::Ca_up_max = 5;
const double DiFrancescoNoble1985OdeSystem::pf = 0.0007;
const double DiFrancescoNoble1985OdeSystem::Vecs = 0.05;


/*Constructor*/
DiFrancescoNoble1985OdeSystem::DiFrancescoNoble1985OdeSystem(
        boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
        boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : AbstractCardiacCell(pSolver, 16, 0, pIntracellularStimulus)
{
    mpSystemInfo = OdeSystemInformation<DiFrancescoNoble1985OdeSystem>::Instance();

    //compute the fixed parameters
    Vcell= 3.141592654*pow(radius,2)*length;
    Vi=Vcell*(1.0-V_e_ratio);
    Vup=0.05*Vi;
    Vrel=0.02*Vi;
    Ve=Vi/(1.0-Vecs);
    RToNF=R*T/F;

    //initialise
    Init();
}

/**
 * Destructor
 */
DiFrancescoNoble1985OdeSystem::~DiFrancescoNoble1985OdeSystem(void)
{
}

/**
 * Fill in a vector representing the RHS of the LuoRudyIModel1991OdeSystem system
 * of Odes at each time step, y' = [y1' ... yn'].
 * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
 *
 * @param time  the current time, in milliseconds
 * @param rY  current values of the state variables
 * @param rDY  to be filled in with derivatives
 */
void DiFrancescoNoble1985OdeSystem::EvaluateYDerivatives(double time,
                                                      const std::vector<double> &rY,
                                                      std::vector<double> &rDY)
{
    double Vm = rY[0];
    double y_gate = rY[1];
    double s_gate = rY[2];
    double m_gate = rY[3];
    double h_gate = rY[4];
    double d_gate= rY[5];
    double f_gate = rY[6];
    double f2_gate = rY[7];
    double p = rY[8];
    double Nai = rY[9];
    double Cai = rY[10];
    double Kc = rY[11];
    double Ki = rY[12];
    double x_gate = rY[13];
    double Ca_up = rY[14];
    double Ca_rel = rY[15];
    VerifyStateVariables();

   ////////////////////////////////////////////////////////
   // gating variables and ionic fluxes have been scaled by a factor of 1,000 to convert them into milliseconds
   ////////////////////////////////////////////////////////
    double Emh=RToNF*log((Nao+0.12*Kc)/(Nai+0.12*Ki));
    double E_na=RToNF*log(Nao/Nai);
    double E_k=RToNF*log(Kc/Ki);
    double E_ca=0.5*RToNF*log(Cao/Cai);
    //i_f parameters
    double i_fk=(Kc/(Kc+Kmf))*(g_fk*(Vm-E_k));
    double i_fna=(Kc/(Kc+Kmf))*g_fna*(Vm-E_na);
    double alpha_y=0.001*0.05*exp(-0.067*(Vm+52.0-10.0));
    double beta_y;

   if(fabs(Vm+52.0-10.0)<=0.0001) //i.e. if Vm = (-52)
    {
        #define COVERAGE_IGNORE
        beta_y=2.5*0.001;
        #undef COVERAGE_IGNORE
    }
    else
    {
        beta_y=0.001*(Vm+52.0-10.0)/(1.0-exp(-0.2*(Vm+52.0-10.0)));
    }
    double y_gate_prime=alpha_y*(1.0-y_gate)-beta_y*y_gate;

    //i_k parameters
    double I_K=i_kmax*(Ki-Kc*exp(-Vm/RToNF))/140.0;
    double alpha_x=0.001*0.5*exp(0.0826*(Vm+50.0))/(1.0+exp(0.057*(Vm+50.0)));
    double beta_x=0.001*1.3*exp(-0.06*(Vm+20.0))/(1.0+exp(-0.04*(Vm+20.0)));
    double x_gate_prime=alpha_x*(1.0-x_gate)-beta_x*x_gate;

    //i_to parameters
    double alpha_s=0.001*0.033*exp(-Vm/17.0);
    double beta_s=0.001*33.0/(1.0+exp(-(Vm+10.0)/8.0));
    double s_gate_prime=alpha_s*(1.0-s_gate)-beta_s*s_gate;

    // i_na parameters
    double alpha_m;
    if(fabs(Vm+41.0)<=0.00001) //i.e. if Vm = -41
    {
      #define COVERAGE_IGNORE
      alpha_m=2000*0.001;
      #undef COVERAGE_IGNORE
    }
    else
    {
        alpha_m=0.001*200.0*(Vm+41.0)/(1.0-exp(-0.1*(Vm+41.0)));
    }

    double beta_m=0.001*8000*exp(-0.056*(Vm+66.0));

    double m_gate_prime=alpha_m*(1.0-m_gate)-beta_m*m_gate;

    double alpha_h=0.001*20*exp(-0.125*(Vm+75.0));
    double beta_h=0.001*2000.0/(320.0*exp(-0.1*(Vm+75.0))+1.0);
    double h_gate_prime=alpha_h*(1.0-h_gate)-beta_h*h_gate;

    //d gate
    double alpha_d;
    if(fabs(Vm+24.0-5.0)<=0.0001) //i.e. if Vm = -24
    {
        #define COVERAGE_IGNORE
        alpha_d=120*0.001;
        #undef COVERAGE_IGNORE
    }
    else
    {
        alpha_d=0.001*30*(Vm+24.0-5.0)/(1.0-exp(-(Vm+24.0-5.0)/4.0));
    }

    double beta_d;
    if(fabs(Vm+24.0-5.0)<=0.0001) //i.e. if Vm = -24
    {
        #define COVERAGE_IGNORE
        beta_d=120*0.001;
        #undef COVERAGE_IGNORE
    }
    else
    {
        beta_d=0.001*12.0*(Vm+24.0-5.0)/(exp((Vm+24.0-5.0)/10.0)-1.0);
    }
    double d_gate_prime=alpha_d*(1.0-d_gate)-beta_d*d_gate;

    //f gate.
    double alpha_f;
    if(fabs(Vm+34)<=0.0001)
    {
        #define COVERAGE_IGNORE
        alpha_f=25*0.001;
        #undef COVERAGE_IGNORE
    }
    else
    {
        alpha_f=0.001*6.25*(Vm+34.0)/(exp((Vm+34.0)/4.0)-1.0);
    }
    double beta_f=0.001*50.0/(1.0+exp(-(Vm+34.0)/4.0));

    double f_gate_prime=alpha_f*(1.0-f_gate)-beta_f*f_gate;

    //f2_gate
    double alpha_f2=0.001*5.0;
    double beta_f2=0.001*Cai*alpha_f2/Kmf2;
    double f2_gate_prime=alpha_f2-(f2_gate*(alpha_f2+beta_f2));

    //p_gate
    double alpha_p=0.001*(0.625*(Vm+34.0))/(exp((Vm+34.0)/4.0)-1.0);
    double beta_p=0.001*5.0/(1.0+exp(-(Vm+34.0)/4.0));
    double p_prime=alpha_p*(1.0-p)-beta_p*p;


      //i_na
    double i_na=g_na*pow(m_gate,3)*h_gate*(Vm-Emh);

    //i_f
    double i_f=y_gate*(i_fk+i_fna);
    //i_k
    double i_k=x_gate*I_K;
    //i_k1
    double i_k1=(g_k1*Kc/(Kc+Km1))*(Vm-E_k)/(1.0+exp((Vm-E_k+10)*2.0/RToNF));
    //i_to
    double i_to=((((s_gate*g_to*Cai*(0.2+Kc/(Kmto+Kc)))/(Km_Ca+Cai))*(Vm+10.0))/(1.0-exp(-0.2*(Vm+10.0))))*((Ki*exp((Vm*0.5)/RToNF))-(Kc*exp(-0.5*Vm/RToNF)));
    //i_na_b
    double i_na_b=g_na_b*(Vm-E_na);
    //i_p
    double i_p=((I_P*Kc)/(KmK+Kc))*Nai/(Nai+KmNa);
    //i_naca
    double i_naca=(k_naca*((Cao*(exp((gamma*Vm*(n_naca-2))/(RToNF)))*(pow(Nai,n_naca)))-(Cai*(exp(((gamma-1.0)*(n_naca-2.0)*Vm)/(RToNF)))*(pow(Nao,n_naca)))))/((1.0+(d_naca*((Cai*pow(Nao,n_naca))+(Cao*pow(Nai,n_naca)))))*(1.0+(Cai/0.0069)));
    //i_b_ca
    double i_ca_b=g_ca_b*(Vm-E_ca);

    //i_si
    double i_sica=((4.0*P_si*(Vm-50.0))/(RToNF*(1.0-(exp((-2.0*(Vm-50.0))/RToNF)))))*((Cai*exp(100.0/RToNF))-(Cao*exp(-2.0*(Vm-50.0)/RToNF)))*d_gate*f_gate*f2_gate;
    double i_sik=((0.01*P_si*(Vm-50.0))/(RToNF*(1.0-(exp((-(Vm-50.0))/RToNF)))))*((Ki*exp(50.0/RToNF))-(Kc*exp(-(Vm-50.0)/RToNF)))*d_gate*f_gate*f2_gate;
    double i_sina=((0.01*P_si*(Vm-50.0))/(RToNF*(1.0-(exp((-(Vm-50.0))/RToNF)))))*((Nai*exp(50.0/RToNF))-(Nao*exp(-(Vm-50.0)/RToNF)))*d_gate*f_gate*f2_gate;
    double i_si=i_sik+i_sica+i_sina;

    //Na ionic concentrations
    double Nai_prime= -0.001*(i_na+i_na_b+i_fna+i_sina+3*i_p+(n_naca*i_naca)/(n_naca-2.0))/(Vi*F);//mM/ms

    //Ca ionic concentrations
    double i_up=(2*F*Vi/(tau_up*Ca_up_max))*Cai*(Ca_up_max-Ca_up);//nA
    double i_tr=(2*F*Vrel/tau_rep)*p*(Ca_up-Ca_rel);//nA
    double i_rel=(2*F*Vrel/tau_rel)*Ca_rel*(pow(Cai,rCa))/(pow(Cai,rCa)+pow(KmCa,rCa));//nA

    double Ca_up_prime=0.001*(i_up-i_tr)/(2*Vup*F);//mM/ms
    double Ca_rel_prime=0.001*(i_tr-i_rel)/(2*Vrel*F);//mM/ms
    double Cai_prime=-0.001*(i_sica+i_ca_b+i_up-i_rel-2*(i_naca/(n_naca-2.0)))/(2*Vi*F);//mM/ms

    //K ionic concentrations
    double i_mk=(i_k1+i_k+i_fk+i_sik+i_to-2*i_p);//nA
    double Kc_prime=i_mk*0.001/(F*Ve)-pf*(Kc-Kb);//mM/ms
    double Ki_prime=(-i_mk*0.001/(Vi*F));//mM/ms

    //stimulus current
    double i_stim = GetStimulus(time);

    //calculate dV
    double Vm_prime = (-1.0/membrane_C)*(i_f+
                                         i_k+
                                         i_k1+
                                         i_to+
                                         i_na_b+
                                         i_p+
                                         i_naca+
                                         i_ca_b+
                                         i_na+
                                         i_si+
                                         i_stim);//mV/ms
    // do not update voltage if the mSetVoltageDerivativeToZero flag has been set
    if (mSetVoltageDerivativeToZero)
    {
        Vm_prime = 0;
    }

    rDY[0] = Vm_prime;
    rDY[1] = y_gate_prime;
    rDY[2] = s_gate_prime;
    rDY[3] = m_gate_prime;
    rDY[4] = h_gate_prime;
    rDY[5] = d_gate_prime;
    rDY[6] = f_gate_prime;
    rDY[7] = f2_gate_prime;
    rDY[8] = p_prime;
    rDY[9] = Nai_prime;
    rDY[10]= Cai_prime;
    rDY[11] = Kc_prime;
    rDY[12] = Ki_prime;
    rDY[13] = x_gate_prime;
    rDY[14] = Ca_up_prime;
    rDY[15] = Ca_rel_prime;
}


double DiFrancescoNoble1985OdeSystem::GetIIonic()
{

    double Vm = mStateVariables[0];
    double y_gate = mStateVariables[1];
    double s_gate = mStateVariables[2];
    double m_gate = mStateVariables[3];
    double h_gate = mStateVariables[4];
    double d_gate= mStateVariables[5];
    double f_gate = mStateVariables[6];
    double f2_gate = mStateVariables[7];
    //double p = mStateVariables[8];
    double Nai = mStateVariables[9];
    double Cai = mStateVariables[10];
    double Kc = mStateVariables[11];
    double Ki = mStateVariables[12];
    double x_gate = mStateVariables[13];
    //double Ca_up = mStateVariables[14];
    //double Ca_rel = mStateVariables[15];

    /*
     * Compute the DiFrancescoNoble1985OdeSystem model
     */
    double Emh=RToNF*log((Nao+0.12*Kc)/(Nai+0.12*Ki));
    double E_na=RToNF*log(Nao/Nai);
    double E_k=RToNF*log(Kc/Ki);
    double E_ca=0.5*RToNF*log(Cao/Cai);

    //i_na
    double i_na=g_na*pow(m_gate,3)*h_gate*(Vm-Emh);
    //i_f
    double i_fk=(Kc/(Kc+Kmf))*(g_fk*(Vm-E_k));
    double i_fna=(Kc/(Kc+Kmf))*g_fna*(Vm-E_na);
    double i_f=y_gate*(i_fk+i_fna);
    //i_k
    double I_K=i_kmax*(Ki-Kc*exp(-Vm/RToNF))/140.0;
    double i_k=x_gate*I_K;
    //i_k1
    double i_k1=(g_k1*Kc/(Kc+Km1))*(Vm-E_k)/(1.0+exp((Vm-E_k+10)*2.0/RToNF));
    //i_to
    double i_to=((((s_gate*g_to*Cai*(0.2+Kc/(Kmto+Kc)))/(Km_Ca+Cai))*(Vm+10.0))/(1.0-exp(-0.2*(Vm+10.0))))*((Ki*exp((Vm*0.5)/RToNF))-(Kc*exp(-0.5*Vm/RToNF)));
    //i_na_b
    double i_na_b=g_na_b*(Vm-E_na);
    //i_p
    double i_p=((I_P*Kc)/(KmK+Kc))*Nai/(Nai+KmNa);
    //i_naca
    double i_naca=(k_naca*((Cao*(exp((gamma*Vm*(n_naca-2))/(RToNF)))*(pow(Nai,n_naca)))-(Cai*(exp(((gamma-1.0)*(n_naca-2.0)*Vm)/(RToNF)))*(pow(Nao,n_naca)))))/((1.0+(d_naca*((Cai*pow(Nao,n_naca))+(Cao*pow(Nai,n_naca)))))*(1.0+(Cai/0.0069)));
    //i_b_ca
    double i_ca_b=g_ca_b*(Vm-E_ca);
    //i_si
    double i_sica=((4.0*P_si*(Vm-50.0))/(RToNF*(1.0-(exp((-2.0*(Vm-50.0))/RToNF)))))*((Cai*exp(100.0/RToNF))-(Cao*exp(-2.0*(Vm-50.0)/RToNF)))*d_gate*f_gate*f2_gate;
    double i_sik=((0.01*P_si*(Vm-50.0))/(RToNF*(1.0-(exp((-(Vm-50.0))/RToNF)))))*((Ki*exp(50.0/RToNF))-(Kc*exp(-(Vm-50.0)/RToNF)))*d_gate*f_gate*f2_gate;
    double i_sina=((0.01*P_si*(Vm-50.0))/(RToNF*(1.0-(exp((-(Vm-50.0))/RToNF)))))*((Nai*exp(50.0/RToNF))-(Nao*exp(-(Vm-50.0)/RToNF)))*d_gate*f_gate*f2_gate;
    double i_si=i_sik+i_sica+i_sina;

    double i_ionic =   (i_f+
                         i_k+
                         i_k1+
                         i_to+
                         i_na_b+
                         i_p+
                         i_naca+
                         i_ca_b+
                         i_na+
                         i_si); /*this is in nA*/

    assert(!std::isnan(i_ionic));

    /*
     * The return value has to be scaled to match the units required by the mono/bidomain equations.
     * The cell model ionic current is in nano Amps, we require micro Amps/cm^2.
     * The estimate of the cell area is obtained by observing that Cm in the cell model and Cm in the bidomain equation are conceptually the same thing.
     * The Cm in the bidomain equation is expressed in capacitance units per area.
     * An estimate of the cell area is then the ratio of the two values of Cm.
     *
     */

    double i_ionic_in_microA = i_ionic*pow(10,-3);
    double estimated_cell_surface_in_cm_square = 0.075 / HeartConfig::Instance()->GetCapacitance();
    double i_ionic_in_microA_per_cm2=i_ionic_in_microA / estimated_cell_surface_in_cm_square;
    return i_ionic_in_microA_per_cm2;

}

void DiFrancescoNoble1985OdeSystem::VerifyStateVariables()
{
//#ifndef NDEBUG
    const std::vector<double>& rY = rGetStateVariables();

    const double Vm = rY[0];
    const double y_gate = rY[1];
    const double r_gate = rY[2];
    const double m_gate = rY[3];
    const double h_gate = rY[4];
    const double d_gate= rY[5];
    const double f_gate = rY[6];
    const double f2_gate = rY[7];
    const double p = rY[8];
    const double Nai = rY[9];
    const double Cai = rY[10];
    const double Kc = rY[11];
    const double Ki = rY[12];
    const double x_gate = rY[13];
    const double Ca_up = rY[14];
    const double Ca_rel = rY[15];
    #define COVERAGE_IGNORE
    if (!(200>=Vm && Vm>=-200))
    {
        EXCEPTION(DumpState("Vm is really out of range!"));
    }
    if (!(0.0<=y_gate && y_gate<=1.0))
    {
        EXCEPTION(DumpState("y gate has gone out of range. Check model parameters, for example spatial stepsize"));
    }
    if (!(0.0<=r_gate && r_gate<=1.0))
    {
        EXCEPTION(DumpState("r gate has gone out of range. Check model parameters, for example spatial stepsize"));
    }
    if (!(0.0<=m_gate && m_gate<=1.0))
    {
        EXCEPTION(DumpState("m gate  has gone out of range. Check model parameters, for example spatial stepsize"));
    }
    if (!(0.0<=h_gate && h_gate<=1.0))
    {
        EXCEPTION(DumpState("h gate  has gone out of range. Check model parameters, for example spatial stepsize"));
    }
    if (!(0.0<=d_gate && d_gate<=1.0))
    {
        EXCEPTION(DumpState("d gate  has gone out of range. Check model parameters, for example spatial stepsize"));
    }
    if (!(0.0<=f_gate && f_gate<=1.0))
    {
        EXCEPTION(DumpState("f gate for the plateau current has gone out of range. Check model parameters, for example spatial stepsize"));
    }
    if (!(0.0<=f2_gate && f2_gate<=1.0))
    {
        EXCEPTION(DumpState("f2 gate has gone out of range. Check model parameters, for example spatial stepsize"));
    }
    if (!(0.0<=p && p<=1.0))
    {
        EXCEPTION(DumpState("p gate has gone out of range. Check model parameters, for example spatial stepsize"));
    }
    if (!(Nai>0.0))
    {
        EXCEPTION(DumpState("intarcellular Na concentration is negative!"));
    }
    if (!(Cai>0.0))
    {
        EXCEPTION(DumpState("intarcellular Ca concentration is negative!"));
    }
    if (!(Kc>0.0))
    {
        EXCEPTION(DumpState("extracellular K concentration is negative!"));
    }
    if (!(Ki>0.0))
    {
        EXCEPTION(DumpState("intarcellular K concentration is negative!"));
    }
    if (!(0.0<=x_gate && x_gate<=1.0))
    {
        EXCEPTION(DumpState("x gate  has gone out of range. Check model parameters, for example spatial stepsize"));
    }
    if (!(Ca_up>0.0))
    {
        EXCEPTION(DumpState("Ca_up concentration is negative!"));
    }
    if (!(Ca_rel>0.0))
    {
        EXCEPTION(DumpState("Ca_rel concentration is negative!"));
    }
    #undef COVERAGE_IGNORE
//#endif
}

template<>
void OdeSystemInformation<DiFrancescoNoble1985OdeSystem>::Initialise(void)
{
    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("mV");
    this->mInitialConditions.push_back(-87.01);

    this->mVariableNames.push_back("y_gate");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.2);

    this->mVariableNames.push_back("s_gate");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("m_gate");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("h_gate");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.8);

    this->mVariableNames.push_back("d_gate");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.005);

    this->mVariableNames.push_back("f_gate");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("f2_gate");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("p");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("Nai");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(8.0);

    this->mVariableNames.push_back("Cai");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.00005);

    this->mVariableNames.push_back("Kc");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(4.0);

    this->mVariableNames.push_back("Ki");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(140.0);

    this->mVariableNames.push_back("x_gate");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("Ca_up");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(2.0);

    this->mVariableNames.push_back("Ca_rel");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(1.0);

    this->mInitialised = true;
}
