#ifndef _MONODOMAINPARABOLICPDE_HPP_
#define _MONODOMAINPARABOLICPDE_HPP_

//#include <cmath>
#include <iostream>
//#include <fstream>
#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"
#include "RegularStimulus.hpp"

#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "AdamsBashforthIvpOdeSolver.hpp"

#include "OdeSolution.hpp"

//#include "ColumnDataWriter.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

#include "AbstractLinearParabolicPde.hpp"


/**
 * MonodomainPde class.
 * 
 * A monodomain PDE which deals with some single cell model (e.g. Luo-Rudy) 
 * 
 * Monodmain equation is of the form:
 * c(x) du/dt = a/(2*Rm) *Grad.(Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)
 * 
 */



typedef std::vector<double> odeVariablesType;

//template <int SPACE_DIM>
class MonodomainPde : public AbstractLinearParabolicPde<1>
{
    private:
    
        friend class TestMonodomainPde;
        
        double mSmallTimeStep;
        double mBigTimeStep;
        AbstractIvpOdeSolver *mpOdeSolver;
        int mNumNodes;
        
        AbstractStimulusFunction*                mpZeroStimulus;

        std::vector<odeVariablesType>            mOdeVarsAtNode;
        std::vector<AbstractStimulusFunction* >  mStimulusAtNode;
        std::vector<bool>                        mOdeSolvedAtNode;      

        double mTime;                  

    public:
    
    //Constructor
    MonodomainPde(int numNodes, double bigTimeStep, AbstractIvpOdeSolver *pOdeSolver, double smallTimeStep)
    {
        assert(smallTimeStep <= bigTimeStep);
        assert(numNodes > 0);
        
        mNumNodes=numNodes;
        mBigTimeStep=bigTimeStep;
        mpOdeSolver=pOdeSolver;
        mSmallTimeStep=smallTimeStep;
     
        mTime = 0;
        
        mOdeVarsAtNode.resize(mNumNodes);
        mOdeSolvedAtNode.resize(mNumNodes);
        mStimulusAtNode.resize(mNumNodes);
        
        mpZeroStimulus = new InitialStimulus(0, 0); 
                        
        for(int i=0; i<numNodes; i++)
        {   
            mOdeSolvedAtNode[i] = false;
            mStimulusAtNode[i] = mpZeroStimulus;
        }        
        
    }
    
    
    double ComputeLinearSourceTerm(Point<1> x)
    {
        // use ComputeLinearSourceTermAtNode instead
        assert(0);
    }
    
    double ComputeNonlinearSourceTerm(Point<1> x, double u)
    {
        // use ComputeNonlinearSourceTermAtNode instead
        assert(0);
    }
    
    MatrixDouble ComputeDiffusionTerm(Point<1> x)
    {
        return MatrixDouble::Identity(1);
    }
    
    double ComputeNonlinearSourceTermAtNode(const Node<1>& node, double voltage)
    {
        int index = node.GetIndex();
        if( !mOdeSolvedAtNode[ index ] )
        {
            LuoRudyIModel1991OdeSystem* pLr91OdeSystem = new LuoRudyIModel1991OdeSystem( mStimulusAtNode[ index ] );
            
            mOdeVarsAtNode[index][4] = voltage;
            OdeSolution solution = mpOdeSolver->Solve(pLr91OdeSystem, mTime, mTime+mBigTimeStep, mSmallTimeStep, mOdeVarsAtNode[ index ]);
            
            mOdeVarsAtNode[ index ] = solution.mSolutions[ solution.mSolutions.size()-1 ];
            mOdeSolvedAtNode[ index ] = true;
        }
        
        return mStimulusAtNode[index]->GetStimulus(mTime+mBigTimeStep) + GetIIonic( mOdeVarsAtNode[ index ], voltage );
    }
    
    
    double ComputeLinearSourceTermAtNode(const Node<1>& node)
    {   
        return 0;
    }
    
    
    double ComputeDuDtCoefficientFunction(Point<1> x)
    {
        return 1;
    }
    
    // MonodomainPde class needs to store gating variables and time
    void SetUniversalInitialConditions(odeVariablesType initialConditions)
    {
        for(int i=0; i<mNumNodes; i++)
        {
            mOdeVarsAtNode[i] = initialConditions;
        }
    }
    
    
    void SetStimulusFunctionAtNode(int nodeIndex, AbstractStimulusFunction* pStimulus)
    {
        mStimulusAtNode[ nodeIndex ] = pStimulus;        
    }
    

    // Reset function : reset bools and set new time
    // void Reset();
    
        


    
    double GetIIonic(odeVariablesType odeVars, double voltage)
    {
       double fast_sodium_current_h_gate_h              = odeVars[0];
       double fast_sodium_current_j_gate_j              = odeVars[1];
       double fast_sodium_current_m_gate_m              = odeVars[2];
       double intracellular_calcium_concentration_Cai   = odeVars[3];
       double membrane_V                                = odeVars[4];  //voltage;
       double slow_inward_current_d_gate_d              = odeVars[5];
       double slow_inward_current_f_gate_f              = odeVars[6];
       double time_dependent_potassium_current_X_gate_X = odeVars[7];
    
//    double membrane_V,
//                    double fast_sodium_current_m_gate_m,
//                    double fast_sodium_current_h_gate_h,
//                    double fast_sodium_current_j_gate_j,
//                    double slow_inward_current_d_gate_d,
//                    double slow_inward_current_f_gate_f,
//                    double time_dependent_potassium_current_X_gate_X,
//                    double intracellular_calcium_concentration_Cai
   
    // Define some constants
                    
       double membrane_F = 96485.0;
       double membrane_R = 8314;
       double membrane_T = 310.0;
       
       double background_current_E_b = -59.87;
       double background_current_g_b = 0.03921;
       
       double fast_sodium_current_g_Na = 23.0;
       double ionic_concentrations_Ki = 145.0;
       double ionic_concentrations_Ko = 5.4;
       double ionic_concentrations_Nai = 18.0;
       double ionic_concentrations_Nao = 140.0;
       
       double fast_sodium_current_E_Na = ((membrane_R * membrane_T) / membrane_F) * 
                                  log(ionic_concentrations_Nao / ionic_concentrations_Nai);
       
       double plateau_potassium_current_g_Kp = 0.0183;
       double time_dependent_potassium_current_PR_NaK = 0.01833;
       
       
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
    
       double fast_sodium_current_h_gate_h_prime = fast_sodium_current_h_gate_alpha_h*(1.0-fast_sodium_current_h_gate_h)-fast_sodium_current_h_gate_beta_h*fast_sodium_current_h_gate_h;
    
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
    
       double fast_sodium_current_j_gate_j_prime = fast_sodium_current_j_gate_alpha_j*(1.0-fast_sodium_current_j_gate_j)-fast_sodium_current_j_gate_beta_j*fast_sodium_current_j_gate_j;
       double fast_sodium_current_m_gate_alpha_m = 0.32*(membrane_V+47.13)/(1.0-exp(-0.1*(membrane_V+47.13)));
       double fast_sodium_current_m_gate_beta_m = 0.08*exp(-membrane_V/11.0);
       double fast_sodium_current_m_gate_m_prime = fast_sodium_current_m_gate_alpha_m*(1.0-fast_sodium_current_m_gate_m)-fast_sodium_current_m_gate_beta_m*fast_sodium_current_m_gate_m;
       double fast_sodium_current_i_Na = fast_sodium_current_g_Na*pow(fast_sodium_current_m_gate_m, 3.0)*fast_sodium_current_h_gate_h*fast_sodium_current_j_gate_j*(membrane_V-fast_sodium_current_E_Na);
     
       double slow_inward_current_d_gate_alpha_d = 0.095*exp(-0.01*(membrane_V-5.0))/(1.0+exp(-0.072*(membrane_V-5.0)));
       double slow_inward_current_d_gate_beta_d = 0.07*exp(-0.017*(membrane_V+44.0))/(1.0+exp(0.05*(membrane_V+44.0)));
       double slow_inward_current_d_gate_d_prime = slow_inward_current_d_gate_alpha_d*(1.0-slow_inward_current_d_gate_d)-slow_inward_current_d_gate_beta_d*slow_inward_current_d_gate_d;
       
       double slow_inward_current_f_gate_alpha_f = 0.012*exp(-0.008*(membrane_V+28.0))/(1.0+exp(0.15*(membrane_V+28.0)));
       double slow_inward_current_f_gate_beta_f = 0.0065*exp(-0.02*(membrane_V+30.0))/(1.0+exp(-0.2*(membrane_V+30.0)));
       double slow_inward_current_f_gate_f_prime = slow_inward_current_f_gate_alpha_f*(1.0-slow_inward_current_f_gate_f)-slow_inward_current_f_gate_beta_f*slow_inward_current_f_gate_f;
       
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
          time_dependent_potassium_current_Xi_gate_Xi = 1.0;
       }
       
       double time_dependent_potassium_current_X_gate_alpha_X = 0.0005*exp(0.083*(membrane_V+50.0))/(1.0+exp(0.057*(membrane_V+50.0)));
       double time_dependent_potassium_current_X_gate_beta_X = 0.0013*exp(-0.06*(membrane_V+20.0))/(1.0+exp(-0.04*(membrane_V+20.0)));
       double time_dependent_potassium_current_X_gate_X_prime = time_dependent_potassium_current_X_gate_alpha_X*(1.0-time_dependent_potassium_current_X_gate_X)-time_dependent_potassium_current_X_gate_beta_X*time_dependent_potassium_current_X_gate_X;
     
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
       return i_ionic;
    }
};

#endif //_MONODOMAINPARABOLICPDE_HPP_
