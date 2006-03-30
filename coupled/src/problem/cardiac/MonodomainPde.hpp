#ifndef _MONODOMAINPDE_HPP_
#define _MONODOMAINPDE_HPP_

#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"
#include "OdeSolution.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MatrixDouble.hpp"
#include "AbstractCoupledPde.hpp"



const double rMyo = 150;                                // myoplasmic resistance, ohm*cm
const double rG = 1.5;                                  // gap junction resistance, ohm*cm^2
const double RADIUS = 0.00011;                          // radius of cell, cm
const double LENGTH = 0.01;                             // length of cell, cm
//const double BETA = 2*(RADIUS+LENGTH)/(RADIUS*LENGTH);  // surface to volume ratio
const double rA = rMyo + rG / LENGTH;//* BETA;
//const double DIFFUSION_CONST = 0.5*RADIUS/(2*rA);
//const double DIFFUSION_CONST = 0.0;
//const double DIFFUSION_CONST = 0.0005;

//memfem:
//const double DIFFUSION_CONST = 0.000019;
const double BETA = 0.00014;



/**
 * MonodomainPde class.
 * 
 * A monodomain PDE which deals with some single cell model (e.g. Luo-Rudy) 
 * 
 * Monodmain equation is of the form:
 * c(x) du/dt = a/(2*Rm) *Grad.(Grad(u))  +  LinearSourceTerm(x)  +  NonlinearSourceTerm(x, u)
 * 
 */

typedef std::vector<double> odeVariablesType;



template <int SPACE_DIM>
class MonodomainPde : public AbstractCoupledPde<SPACE_DIM>
{
    private:
        friend class TestMonodomainPde;

        AbstractStimulusFunction* mpZeroStimulus;
        double mDiffusionCoefficient;

        /** The ODE system at each cell. Distributed. */
        std::vector<LuoRudyIModel1991OdeSystem*> mOdeSystemsDistributed;

    public:
    
    //Constructor
    MonodomainPde(int numNodes, AbstractIvpOdeSolver *pOdeSolver, 
                 double tStart, double bigTimeStep, double smallTimeStep) :
    AbstractCoupledPde<SPACE_DIM>(numNodes, pOdeSolver, 
                  tStart,  bigTimeStep,  smallTimeStep)          
    {
        // Initialise as zero stimulus everywhere.
        mpZeroStimulus = new InitialStimulus(0, 0); 

        // Initialise the diffusion coefficient
        
        mDiffusionCoefficient = 0.0005;

        // Create all our ODE systems, with a zero stimulus
        int lo=this->mOwnershipRangeLo;
        int hi=this->mOwnershipRangeHi;
        mOdeSystemsDistributed.reserve(hi-lo);
        for (int local_index=0; local_index < hi-lo; local_index++)
        {
            mOdeSystemsDistributed.push_back(new LuoRudyIModel1991OdeSystem(mpZeroStimulus));
        }
        // Set default initial conditions everywhere
        LuoRudyIModel1991OdeSystem ode_system(mpZeroStimulus);
        this->SetUniversalInitialConditions(ode_system.GetInitialConditions());
    }

    ~MonodomainPde(void)
    {
        delete mpZeroStimulus;
        
        // Free ODE systems
        for (unsigned i=0; i<mOdeSystemsDistributed.size(); i++)
        {
            delete mOdeSystemsDistributed[i];
        }
    }
    
    /**
     * Set the diffusion coefficient
     */
    void SetDiffusionCoefficient(const double& rDiffusionCoefficient)
    {
        mDiffusionCoefficient = rDiffusionCoefficient;
    }
    
    /**
     * This should not be called; use 
     * ComputeLinearSourceTermAtNode instead
     */
    double ComputeLinearSourceTerm(Point<SPACE_DIM> )
    {
        assert(0);
	    return 0.0;
    }
    
    /**
     * This should not be called; use 
     * ComputeNonlinearSourceTermAtNode instead
     */
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> , double )
    {
        assert(0);
	    return 0.0;
    }

        
    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> )
    {
        return  mDiffusionCoefficient * MatrixDouble::Identity(SPACE_DIM);
    }
    
    /** ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double voltage)
     * 
     *  Main function is this class:
     *  computeNonlinearSourceTerm first checks to see if the ode set of equations have been
     *  solved for in this timestep. If not, it integrates the odes over the timestep, and uses
     *  the new results for the gating variables, together with the OLD voltage, to calculate and
     *  return the ionic current.
     */
    double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double )
    {
        int index = node.GetIndex();
        return this->solutionCacheReplicated[index];
    }
    
    
    double ComputeLinearSourceTermAtNode(const Node<SPACE_DIM>& )
    {   
        return 0;
    }
    
    // Capacitance = 1
    double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> )
    {
        return 1;
    }
    
    
    /**
     * Set given stimulus function at a particular node.
     * 
     * @param nodeIndex  Global index specifying the node to set the stimulus at.
     * @param pStimulus  Pointer to the stimulus object to use.
     */
    void SetStimulusFunctionAtNode(int nodeGlobalIndex, AbstractStimulusFunction* pStimulus)
    {
        if (nodeGlobalIndex >= this->mOwnershipRangeLo && nodeGlobalIndex < this->mOwnershipRangeHi)
        {
            int local_index = nodeGlobalIndex - this->mOwnershipRangeLo;
            mOdeSystemsDistributed[local_index]->SetStimulusFunction(pStimulus);
        }
    }
    

    /**
     * This function informs the class that the current pde timestep is over,
     * so time is advanced.
     */
    void ResetAsUnsolvedOdeSystem()
    {
        this->mTime += this->mBigTimeStep;
    }
    

    /**
     * Calculate the ionic current, using the value of the gating variables
     * at time t+dt, but using the old voltage at time t
     * 
     * \todo Add a method to the ODE system object to retrieve the last
     * calculated value for this?
     */
    double GetIIonic(odeVariablesType odeVars)
    {
       double fast_sodium_current_h_gate_h              = odeVars[0];
       double fast_sodium_current_j_gate_j              = odeVars[1];
       double fast_sodium_current_m_gate_m              = odeVars[2];
       double intracellular_calcium_concentration_Cai   = odeVars[3];

       double membrane_V                                = odeVars[4];

       double slow_inward_current_d_gate_d              = odeVars[5];
       double slow_inward_current_f_gate_f              = odeVars[6];
       double time_dependent_potassium_current_X_gate_X = odeVars[7];
  
   
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
    
       /*
        * double fast_sodium_current_h_gate_h_prime = fast_sodium_current_h_gate_alpha_h*(1.0-fast_sodium_current_h_gate_h)-fast_sodium_current_h_gate_beta_h*fast_sodium_current_h_gate_h;
        */
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
    
       /*
        * double fast_sodium_current_j_gate_j_prime = fast_sodium_current_j_gate_alpha_j*(1.0-fast_sodium_current_j_gate_j)-fast_sodium_current_j_gate_beta_j*fast_sodium_current_j_gate_j;
        * double fast_sodium_current_m_gate_alpha_m = 0.32*(membrane_V+47.13)/(1.0-exp(-0.1*(membrane_V+47.13)));
        * double fast_sodium_current_m_gate_beta_m = 0.08*exp(-membrane_V/11.0);
        * double fast_sodium_current_m_gate_m_prime = fast_sodium_current_m_gate_alpha_m*(1.0-fast_sodium_current_m_gate_m)-fast_sodium_current_m_gate_beta_m*fast_sodium_current_m_gate_m;
       */
       double fast_sodium_current_i_Na = fast_sodium_current_g_Na*pow(fast_sodium_current_m_gate_m, 3.0)*fast_sodium_current_h_gate_h*fast_sodium_current_j_gate_j*(membrane_V-fast_sodium_current_E_Na);
       /* 
        * double slow_inward_current_d_gate_alpha_d = 0.095*exp(-0.01*(membrane_V-5.0))/(1.0+exp(-0.072*(membrane_V-5.0)));
        * double slow_inward_current_d_gate_beta_d = 0.07*exp(-0.017*(membrane_V+44.0))/(1.0+exp(0.05*(membrane_V+44.0)));
        * double slow_inward_current_d_gate_d_prime = slow_inward_current_d_gate_alpha_d*(1.0-slow_inward_current_d_gate_d)-slow_inward_current_d_gate_beta_d*slow_inward_current_d_gate_d;
        * double slow_inward_current_f_gate_alpha_f = 0.012*exp(-0.008*(membrane_V+28.0))/(1.0+exp(0.15*(membrane_V+28.0)));
        * double slow_inward_current_f_gate_beta_f = 0.0065*exp(-0.02*(membrane_V+30.0))/(1.0+exp(-0.2*(membrane_V+30.0)));
        * double slow_inward_current_f_gate_f_prime = slow_inward_current_f_gate_alpha_f*(1.0-slow_inward_current_f_gate_f)-slow_inward_current_f_gate_beta_f*slow_inward_current_f_gate_f;
        */
       double slow_inward_current_E_si = 7.7-13.0287*log(intracellular_calcium_concentration_Cai);
       double slow_inward_current_i_si = 0.09*slow_inward_current_d_gate_d*slow_inward_current_f_gate_f*(membrane_V-slow_inward_current_E_si);
       /*
        * double intracellular_calcium_concentration_Cai_prime = -1e-4*slow_inward_current_i_si+0.07*(1e-4-intracellular_calcium_concentration_Cai);
       */
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
       /*
        * double time_dependent_potassium_current_X_gate_alpha_X = 0.0005*exp(0.083*(membrane_V+50.0))/(1.0+exp(0.057*(membrane_V+50.0)));
        * double time_dependent_potassium_current_X_gate_beta_X = 0.0013*exp(-0.06*(membrane_V+20.0))/(1.0+exp(-0.04*(membrane_V+20.0)));
        * double time_dependent_potassium_current_X_gate_X_prime = time_dependent_potassium_current_X_gate_alpha_X*(1.0-time_dependent_potassium_current_X_gate_X)-time_dependent_potassium_current_X_gate_beta_X*time_dependent_potassium_current_X_gate_X;
        */
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
    
    virtual void PrepareForAssembleSystem(Vec currentSolution)
    {
        AbstractCoupledPde<SPACE_DIM>::PrepareForAssembleSystem(currentSolution);
        //std::cout<<"MonodomainPde::PrepareForAssembleSystem\n";

        double *p_current_solution;
        VecGetArray(currentSolution, &p_current_solution);
        int lo=this->mOwnershipRangeLo;
        int hi=this->mOwnershipRangeHi;
        double time=this->mTime;
        double small_time_step=this->mSmallTimeStep;
        double big_time_step=this->mBigTimeStep;
        
    	for (int local_index=0; local_index < hi-lo; local_index++)
    	{
            int global_index = local_index + lo;
    		LuoRudyIModel1991OdeSystem* pLr91OdeSystem = mOdeSystemsDistributed[local_index];
     		
            // overwrite the voltage with the input value
            this->mOdeVarsAtNode[local_index][4] = p_current_solution[local_index]; 
            
            // solve            
            OdeSolution solution =
                this->mpOdeSolver->Solve(pLr91OdeSystem,
                                         time, 
                                         time + big_time_step,
                                         small_time_step,
                                         this->mOdeVarsAtNode[ local_index ]);
                                         
             // extract solution at end time and save in the store 
            this->mOdeVarsAtNode[ local_index ] = solution.mSolutions[ solution.mSolutions.size()-1 ];

            // this tests variables are in the correct range (and maybe resets some if they are)
            mOdeSystemsDistributed[local_index]->VerifyVariables( this->mOdeVarsAtNode[ local_index ] );                          
  
            double Itotal = pLr91OdeSystem->GetStimulus(time + big_time_step)
                            + GetIIonic( this->mOdeVarsAtNode[ local_index ]);
        
		    this->solutionCacheReplicated[global_index] = - Itotal;
        }
        
        AbstractCoupledPde<SPACE_DIM>::ReplicateSolutionCache();
     }
};

#endif //_MONODOMAINPDE_HPP_
