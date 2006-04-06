#ifndef _MONODOMAINPDEFITZHUGHNAGUMO_HPP_
#define _MONODOMAINPDEFITZHUGHNAGUMO_HPP_

#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"
#include "OdeSolution.hpp"
#include "FitzHughNagumo1961OdeSystem.hpp"
#include "MatrixDouble.hpp"
#include "AbstractCoupledPde.hpp"



/**
 * MonodomainPde class.
 * 
 * A monodomain PDE which deals with some single cell model (e.g. FHN) 
 * 
 * Monodmain equation is of the form:
 * dv/dt = Grad.(Grad(v))  +  NonlinearSourceTerm(w, v)
 * 
 */

typedef std::vector<double> odeVariablesType;



template <int SPACE_DIM>
class MonodomainPdeFitzHughNagumo : public AbstractCoupledPde<SPACE_DIM>
{
private:
    
    
    /** Default stimulus function
     */
    AbstractStimulusFunction* mpZeroStimulus;
    
     /** Vector of pointers to an ODE system object for each node
      */
    std::vector<FitzHughNagumo1961OdeSystem*> mOdeSystemsDistributed;
   
    
    
public:
    
    /**
     * Constructor. Initialises state vectors, associates a default
     * zero stimulus function with each node, and creates an ODE
     * system for each node.
     *
     * @param numNodes The number of nodes in the mesh
     * @param pOdeSolver Pointer to an ODE solver to use
     * @param tStart Start time for the simulation
     * @param bigTimeStep Time duration to simulate the ODE systems for each time
     *    this is requested
     * @param smallTimeStep Time step the ODE solver should use
     */
    MonodomainPdeFitzHughNagumo(int numNodes, AbstractIvpOdeSolver *pOdeSolver,
				double tStart, double bigTimeStep, double smallTimeStep):
    AbstractCoupledPde<SPACE_DIM>(numNodes, pOdeSolver, 
                  tStart,  bigTimeStep,  smallTimeStep)          
    {
        int lo=this->mOwnershipRangeLo;
        int hi=this->mOwnershipRangeHi;
		mOdeSystemsDistributed.reserve(hi-lo);
        
        // Initialise as zero stimulus everywhere.
        mpZeroStimulus = new InitialStimulus(0, 0); 
                        
        for (int i=0; i<numNodes; i++)
        {
	        mOdeSystemsDistributed.push_back(new FitzHughNagumo1961OdeSystem(pOdeSolver,mpZeroStimulus, smallTimeStep));
        }        
    }


    ~MonodomainPdeFitzHughNagumo(void)
    {
		delete mpZeroStimulus;
		for (unsigned i=0; i<mOdeSystemsDistributed.size(); i++)
		{
	    	delete mOdeSystemsDistributed[i];
		}
    }
    
    // This should not be called as it is a virtual function, use ComputeLinearSourceTermAtNode instead
    double ComputeLinearSourceTerm(Point<SPACE_DIM> )
    {
        assert(0);
	    return 0.0;
    }
    
    // This should not be called as it is a virtual function, use ComputeNonlinearSourceTermAtNode instead
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> , double )
    {
        assert(0);
	    return 0.0; 
    }

        
    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> )
    {
        return 10 * MatrixDouble::Identity(SPACE_DIM);
    }
    
    /**
     * Main method in this class.
     * First checks to see if the ode set of equations have been solved for
     * in this timestep. If not, it integrates the odes over the timestep,
     * and uses the new results for the gating variables, together with the
     * OLD voltage, to calculate and return the ionic current.
     */
    double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double )
    {
        int index = node.GetIndex();
       
        return AbstractCoupledPde<SPACE_DIM>::solutionCacheReplicated[index];
   
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
        AbstractCoupledPde<SPACE_DIM>::mTime += AbstractCoupledPde<SPACE_DIM>::mBigTimeStep;
          
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
	double membrane_V        =  odeVars[0];
	double recovery_variable =  odeVars[1];
     
	// Define some constants
		
	double mAlpha = 0.14;//-0.08; // Typical values between 0.10 and 0.15
   		
	// Define nonlinear source term                  
	double nonlinearSourceTerm = membrane_V*(membrane_V-mAlpha)*(1-membrane_V)
	    - recovery_variable;
 
	return nonlinearSourceTerm;
    }


   
	virtual void PrepareForAssembleSystem(Vec currentSolution)
  	{
  		AbstractCoupledPde<SPACE_DIM>::PrepareForAssembleSystem(currentSolution);
      	
     	double *p_current_solution;
        VecGetArray(currentSolution, &p_current_solution);
     	
        for (int local_index=0; local_index<AbstractCoupledPde<SPACE_DIM>::mOwnershipRangeHi-AbstractCoupledPde<SPACE_DIM>::mOwnershipRangeLo; local_index++)
     	{
            int global_index = local_index + AbstractCoupledPde<SPACE_DIM>::mOwnershipRangeLo;
            FitzHughNagumo1961OdeSystem* pFitzHughNagumoOdeSystem = mOdeSystemsDistributed[local_index];
        
            
            // overwrite the voltage with the input value
            AbstractCoupledPde<SPACE_DIM>::mOdeVarsAtNode[local_index][0] = p_current_solution[local_index];             
            
	    	// solve            
            this->mpOdeSolver->Solve(pFitzHughNagumoOdeSystem,
                                     this->mOdeVarsAtNode[ local_index ],
                                     this->mTime, 
                                     this->mTime+this->mBigTimeStep, 
                                     this->mSmallTimeStep,
                                     this->mSmallTimeStep);
            
       		double Itotal = pFitzHughNagumoOdeSystem->GetStimulus(this->mTime+this->mBigTimeStep) +
	    					GetIIonic( this->mOdeVarsAtNode[ local_index ] );
	    	this->solutionCacheReplicated[global_index] = -Itotal;
         }
         
         AbstractCoupledPde<SPACE_DIM>::ReplicateSolutionCache();
  	}
};	        
#endif //_MONODOMAINPDEFITZHUGHNAGUMO_HPP_

