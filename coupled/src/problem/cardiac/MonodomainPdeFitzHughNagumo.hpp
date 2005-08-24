#ifndef _MONODOMAINPDEFITZHUGHNAGUMO_HPP_
#define _MONODOMAINPDEFITZHUGHNAGUMO_HPP_

#include <iostream>
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
#include "FitzHughNagumo1961OdeSystem.hpp"
#include "AbstractLinearParabolicPde.hpp"
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
    
    
    // Default stimulus function
    AbstractStimulusFunction*                mpZeroStimulus;
    
    
    // Stimulus function applied to each node
    std::vector<AbstractStimulusFunction* >  mStimulusAtNode;
    
     // Vector of pointers to an ODE system object for each node
    std::vector<FitzHughNagumo1961OdeSystem*> mOdeSystems;
   
    
    
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
              
        mStimulusAtNode.resize(AbstractCoupledPde<SPACE_DIM>::mNumNodes);
		mOdeSystems.reserve(AbstractCoupledPde<SPACE_DIM>::mNumNodes);
        
        // Initialise as zero stimulus everywhere.
        mpZeroStimulus = new InitialStimulus(0, 0); 
                        
        for (int i=0; i<numNodes; i++)
        {   
            mStimulusAtNode[i] = mpZeroStimulus;
	    mOdeSystems.push_back(new FitzHughNagumo1961OdeSystem(mStimulusAtNode[i]));
        }        
    }


    ~MonodomainPdeFitzHughNagumo(void)
    {
		delete mpZeroStimulus;
		for (int i=0; i<AbstractCoupledPde<SPACE_DIM>::mNumNodes; i++)
		{
	    	delete mOdeSystems[i];
		}
    }
    
    // This should not be called as it is a virtual function, use ComputeLinearSourceTermAtNode instead
    double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
    {
        assert(0);
	return 0.0;
    }
    
    // This should not be called as it is a virtual function, use ComputeNonlinearSourceTermAtNode instead
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
        assert(0);
	return 0.0; 
    }

        
    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x)
    {
        return 10 * MatrixDouble::Identity(SPACE_DIM);
    }
    
    /** double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double voltage)
     * 
     * Main method in this class.
     * First checks to see if the ode set of equations have been solved for
     * in this timestep. If not, it integrates the odes over the timestep,
     * and uses the new results for the gating variables, together with the
     * OLD voltage, to calculate and return the ionic current.
     */
    double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double voltage)
    {
        int index = node.GetIndex();
       
        return AbstractCoupledPde<SPACE_DIM>::solutionCache[index];
   
    }
    
    
    double ComputeLinearSourceTermAtNode(const Node<SPACE_DIM>& node)
    {   
        return 0;
    }
    
    // Capacitance = 1
    double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> x)
    {
        return 1;
    }
    
    
    // Apply same initial conditions to each node in the mesh
    void SetUniversalInitialConditions(odeVariablesType initialConditions)
    {
        for(int i=0; i<AbstractCoupledPde<SPACE_DIM>::mNumNodes; i++)
        {
            AbstractCoupledPde<SPACE_DIM>::mOdeVarsAtNode[i] = initialConditions;
        }
    }
    
    
    // Set given stimulus function to a particular node
    void SetStimulusFunctionAtNode(int nodeIndex, AbstractStimulusFunction* pStimulus)
    {
        mStimulusAtNode[ nodeIndex ] = pStimulus;        
    }
    

    // This function informs the class that the current pde timestep is over,
    //  so the odes are reset as being unsolved.
    void ResetAsUnsolvedOdeSystem()
    {
        AbstractCoupledPde<SPACE_DIM>::mTime += AbstractCoupledPde<SPACE_DIM>::mBigTimeStep;
          
    }
    
    
    odeVariablesType GetOdeVarsAtNode( int index )
    {
        return AbstractCoupledPde<SPACE_DIM>::mOdeVarsAtNode[index];
    }        



    // Calculate the ionic current, using the value of the gating variables
    //  at time t+dt, but using the old voltage at time t. 
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
     	
     	double *currentSolutionArray;
        int ierr = VecGetArray(currentSolution, &currentSolutionArray);
     	
     	for (int index=0; index<AbstractCoupledPde<SPACE_DIM>::mNumNodes; index++)
     	{
    	    FitzHughNagumo1961OdeSystem* pFitzHughNagumoOdeSystem = mOdeSystems[index];
        
            
            // overwrite the voltage with the input value
            AbstractCoupledPde<SPACE_DIM>::mOdeVarsAtNode[index][0] = currentSolutionArray[index];             
            
	    	// solve            
            OdeSolution solution = AbstractCoupledPde<SPACE_DIM>::mpOdeSolver->Solve(
                        pFitzHughNagumoOdeSystem, AbstractCoupledPde<SPACE_DIM>::mTime, 
                        AbstractCoupledPde<SPACE_DIM>::mTime+AbstractCoupledPde<SPACE_DIM>::mBigTimeStep, 
                        AbstractCoupledPde<SPACE_DIM>::mSmallTimeStep, 
                        AbstractCoupledPde<SPACE_DIM>::mOdeVarsAtNode[ index ]);
                    
            // extract solution at end time and save in the store 
            for(int j=0; j < 2; j++)
            {
                AbstractCoupledPde<SPACE_DIM>::mOdeVarsAtNode[ index ][j] = solution.mSolutions[ solution.mSolutions.size()-1 ][j];
            }
            
       		double Itotal = mStimulusAtNode[index]->GetStimulus(AbstractCoupledPde<SPACE_DIM>::mTime+AbstractCoupledPde<SPACE_DIM>::mBigTimeStep) +
	    					GetIIonic( AbstractCoupledPde<SPACE_DIM>::mOdeVarsAtNode[ index ] );
	    	AbstractCoupledPde<SPACE_DIM>::solutionCache[index] = -Itotal;
         }
  	}
};	        
#endif //_MONODOMAINPDEFITZHUGHNAGUMO_HPP_

