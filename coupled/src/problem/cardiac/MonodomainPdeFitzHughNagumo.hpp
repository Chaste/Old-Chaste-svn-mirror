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
class MonodomainPdeFitzHughNagumo : public AbstractLinearParabolicPde<SPACE_DIM>
{
private:
    // friend class TestMonodomainPdeFitzHughNagumo;

    /// timestep used in the ode solvers        
    double mSmallTimeStep;
    
    /// timestep used by the pde solver
    double mBigTimeStep;
    
    AbstractIvpOdeSolver *mpOdeSolver;
    
    /// number of nodes in the mesh 
    int mNumNodes;
    
    // Default stimulus function
    AbstractStimulusFunction*                mpZeroStimulus;
    
    /** mOdeVarsAtNode[i] is a vector of the current values of the
     *  voltage, gating variables, intracellular Calcium concentration at node 
     *  i. The voltage returned by the ode solver is not used later since the pde
     *  solves for the voltage.
     */
    std::vector<odeVariablesType>            mOdeVarsAtNode;
    
    /** Stimulus function applied to each node
     */
    std::vector<AbstractStimulusFunction* >  mStimulusAtNode;
    
    /**
     * Vector of pointers to an ODE system object for each node
     */
    std::vector<FitzHughNagumo1961OdeSystem*> mOdeSystems;

    /** boolean stating whether the gating variables have been solved for at this node
     *  yet
     */
    std::vector<bool>                        mOdeSolvedAtNode;      
    
    double mTime;                  
    
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
				double tStart, double bigTimeStep, double smallTimeStep)
    {
        assert(smallTimeStep < bigTimeStep + 1e-10);
        assert(numNodes > 0);
        
        mNumNodes = numNodes;
        mBigTimeStep = bigTimeStep;
        mpOdeSolver = pOdeSolver;
        mSmallTimeStep = smallTimeStep;
     
        mTime = tStart;
        
        mOdeVarsAtNode.resize(mNumNodes);
        mOdeSolvedAtNode.resize(mNumNodes);
        mStimulusAtNode.resize(mNumNodes);
	mOdeSystems.reserve(mNumNodes);
        
        // Initialise as zero stimulus everywhere.
        mpZeroStimulus = new InitialStimulus(0, 0); 
                        
        for (int i=0; i<numNodes; i++)
        {   
            mOdeSolvedAtNode[i] = false;
            mStimulusAtNode[i] = mpZeroStimulus;
	    mOdeSystems.push_back(new FitzHughNagumo1961OdeSystem(mStimulusAtNode[i]));
        }        
    }

    /**
     * Destructor to free allocated memory
     */
    ~MonodomainPdeFitzHughNagumo(void)
    {
	delete mpZeroStimulus;
	for (int i=0; i<mNumNodes; i++)
	{
	    delete mOdeSystems[i];
	}
    }
    
    /** This should not be called, use ComputeLinearSourceTermAtNode instead
     */
    double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
    {
        assert(0);
	return 0.0; // Avoid compiler warning
    }
    
    /** This should not be called, use ComputeNonlinearSourceTermAtNode instead
     */
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
        assert(0);
	return 0.0; // Avoid compiler warning
    }

        
    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x)
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
    double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double voltage)
    {
        int index = node.GetIndex();
        if( !mOdeSolvedAtNode[ index ] )
        {
            FitzHughNagumo1961OdeSystem* pFitzHughNagumoOdeSystem = mOdeSystems[index];
            
            // overwrite the voltage with the input value
            mOdeVarsAtNode[index][0] = voltage; 
            
            if (0) // fabs(mTime+mBigTimeStep - 0.5) < 1e-4 )
            {
                std::cout << "\n\n--------before-------\n\n";
                std::cout << "t = " << mTime+mBigTimeStep << "\n";
                std::cout << "index = " << index << "\n";
                std::cout << "stim = " << mStimulusAtNode[index]->GetStimulus(mTime) << "\n";
                for(int j = 0;j<2; j++)
                {
                    std::cout << mOdeVarsAtNode[index][j] << "\n";
                }
            }
            
	    // solve            
            OdeSolution solution = mpOdeSolver->Solve(pFitzHughNagumoOdeSystem, mTime, mTime+mBigTimeStep, mSmallTimeStep, mOdeVarsAtNode[ index ]);
            
            if (0) // fabs(mTime+mBigTimeStep - 0.5) < 1e-4 )
            {
                std::cout << "\n\n--------after-------\n";
                std::cout << " t = " << mTime+mBigTimeStep << "\n";
                std::cout << " index = " << index << "\n";
                std::cout << " stim = " << mStimulusAtNode[index]->GetStimulus(mTime) << "\n";
                for(int j=0;j<2; j++)
                {
//                     std::cout << solution.mSolutions[j][0] << "\n";
                    std::cout << " " << solution.mSolutions[ solution.mSolutions.size()-1 ][j] << "\n";
                }
            }
                    
            // extract solution at end time and save in the store 
            //std::cout << "\n\n--------really after-------\n";
            for(int j=0; j < 2; j++)
            {
                mOdeVarsAtNode[ index ][j] = solution.mSolutions[ solution.mSolutions.size()-1 ][j];
                //std::cout << mOdeVarsAtNode[index][j] << "\n";
            }
            mOdeSolvedAtNode[ index ] = true;
           
        }
        
        double Itotal = mStimulusAtNode[index]->GetStimulus(mTime+mBigTimeStep) +
	    GetIIonic( mOdeVarsAtNode[ index ], voltage );
        
        return -Itotal;
    }
    
    
    double ComputeLinearSourceTermAtNode(const Node<SPACE_DIM>& node)
    {   
        return 0;
    }
    
    /** Capacitance = 1
     */
    double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> x)
    {
        return 1;
    }
    
    
    /** Apply same initial conditions to each node in the mesh
     */
    void SetUniversalInitialConditions(odeVariablesType initialConditions)
    {
        for(int i=0; i<mNumNodes; i++)
        {
            mOdeVarsAtNode[i] = initialConditions;
        }
    }
    
    
    /** Set given stimulus function to a particular node
     */
    void SetStimulusFunctionAtNode(int nodeIndex, AbstractStimulusFunction* pStimulus)
    {
        mStimulusAtNode[ nodeIndex ] = pStimulus;        
    }
    

    /** This function informs the class that the current pde timestep is over,
     *  so the odes are reset as being unsolved.
     */
    void ResetAsUnsolvedOdeSystem()
    {
        mTime += mBigTimeStep;
        
        for(int i=0; i<mNumNodes; i++)
        {
            mOdeSolvedAtNode[i] = false;
        }        
    }
    
    
    odeVariablesType GetOdeVarsAtNode( int index )
    {
        return mOdeVarsAtNode[index];
    }        



    /** Calculate the ionic current, using the value of the gating variables
     *  at time t+dt, but using the old voltage at time t.
     */    
    double GetIIonic(odeVariablesType odeVars, double voltage)
    {
	// \todo ignore the voltage returned by the ode system solver ??
	double membrane_V        =  odeVars[0]; // or use voltage;
	double recovery_variable =  odeVars[1];
     
	// Define some constants
		
	double mAlpha = 0.14;//-0.08; // Typical values between 0.10 and 0.15
   		
	// Define nonlinear source term                  
	double nonlinearSourceTerm = membrane_V*(membrane_V-mAlpha)*(1-membrane_V)
	    - recovery_variable;
 
	return nonlinearSourceTerm;
    }
}; 


#endif //_MONODOMAINPDEFITZHUGHNAGUMO_HPP_

