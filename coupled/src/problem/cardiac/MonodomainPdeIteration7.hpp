#ifndef MONODOMAINPDEITERATION7_HPP_
#define MONODOMAINPDEITERATION7_HPP_
#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"
#include "OdeSolution.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MatrixDouble.hpp"
#include "AbstractCoupledPdeIteration7.hpp"



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
 * MonodomainPdeIteration7 class.
 * 
 * A monodomain PDE which deals with some single cell model (e.g. Luo-Rudy) 
 * 
 * Monodmain equation is of the form:
 * c(x) du/dt = a/(2*Rm) *Grad.(Grad(u))  +  LinearSourceTerm(x)  +  NonlinearSourceTerm(x, u)
 * 
 */

typedef std::vector<double> odeVariablesType;



template <int SPACE_DIM>
class MonodomainPdeIteration7 : public AbstractCoupledPdeIteration7<SPACE_DIM>
{
    private:
        friend class TestMonodomainPde;

        double mDiffusionCoefficient;

        /** The vector of cells. Distributed. */
        //\todo parallelise cells
        std::vector< AbstractCardiacCell* > mCells;

    public:
    
    //Constructor     
    MonodomainPdeIteration7(std::vector< AbstractCardiacCell* > cells, double tStart, double bigTimeStep) :
        AbstractCoupledPdeIteration7<SPACE_DIM>(cells.size(), tStart,  bigTimeStep)          
     {
        // Initialise the diffusion coefficient
        mDiffusionCoefficient = 0.0005;

        // store cells        
        mCells = cells;
     }

    ~MonodomainPdeIteration7(void)
    {
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
     * This function informs the class that the current pde timestep is over,
     * so time is advanced.
     */
    void ResetAsUnsolvedOdeSystem()
    {
        this->mTime += this->mBigTimeStep;
    }
    


    virtual void PrepareForAssembleSystem(Vec currentSolution)
    {
        AbstractCoupledPdeIteration7 <SPACE_DIM>::PrepareForAssembleSystem(currentSolution);
        //std::cout<<"MonodomainPdeIteration7::PrepareForAssembleSystem\n";

        double *p_current_solution;
        VecGetArray(currentSolution, &p_current_solution);
        int lo=this->mOwnershipRangeLo;
        int hi=this->mOwnershipRangeHi;
        double time=this->mTime;

        double big_time_step=this->mBigTimeStep;
        
        for (int local_index=0; local_index < hi-lo; local_index++)
        {
            int global_index = local_index + lo;
            
//            LuoRudyIModel1991OdeSystem* pLr91OdeSystem = mOdeSystemsDistributed[local_index];
            
            // overwrite the voltage with the input value
//            this->mOdeVarsAtNode[local_index][4] = p_current_solution[local_index]; 
  
            mCells[local_index]->SetVoltage( p_current_solution[local_index] );
            
            // solve            
//            this->mpOdeSolver->Solve(pLr91OdeSystem,
//                                     this->mOdeVarsAtNode[ local_index ],
//                                     time, 
//                                     time + big_time_step,
//                                     small_time_step);

            mCells[local_index]->Compute(time, time+big_time_step);

            // this tests variables are in the correct range (and maybe resets some if they are)
//            mOdeSystemsDistributed[local_index]->VerifyVariables( this->mOdeVarsAtNode[ local_index ] );                          
  
  
  
  /////// NEED TO BRING THIS BACK
  //          mCells[local_index]->VerifyVariables( mCells[local_index]->GetStateVariables() );
  
//            double Itotal = pLr91OdeSystem->GetStimulus(time + big_time_step)
//                            + GetIIonic( this->mOdeVarsAtNode[ local_index ]);
            
            double Itotal =   mCells[local_index]->GetStimulus(time + big_time_step) 
                            + mCells[local_index]->GetIIonic();
          
            this->solutionCacheReplicated[global_index] = - Itotal;
        }
        
        AbstractCoupledPdeIteration7 <SPACE_DIM>::ReplicateSolutionCache();
     }
};

#endif /*MONODOMAINPDEITERATION7_HPP_*/
