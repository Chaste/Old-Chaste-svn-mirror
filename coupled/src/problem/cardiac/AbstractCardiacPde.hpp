#ifndef ABSTRACTCARDIACPDE_HPP_
#define ABSTRACTCARDIACPDE_HPP_


#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "MatrixDouble.hpp"
#include "AbstractCoupledPde.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "AbstractCardiacCell.hpp"

template <int SPACE_DIM>
class AbstractCardiacPde : public AbstractCoupledPde<SPACE_DIM>
{
private:
    /** The vector of cells. Distributed. */
    std::vector< AbstractCardiacCell* > mCellsDistributed;


public:
    AbstractCardiacPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, double tStart, double pdeTimeStep) 
      :  AbstractCoupledPde<SPACE_DIM>(pCellFactory->GetNumberOfNodes(), tStart, pdeTimeStep)
    {
        int lo=this->mOwnershipRangeLo;
        int hi=this->mOwnershipRangeHi;
        
        mCellsDistributed.resize(hi-lo);
        
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            mCellsDistributed[local_index] = pCellFactory->CreateCardiacCellForNode(global_index);
        }        
        pCellFactory->FinaliseCellCreation(&mCellsDistributed, lo, hi);        
    }


    ~AbstractCardiacPde()
     {
        int lo=this->mOwnershipRangeLo;
        int hi=this->mOwnershipRangeHi;
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            delete mCellsDistributed[local_index];
        }
     }


    AbstractCardiacCell* GetCardiacCell( int globalIndex )
    {
        if (!(this->mOwnershipRangeLo <= globalIndex && globalIndex < this->mOwnershipRangeHi)) {
            std::cout << "i " << globalIndex << " lo " << this->mOwnershipRangeLo <<
                " hi " << this->mOwnershipRangeHi << std::endl;
        }
        assert(this->mOwnershipRangeLo <= globalIndex && globalIndex < this->mOwnershipRangeHi);
        return mCellsDistributed[globalIndex-this->mOwnershipRangeLo];
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
        AbstractCoupledPde <SPACE_DIM>::PrepareForAssembleSystem(currentSolution);
        //std::cout<<"MonodomainPde::PrepareForAssembleSystem\n";

        double *p_current_solution;
        VecGetArray(currentSolution, &p_current_solution);
        unsigned lo=this->mOwnershipRangeLo;
        unsigned hi=this->mOwnershipRangeHi;
        double time=this->mTime;

        double big_time_step=this->mBigTimeStep;
        
        for (unsigned global_index=lo; global_index < hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            
            // overwrite the voltage with the input value
            mCellsDistributed[local_index]->SetVoltage( p_current_solution[local_index] );
            
            // solve            
            
            mCellsDistributed[local_index]->Compute(time, time+big_time_step);

            double Itotal =   mCellsDistributed[local_index]->GetStimulus(time + big_time_step) 
                            + mCellsDistributed[local_index]->GetIIonic();
          
            this->mSolutionCacheReplicated[global_index] = - Itotal;

        }
        
        this->mSolutionCacheReplicated.Replicate(lo, hi);
     }


};
#endif /*ABSTRACTCARDIACPDE_HPP_*/
