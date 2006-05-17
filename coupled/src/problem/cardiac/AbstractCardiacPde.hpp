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
protected:
    /** The vector of cells. Distributed. */
    std::vector< AbstractCardiacCell* > mCellsDistributed;

    ReplicatableVector mIionicCacheReplicated;
    ReplicatableVector mIntracellularStimulusCacheReplicated;
 
    const int mStride;
 
    /** The following are currently only used in Bidomain.
     * PLEASE change this comment when it's no longer true.
     */
    double mSurfaceAreaToVolumeRatio;
    double mCapacitance;
    c_matrix<double, SPACE_DIM, SPACE_DIM> mIntracellularConductivityTensor;

public:
    AbstractCardiacPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, double tStart, double pdeTimeStep, const int stride=1)
      :  AbstractCoupledPde<SPACE_DIM>(pCellFactory->GetNumberOfNodes(), tStart, pdeTimeStep),
         mStride(stride)
    {
        /// \todo: pick good default values;
        mSurfaceAreaToVolumeRatio = 1;
        mCapacitance = 1; 
        double const_intra_conductivity=0.0005;
 
        mIntracellularConductivityTensor.clear();
                
        for(int i=0;i<SPACE_DIM;i++)
        {
            mIntracellularConductivityTensor(i,i) = const_intra_conductivity;
        }
  
        int lo=this->mOwnershipRangeLo;
        int hi=this->mOwnershipRangeHi;

        mCellsDistributed.resize(hi-lo);

        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            mCellsDistributed[local_index] = pCellFactory->CreateCardiacCellForNode(global_index);
        }        
        pCellFactory->FinaliseCellCreation(&mCellsDistributed, lo, hi);        


        mIionicCacheReplicated.resize( pCellFactory->GetNumberOfNodes() );
        mIntracellularStimulusCacheReplicated.resize( pCellFactory->GetNumberOfNodes() );
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


    void SetSurfaceAreaToVolumeRatio(double surfaceAreaToVolumeRatio)
    {
        assert(surfaceAreaToVolumeRatio > 0);
        mSurfaceAreaToVolumeRatio = surfaceAreaToVolumeRatio;
    }
     
    void SetCapacitance(double capacitance)
    {
        assert(capacitance > 0);
        mCapacitance = capacitance;
    }     

    void SetIntracellularConductivityTensor(c_matrix<double, SPACE_DIM, SPACE_DIM> intracellularConductivity)
    {
        mIntracellularConductivityTensor = intracellularConductivity;
    } 

    double GetSurfaceAreaToVolumeRatio()
    {
        return mSurfaceAreaToVolumeRatio;
    }

    double GetCapacitance()
    {
        return mCapacitance;
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> GetIntracellularConductivityTensor()
    {
        return mIntracellularConductivityTensor;
    }

    /** 
     *  Get a pointer to a cell, indexed by the global node index. Should only called by the process
     *  owning the cell though.
     */
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

        double big_timestep=this->mBigTimeStep;
        
        for (unsigned global_index=lo; global_index < hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            
            // overwrite the voltage with the input value
            mCellsDistributed[local_index]->SetVoltage( p_current_solution[mStride*local_index] );
            
            // solve            
            mCellsDistributed[local_index]->Compute(time, time+big_timestep);

            // update the Iionic and stimulus caches
            UpdateCaches(global_index, local_index, time+big_timestep);
        }
        
        ReplicateCaches();
    }

    ReplicatableVector& GetIionicCacheReplicated()
    {
        return mIionicCacheReplicated;
    }

    ReplicatableVector& GetIntracellularStimulusCacheReplicated()
    {
        return mIntracellularStimulusCacheReplicated;
    }

 
    /** 
     *  Update the Iionic and intracellular stimulus caches.
     *  The method is overridden in the BidomainPde to also update the
     *  extracellular stimulus.
     */
    virtual void UpdateCaches(unsigned globalIndex, unsigned localIndex, double nextTime)
    {
        mIionicCacheReplicated[globalIndex] = mCellsDistributed[localIndex]->GetIIonic();
        mIntracellularStimulusCacheReplicated[globalIndex] = mCellsDistributed[localIndex]->GetIntracellularStimulus(nextTime);
    }

    /** 
     *  Replicate the Iionic and intracellular stimulus caches.
     *  The method is overridden in the BidomainPde to also replicate the
     *  extracellular stimulus.
     */
    virtual void ReplicateCaches()
    {
        unsigned lo=this->mOwnershipRangeLo;
        unsigned hi=this->mOwnershipRangeHi;
        
        mIionicCacheReplicated.Replicate(lo, hi);
        mIntracellularStimulusCacheReplicated.Replicate(lo, hi);
    }
};
#endif /*ABSTRACTCARDIACPDE_HPP_*/
