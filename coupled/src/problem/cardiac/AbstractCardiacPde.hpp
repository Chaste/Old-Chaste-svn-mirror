#ifndef ABSTRACTCARDIACPDE_HPP_
#define ABSTRACTCARDIACPDE_HPP_

#include <iostream>
#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "AbstractCoupledPde.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "AbstractCardiacCell.hpp"

/** 
 *  Pde containing common functionality to mono and bidomain pdes.
 */
template <int SPACE_DIM>
class AbstractCardiacPde : public AbstractCoupledPde<SPACE_DIM>
{
protected:

    /**
     *  Parameters used in mono and bidomain simulations
     *  UNITS: surface area to volume ratio: 1/cm
     *         capacitance                 : uF/cm^2
     *         conductivity                : mS/cm
     *  
     *  which means the units of these should be
     *         Iionic                      : uA/cm^2   
     *         Istimuli                    : uA/cm^3   
     */
    double mSurfaceAreaToVolumeRatio;
    double mCapacitance;
    c_matrix<double, SPACE_DIM, SPACE_DIM> mIntracellularConductivityTensor;

    /** The vector of cells. Distributed. */
    std::vector< AbstractCardiacCell* > mCellsDistributed;

    /** 
     *  Caches containing all the ionic and stimulus currents for each node, 
     *  replicated over all processes
     */
    ReplicatableVector mIionicCacheReplicated;
    ReplicatableVector mIntracellularStimulusCacheReplicated;
 
    /** 
     *  Constant set to 1 in monodomain and 2 in bidomain. Used when accessing
     *  the voltage components in the solution vector (because the solution vector
     *  is of the form (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N), where V_j is 
     *  the voltage at node j and phi_j is the extracellular potential at node j.
     */  
    const int mStride;
 

public:
    AbstractCardiacPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, double pdeTimeStep, const int stride=1)
      :  AbstractCoupledPde<SPACE_DIM>(pCellFactory->GetNumberOfCells(), pdeTimeStep),
         mStride(stride)
    {
        // Reference: Trayanova (2002 - "Look inside the heart")
        mSurfaceAreaToVolumeRatio = 1400;            // 1/cm
        mCapacitance = 1.0;                          // uF/cm^2   
        double const_intra_conductivity = 1.75;      // mS/cm (Averaged)
 
        // Old parameter values used:
        //mSurfaceAreaToVolumeRatio = 1;   
        //mCapacitance = 1;                
        //double const_intra_conductivity = 0.0005;    

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


        mIionicCacheReplicated.resize( pCellFactory->GetNumberOfCells() );
        mIntracellularStimulusCacheReplicated.resize( pCellFactory->GetNumberOfCells() );
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
        #ifndef NDEBUG
            if (!(this->mOwnershipRangeLo <= globalIndex && globalIndex < this->mOwnershipRangeHi)) 
            {
                std::cout << "i " << globalIndex << " lo " << this->mOwnershipRangeLo <<
                    " hi " << this->mOwnershipRangeHi << std::endl;
            }
        #endif
        assert(this->mOwnershipRangeLo <= globalIndex && globalIndex < this->mOwnershipRangeHi);
        return mCellsDistributed[globalIndex-this->mOwnershipRangeLo];
    }
       
  
    virtual void PrepareForAssembleSystem(Vec currentSolution, double time)
    {
        AbstractCoupledPde <SPACE_DIM>::PrepareForAssembleSystem(currentSolution, time);

        double *p_current_solution;
        VecGetArray(currentSolution, &p_current_solution);
        unsigned lo=this->mOwnershipRangeLo;
        unsigned hi=this->mOwnershipRangeHi;
        
//        double time=this->mTime;

        double big_timestep=this->mBigTimeStep;
        
        for (unsigned global_index=lo; global_index < hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            
            // overwrite the voltage with the input value
            mCellsDistributed[local_index]->SetVoltage( p_current_solution[mStride*local_index] );
            
            
            try
            {
                
                // solve
                // Note: Voltage should not be updated. GetIIonic will be called later
                // and needs the old voltage. The voltage will be updated from the pde.            
                mCellsDistributed[local_index]->ComputeExceptVoltage(time, time+big_timestep);
            } catch (Exception &e)
            {
               ReplicateException(true);
               throw e;
            }

            // update the Iionic and stimulus caches
            UpdateCaches(global_index, local_index, time+big_timestep);
        }
        VecRestoreArray(currentSolution, &p_current_solution);
        
        ReplicateException(false);
         
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
     
     void ReplicateException(bool flag)
     {
        int my_error = (int) flag;
        int anyones_error;
        MPI_Allreduce(&my_error, &anyones_error, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
        if (flag)
        {
            // Return control to exception thrower
            return;
        }
        if (anyones_error)
        {
            throw Exception("Another process threw an exception in PrepareForAssembleSystem");
        }
     }
};
#endif /*ABSTRACTCARDIACPDE_HPP_*/
