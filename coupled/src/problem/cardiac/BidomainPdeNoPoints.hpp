#ifndef BIDOMAINPDE_HPP_
#define BIDOMAINPDE_HPP_

#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
//#include "AbstractCardiacPde.hpp"

#include <iostream>
#include "AbstractStimulusFunction.hpp"
#include "ReplicatableVector.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "AbstractCardiacCell.hpp"

/**
 * BidomainPde class.
 *
 * The bidmain equation is of the form:
 *
 * A_m ( C_m d(V_m)/dt + I_ion ) = div ( sigma_i grad( V_m + phi_e ) ) + I_intra_stim
 *   and
 * div ( (sigma_i + sigma_e) grad phi_e    +   sigma_i (grad V_m) )  + I_extra_stim
 *
 * where V_m is the trans-membrane potential = phi_i - phi_e            (mV),
 *       phi_i is the intracellular potential                           (mV),
 *       phi_e is the intracellular potential                           (mV),
 * and   A_m is the surface area to volume ratio of the cell membrane   (1/cm),
 *       C_m is the membrane capacitance                                (uF/cm^2),
 *       sigma_i is the intracellular conductivity tensor               (mS/cm),
 *       sigma_e is the intracellular conductivity tensor               (mS/cm),
 * and   I_ion is the ionic current                                     (uA/cm^2),
 *       I_intra_stim is the internal stimulus                          (uA/cm^3),
 *       I_extra_stim is the external stimulus (a shock)                (uA/cm^3).
 */
template <int SPACE_DIM>
class BidomainPde //: public AbstractCardiacPde<SPACE_DIM>
{
protected:
    //
    // From AbstractCardiacPde
    //
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
    const unsigned mStride;
    
    // number of nodes in the mesh
    unsigned mNumNodes;
    
    // Lowest value of index that this part of the global object stores
    unsigned mOwnershipRangeLo;
    
    // One more than the local highest index
    unsigned mOwnershipRangeHi;
    
private:
    c_matrix<double, SPACE_DIM, SPACE_DIM> mExtracellularConductivityTensor;
    ReplicatableVector mExtracellularStimulusCacheReplicated;
    
public:
    //Constructor
    BidomainPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            :  mStride(2)
    {
        //
        // From AbstractCardiacPde
        //
        mNumNodes = pCellFactory->GetNumberOfCells();
        
        // Create a temporary PETSc vector and use the ownership range of
        // the PETSc vector to size our C++ vectors
        Vec tempVec;
        VecCreate(PETSC_COMM_WORLD, &tempVec);
        VecSetSizes(tempVec, PETSC_DECIDE, mNumNodes);
        VecSetFromOptions(tempVec);
        PetscInt temp_lo, temp_hi;
        VecGetOwnershipRange(tempVec, &temp_lo, &temp_hi);
        mOwnershipRangeLo=(unsigned) temp_lo;
        mOwnershipRangeHi=(unsigned) temp_hi;
        VecDestroy(tempVec); // vector no longer needed
        
        
        // Reference: Trayanova (2002 - "Look inside the heart")
        mSurfaceAreaToVolumeRatio = 1400;            // 1/cm
        mCapacitance = 1.0;                          // uF/cm^2
        double const_intra_conductivity = 1.75;      // mS/cm (Averaged)
        
        // Old parameter values used:
        //mSurfaceAreaToVolumeRatio = 1;
        //mCapacitance = 1;
        //double const_intra_conductivity = 0.0005;
        
        mIntracellularConductivityTensor.clear();
        
        for (int i=0;i<SPACE_DIM;i++)
        {
            mIntracellularConductivityTensor(i,i) = const_intra_conductivity;
        }
        
        unsigned lo=this->mOwnershipRangeLo;
        unsigned hi=this->mOwnershipRangeHi;
        
        mCellsDistributed.resize(hi-lo);
        
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            mCellsDistributed[local_index] = pCellFactory->CreateCardiacCellForNode(global_index);
        }
        pCellFactory->FinaliseCellCreation(&mCellsDistributed, lo, hi);
        
        
        mIionicCacheReplicated.resize( pCellFactory->GetNumberOfCells() );
        mIntracellularStimulusCacheReplicated.resize( pCellFactory->GetNumberOfCells() );
    
        //
        // From us
        //
    
        /**
          *  Parameters used in mono and bidomain simulations
          *  UNITS: surface area to volume ratio: 1/cm,
          *         capacitance                 : uF/cm^2,
          *         conductivity                : mS/cm.
          * 
          *  Extracellular conductivity set 7.0 mS/cm (Ref: Trayanova 2002 - "Look inside the heart")
          */
        double const_extra_conductivity = 7.0;
        
        mExtracellularConductivityTensor.clear();
        for (int i=0;i<SPACE_DIM;i++)
        {
            mExtracellularConductivityTensor(i,i) = const_extra_conductivity;
        }
        
        mExtracellularStimulusCacheReplicated.resize( pCellFactory->GetNumberOfCells() );
    }
    
    ~BidomainPde()
    {
        // From AbstractCardiacPde
        unsigned lo=this->mOwnershipRangeLo;
        unsigned hi=this->mOwnershipRangeHi;
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            delete mCellsDistributed[local_index];
        }
    }
    
    void SetExtracellularConductivityTensor(c_matrix<double, SPACE_DIM, SPACE_DIM> extracellularConductivity)
    {
        mExtracellularConductivityTensor = extracellularConductivity;
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> GetExtracellularConductivityTensor()
    {
        return mExtracellularConductivityTensor;
    }
    
    /**
     * The bidomain pde also updates the extracellular stimulus cache
     */
    void UpdateCaches(unsigned globalIndex, unsigned localIndex, double nextTime)
    {
        // First 2 lines from AbstractCardiacPde
        mIionicCacheReplicated[globalIndex] = mCellsDistributed[localIndex]->GetIIonic();
        mIntracellularStimulusCacheReplicated[globalIndex] = mCellsDistributed[localIndex]->GetIntracellularStimulus(nextTime);

        mExtracellularStimulusCacheReplicated[globalIndex] = this->mCellsDistributed[localIndex]->GetExtracellularStimulus(nextTime);
    }
    
    /**
     * The bidomain Pde also replicates the extracellular stimulus cache
     */
    void ReplicateCaches()
    {
        unsigned lo=this->mOwnershipRangeLo;
        unsigned hi=this->mOwnershipRangeHi;
        
        // Next 2 lines from AbstractCardiacPde
        mIionicCacheReplicated.Replicate(lo, hi);
        mIntracellularStimulusCacheReplicated.Replicate(lo, hi);
        
        mExtracellularStimulusCacheReplicated.Replicate(lo, hi);
    }
    
    ReplicatableVector& GetExtracellularStimulusCacheReplicated()
    {
        return mExtracellularStimulusCacheReplicated;
    }
    
    //
    // From AbstractCardiacPde
    //
    void SetSurfaceAreaToVolumeRatio(double surfaceAreaToVolumeRatio)
    {
        if (surfaceAreaToVolumeRatio <= 0)
        {
            EXCEPTION("Surface area to volume ratio should be positive");
        }
        mSurfaceAreaToVolumeRatio = surfaceAreaToVolumeRatio;
    }
    
    void SetCapacitance(double capacitance)
    {
        if (capacitance <= 0)
        {
            EXCEPTION("Capacitance should be positive");
        }
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
    AbstractCardiacCell* GetCardiacCell( unsigned globalIndex )
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
    
    
    /**
     *  SolveCellSystems()
     *  
     *  Integrate the cell ODEs and update ionic current etc for each of the 
     *  cells, between the two times provided.
     * 
     *  NOTE: this used to be PrepareForAssembleSystem, but that method is now
     *  a virtual method in the assemblers not the pdes.
     */
    void SolveCellSystems(Vec currentSolution, double currentTime, double nextTime)
    {
        double *p_current_solution;
        VecGetArray(currentSolution, &p_current_solution);
        unsigned lo=this->mOwnershipRangeLo;
        unsigned hi=this->mOwnershipRangeHi;
        
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
                mCellsDistributed[local_index]->ComputeExceptVoltage(currentTime, nextTime);
            }
            catch (Exception &e)
            {
                ReplicateException(true);
                throw e;
            }
            
            // update the Iionic and stimulus caches
            UpdateCaches(global_index, local_index, nextTime);
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
            EXCEPTION("Another process threw an exception in PrepareForAssembleSystem");
        }
    }
    
    void GetOwnershipRange(unsigned &rLo, unsigned &rHi)
    {
        rLo=mOwnershipRangeLo;
        rHi=mOwnershipRangeHi;
    }
};



#endif /*BIDOMAINPDE_HPP_*/
