#ifndef ABSTRACTCARDIACPDE_HPP_
#define ABSTRACTCARDIACPDE_HPP_

#include <iostream>
#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "ReplicatableVector.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "AbstractCardiacCell.hpp"
#include "VectorPortion.hpp"
#include "DistributedVector.hpp"

/**
 *  Pde containing common functionality to mono and bidomain pdes.
 */


//// OLD NOTE: read this if AbstractPde is brought back
// IMPORTANT NOTE: the inheritance of AbstractPde has to be 'virtual'
// ie "class AbstractCardiacPde : public virtual AbstractPde"
// because AbstractPde will be the top class in a 'dreaded diamond':
//      A
//     / \     A = AbstractPde, B = AbstractCardiac, C = AbtractLinearParabolic (etc)
//    B   C    D = MonodomainPde
//     \ /
//      D
//
// B and C must use virtual inheritence of A in order for D to only contain 1 instance
// of the member variables in A

template <unsigned SPACE_DIM>
class AbstractCardiacPde
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
    const unsigned mStride;
    
    // number of nodes in the mesh
    unsigned mNumNodes;
    
    // Lowest value of index that this part of the global object stores
    unsigned mOwnershipRangeLo;
    
    // One more than the local highest index
    unsigned mOwnershipRangeHi;
    
    
public:
    AbstractCardiacPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, const unsigned stride=1)
            :  mStride(stride)
    {
        mNumNodes = pCellFactory->GetNumberOfCells();
        
        // Reference: Trayanova (2002 - "Look inside the heart")
        mSurfaceAreaToVolumeRatio = 1400;            // 1/cm
        mCapacitance = 1.0;                          // uF/cm^2
        double const_intra_conductivity = 1.75;      // mS/cm (Averaged)
        
        // Old parameter values used:
        //mSurfaceAreaToVolumeRatio = 1;
        //mCapacitance = 1;
        //double const_intra_conductivity = 0.0005;
        
        mIntracellularConductivityTensor = const_intra_conductivity
                                           * identity_matrix<double>(SPACE_DIM);
        
        // Create a temporary PETSc vector and use the ownership range of
        // the PETSc vector to size our C++ vectors
        DistributedVector::SetProblemSize(mNumNodes);
        this->mOwnershipRangeLo = DistributedVector::Begin().Global;
        this->mOwnershipRangeHi = DistributedVector::End().Global;
 
        mCellsDistributed.resize(DistributedVector::End().Global - DistributedVector::Begin().Global);
        
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            mCellsDistributed[index.Local] = pCellFactory->CreateCardiacCellForNode(index.Global);
        }
        pCellFactory->FinaliseCellCreation(&mCellsDistributed,
                                           DistributedVector::Begin().Global,
                                           DistributedVector::End().Global);
        
        mIionicCacheReplicated.resize( pCellFactory->GetNumberOfCells() );
        mIntracellularStimulusCacheReplicated.resize( pCellFactory->GetNumberOfCells() );
    }
    
    
    virtual ~AbstractCardiacPde()
    {
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            delete mCellsDistributed[index.Local];
        }
    }
    
    
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
        if (!(DistributedVector::Begin().Global <= globalIndex && globalIndex < DistributedVector::End().Global))
        {
            #define COVERAGE_IGNORE
            std::cout << "i " << globalIndex << " lo " << DistributedVector::Begin().Global <<
            " hi " << DistributedVector::End().Global << std::endl;
            #undef COVERAGE_IGNORE
        }
#endif
        assert(DistributedVector::Begin().Global <= globalIndex && globalIndex < DistributedVector::End().Global);
        return mCellsDistributed[globalIndex-DistributedVector::Begin().Global];
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
        DistributedVector striped_distributed(currentSolution);
        DistributedVector::Stripe voltage_solution(striped_distributed,0);

        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            mCellsDistributed[index.Local]->SetVoltage( voltage_solution[index] );
            try
            {
                // solve
                // Note: Voltage should not be updated. GetIIonic will be called later
                // and needs the old voltage. The voltage will be updated from the pde.
                mCellsDistributed[index.Local]->ComputeExceptVoltage(currentTime, nextTime);
            }
            catch (Exception &e)
            {
                ReplicateException(true);
                throw e;
            }
            // update the Iionic and stimulus caches
            UpdateCaches(index.Global, index.Local, nextTime);
        }
        striped_distributed.Restore();
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
        unsigned lo=DistributedVector::Begin().Global;
        unsigned hi=DistributedVector::End().Global;
        
        mIionicCacheReplicated.Replicate(lo, hi);
        mIntracellularStimulusCacheReplicated.Replicate(lo, hi);
    }
    
    void ReplicateException(bool flag)
    {
        unsigned my_error = (unsigned) flag;
        unsigned anyones_error;
        MPI_Allreduce(&my_error, &anyones_error, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
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

#endif /*ABSTRACTCARDIACPDE_HPP_*/

