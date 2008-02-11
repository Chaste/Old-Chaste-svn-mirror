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
#include "DistributedVector.hpp"
#include "EventHandler.hpp"
#include "PetscTools.hpp"
#include "ElementwiseConductivityTensors.hpp"

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
    //c_matrix<double, SPACE_DIM, SPACE_DIM> mIntracellularConductivityTensor;
    ElementwiseConductivityTensors<SPACE_DIM> *mpIntracellularConductivityTensors, *mpDefaultIntracellularCondTensors;
    
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
    
    
public:
    AbstractCardiacPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, const unsigned stride=1)
            :  mStride(stride)
    {
        mNumNodes = pCellFactory->GetNumberOfCells();
        
        DistributedVector::SetProblemSize(mNumNodes);
                
        // Reference: Trayanova (2002 - "Look inside the heart")
        mSurfaceAreaToVolumeRatio = 1400;            // 1/cm
        mCapacitance = 1.0;                          // uF/cm^2
        
        // Reference Clerc 1976 
        double long_intra_conductivity = 1.75;      // mS/cm (Averaged)
        double trans_intra_conductivity = 0.19;
        double normal_intra_conductivity = 0.19;
        
        // Old parameter values used:
        //mSurfaceAreaToVolumeRatio = 1;
        //mCapacitance = 1;
        //double const_intra_conductivity = 0.0005;
        
        //mIntracellularConductivityTensor = const_intra_conductivity
        //                                   * identity_matrix<double>(SPACE_DIM);

        mpIntracellularConductivityTensors = new ElementwiseConductivityTensors<SPACE_DIM>;
        mpIntracellularConductivityTensors->SetConstantConductivities(long_intra_conductivity, trans_intra_conductivity, normal_intra_conductivity);
        mpIntracellularConductivityTensors->Init();
        
        // Keep a copy of the pointer to free it at the end (since mpIntracellularConductivityTensors may be changed from outside)
        mpDefaultIntracellularCondTensors = mpIntracellularConductivityTensors;
        
        mCellsDistributed.resize(DistributedVector::End().Global-DistributedVector::Begin().Global);
        
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
        
        delete mpDefaultIntracellularCondTensors;
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
     
    void SetIntracellularConductivityTensors(ElementwiseConductivityTensors<SPACE_DIM>* pIntracellularTensors)
    {        
        mpIntracellularConductivityTensors = pIntracellularTensors;
    }
    
    double GetSurfaceAreaToVolumeRatio()
    {
        return mSurfaceAreaToVolumeRatio;
    }
    
    double GetCapacitance()
    {
        return mCapacitance;
    }
    
    const c_matrix<double, SPACE_DIM, SPACE_DIM>& rGetIntracellularConductivityTensor(unsigned elementIndex)
    {
        return (*mpIntracellularConductivityTensors)[elementIndex];        
    }
    
    /**
     *  Get a pointer to a cell, indexed by the global node index. Should only called by the process
     *  owning the cell though.
     */
    AbstractCardiacCell* GetCardiacCell( unsigned globalIndex )
    {   
        assert(DistributedVector::Begin().Global <= globalIndex &&
               globalIndex < DistributedVector::End().Global);
        return mCellsDistributed[globalIndex - DistributedVector::Begin().Global];
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
        EventHandler::BeginEvent(SOLVE_ODES);
        
        DistributedVector dist_solution(currentSolution);
        DistributedVector::Stripe voltage(dist_solution, 0);
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            // overwrite the voltage with the input value
            mCellsDistributed[index.Local]->SetVoltage( voltage[index] );            
            try
            {
                // solve
                // Note: Voltage should not be updated. GetIIonic will be called later
                // and needs the old voltage. The voltage will be updated from the pde.
                mCellsDistributed[index.Local]->ComputeExceptVoltage(currentTime, nextTime);
            }
            catch (Exception &e)
            {
                PetscTools::ReplicateException(true);
                throw e;
            }
            
            // update the Iionic and stimulus caches
            UpdateCaches(index.Global, index.Local, nextTime);
        }
        EventHandler::EndEvent(SOLVE_ODES);

        PetscTools::ReplicateException(false);

        EventHandler::BeginEvent(COMMUNICATION);
        ReplicateCaches();
        EventHandler::EndEvent(COMMUNICATION);
    }
    
    ReplicatableVector& rGetIionicCacheReplicated()
    {
        return mIionicCacheReplicated;
    }
    
    ReplicatableVector& rGetIntracellularStimulusCacheReplicated()
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
        mIionicCacheReplicated.Replicate(DistributedVector::Begin().Global, DistributedVector::End().Global);
        mIntracellularStimulusCacheReplicated.Replicate(DistributedVector::Begin().Global, DistributedVector::End().Global);
    }
};

#endif /*ABSTRACTCARDIACPDE_HPP_*/

