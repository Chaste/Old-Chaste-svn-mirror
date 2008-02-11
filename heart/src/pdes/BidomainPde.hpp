#ifndef BIDOMAINPDE_HPP_
#define BIDOMAINPDE_HPP_

#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "AbstractCardiacPde.hpp"
#include <boost/numeric/ublas/matrix.hpp>


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
template <unsigned SPACE_DIM>
class BidomainPde : public AbstractCardiacPde<SPACE_DIM>
{

private:
    ElementwiseConductivityTensors<SPACE_DIM> *mpExtracellularConductivityTensors, *mpDefaultExtracellularCondTensors;

    ReplicatableVector mExtracellularStimulusCacheReplicated;
    
public:
    //Constructor
    BidomainPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            :  AbstractCardiacPde<SPACE_DIM>(pCellFactory, 2 /*mStride*/)
    {
    
        /**
          *  Parameters used in mono and bidomain simulations
          *  UNITS: surface area to volume ratio: 1/cm,
          *         capacitance                 : uF/cm^2,
          *         conductivity                : mS/cm.
          * 
          *  (Ref: Trayanova 2002 - "Look inside the heart")
          */
          
          
        // Reference Clerc 1976 
        double long_intra_conductivity = 6.2;      // mS/cm (Averaged)
        double trans_intra_conductivity = 2.4;
        double normal_intra_conductivity = 2.4;
                
        mpExtracellularConductivityTensors = new ElementwiseConductivityTensors<SPACE_DIM>;
        mpExtracellularConductivityTensors->SetConstantConductivities(long_intra_conductivity, trans_intra_conductivity, normal_intra_conductivity);
        mpExtracellularConductivityTensors->Init();
        
        // Keep a copy of the pointer to free it at the end (since mpExtracellularConductivityTensors may be changed from outside)
        mpDefaultExtracellularCondTensors = mpExtracellularConductivityTensors;
                                                   
        
        mExtracellularStimulusCacheReplicated.resize( pCellFactory->GetNumberOfCells() );
    }
    
    ~BidomainPde()
    {
        delete mpDefaultExtracellularCondTensors;
    }
    
    void SetExtracellularConductivityTensors(ElementwiseConductivityTensors<SPACE_DIM>* pIntracellularTensors)
    {        
        mpExtracellularConductivityTensors = pIntracellularTensors;
    }
    
    
    const c_matrix<double, SPACE_DIM, SPACE_DIM>& rGetExtracellularConductivityTensor(unsigned elementIndex)
    {
        return (*mpExtracellularConductivityTensors)[elementIndex];      
    }
    
    /**
     * The bidomain pde also updates the extracellular stimulus cache
     */
    void UpdateCaches(unsigned globalIndex, unsigned localIndex, double nextTime)
    {
        AbstractCardiacPde<SPACE_DIM>::UpdateCaches(globalIndex, localIndex, nextTime);
        mExtracellularStimulusCacheReplicated[globalIndex] = this->mCellsDistributed[localIndex]->GetExtracellularStimulus(nextTime);
    }
    
    /**
     * The bidomain Pde also replicates the extracellular stimulus cache
     */
    void ReplicateCaches()
    {
        AbstractCardiacPde<SPACE_DIM>::ReplicateCaches();
        unsigned lo=DistributedVector::Begin().Global;
        unsigned hi=DistributedVector::End().Global;
        
        mExtracellularStimulusCacheReplicated.Replicate(lo, hi);
    }
    
    ReplicatableVector& rGetExtracellularStimulusCacheReplicated()
    {
        return mExtracellularStimulusCacheReplicated;
    }
};



#endif /*BIDOMAINPDE_HPP_*/
