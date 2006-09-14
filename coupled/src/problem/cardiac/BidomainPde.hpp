#ifndef BIDOMAINPDE_HPP_
#define BIDOMAINPDE_HPP_

#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "AbstractCardiacPde.hpp"


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
class BidomainPde : public AbstractCardiacPde<SPACE_DIM>
{
    c_matrix<double, SPACE_DIM, SPACE_DIM> mExtracellularConductivityTensor;
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
    {}
    
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
        AbstractCardiacPde<SPACE_DIM>::UpdateCaches(globalIndex, localIndex, nextTime);
        mExtracellularStimulusCacheReplicated[globalIndex] = this->mCellsDistributed[localIndex]->GetExtracellularStimulus(nextTime);
    }
    
    /**
     * The bidomain Pde also replicates the extracellular stimulus cache
     */
    void ReplicateCaches()
    {
        AbstractCardiacPde<SPACE_DIM>::ReplicateCaches();
        unsigned lo=this->mOwnershipRangeLo;
        unsigned hi=this->mOwnershipRangeHi;
        
        mExtracellularStimulusCacheReplicated.Replicate(lo, hi);
    }
    
    ReplicatableVector& GetExtracellularStimulusCacheReplicated()
    {
        return mExtracellularStimulusCacheReplicated;
    }
};



#endif /*BIDOMAINPDE_HPP_*/
