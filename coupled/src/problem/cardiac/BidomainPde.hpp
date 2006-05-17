#ifndef BIDOMAINPDE_HPP_
#define BIDOMAINPDE_HPP_

#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "MatrixDouble.hpp"
#include "AbstractCardiacPde.hpp"


/**
 * BidomainPde class.
 * 
 * The bidmain equation is of the form:
 * 
 * \chi ( C_m d(V_m)/dt + I_ion ) = div ( \sigma_i grad( V_m + \phi_e ) + I_i
 * div ( (\sigma_i + \sigma_e) \grad phi_e     +    \sigma_i V_m )  + I_e
 * 
 * where V_m is the trans-membrane potential = \phi_i - \phi_e
 *       \phi_i is the intracellular potential
 *       \phi_e is the intracellular potential
 * and   \chi is the surface area to volume ratio of the cell membrane
 *       C_m is the membrane capacitance
 *       \sigma_i is the intracellular conductivity tensor
 *       \sigma_e is the intracellular conductivity tensor
 * and   I_ion is the ionic current
 *       I_i is the internal stimulus 
 *       I_e is the external stimulus (a shock)
 */

template <int SPACE_DIM>
class BidomainPde : public AbstractCardiacPde<SPACE_DIM>
{
    c_matrix<double, SPACE_DIM, SPACE_DIM> mExtracellularConductivityTensor;
    ReplicatableVector mExtracellularStimulusCacheReplicated;


public:    
    //Constructor     
    BidomainPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, double tStart, double pdeTimeStep) 
       :  AbstractCardiacPde<SPACE_DIM>(pCellFactory, tStart, pdeTimeStep, 2)
    {
        double const_extra_conductivity = 0.0005;

        mExtracellularConductivityTensor.clear();

        for(int i=0;i<SPACE_DIM;i++)
        {
            mExtracellularConductivityTensor(i,i) = const_extra_conductivity;
        }
 
        mExtracellularStimulusCacheReplicated.resize( pCellFactory->GetNumberOfNodes() );
    }

    ~BidomainPde()
    {
    }

    void SetExtracellularConductivityTensor(c_matrix<double, SPACE_DIM, SPACE_DIM> extracellularConductivity)
    {
        mExtracellularConductivityTensor = extracellularConductivity;
    } 

    c_matrix<double, SPACE_DIM, SPACE_DIM> GetExtracellularConductivityTensor()
    {
        return mExtracellularConductivityTensor;
    }

    /** The bidomain Pde also updates the extracellular stimulus cache
     */
    void UpdateCaches(unsigned globalIndex, unsigned localIndex, double nextTime)
    {
        AbstractCardiacPde<SPACE_DIM>::UpdateCaches(globalIndex, localIndex, nextTime);
        mExtracellularStimulusCacheReplicated[globalIndex] = this->mCellsDistributed[localIndex]->GetExtracellularStimulus(nextTime);
    }
  
    /** The bidomain Pde also replicates the extracellular stimulus cache
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

    /**
     * This should not be called, as the bidomain is not of the form
     * of a simple linear parabolic pde
     */
    double ComputeLinearSourceTerm(Point<SPACE_DIM> )
    {
        assert(0);
        return 0.0;
    }
    
    /**
     * This should not be called, as the bidomain is not of the form
     * of a simple linear parabolic pde
     */
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> , double )
    {
        assert(0);
        return 0.0;
    }

    /**
     * This should not be called, as the bidomain is not of the form
     * of a simple linear parabolic pde
     */    
    double ComputeLinearSourceTermAtNode(const Node<SPACE_DIM>& )
    {   
        assert(0);
        return 0;
    }
    
    /**
     * This should not be called, as the bidomain is not of the form
     * of a simple linear parabolic pde
     */
    double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> )
    {
        assert(0);
        return 1;
    }
    
    /**
     * This should not be called, as the bidomain is not of the form
     * of a simple linear parabolic pde
     */
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(Point<SPACE_DIM> )
    {
        assert(0);
        identity_matrix<double> id(SPACE_DIM);
        return 0 * id;
    }
};



#endif /*BIDOMAINPDE_HPP_*/
