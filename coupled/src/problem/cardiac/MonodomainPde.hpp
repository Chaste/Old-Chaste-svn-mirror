#ifndef MONODOMAINPDE_HPP_
#define MONODOMAINPDE_HPP_
#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "AbstractCardiacPde.hpp"
#include "AbstractLinearParabolicPde.hpp"

/*
const double rMyo = 150;                                // myoplasmic resistance, ohm*cm
const double rG = 1.5;                                  // gap junction resistance, ohm*cm^2
const double RADIUS = 0.00011;                          // radius of cell, cm
const double LENGTH = 0.01;                             // length of cell, cm
//const double BETA = 2*(RADIUS+LENGTH)/(RADIUS*LENGTH);  // surface to volume ratio
const double rA = rMyo + rG / LENGTH;// BETA;
//const double DIFFUSION_CONST = 0.5*RADIUS/(2*rA);
//const double DIFFUSION_CONST = 0.0;
//const double DIFFUSION_CONST = 0.0005;

//memfem:
//const double DIFFUSION_CONST = 0.000019;
const double BETA = 0.00014;
*/


/**
 * MonodomainPde class.
 *
 * The monodomain equation is of the form:
 * A (C dV/dt + Iionic) +Istim = Div( sigma_i Grad(V) )
 *
 * where A is the surface area to volume ratio         (1/cm),
 *       C is the capacitance                          (uF/cm^2),
 *       sigma_i is the intracellular conductivity     (mS/cm),
 *       I_ionic is the ionic current                  (uA/cm^2),
 *       I_stim is the intracellular stimulus current  (uA/cm^3).
 *
 * Note that default values of A, C and sigma_i are stored in the parent class
 */
template <int SPACE_DIM>
class MonodomainPde : public AbstractCardiacPde<SPACE_DIM>, public AbstractLinearParabolicPde<SPACE_DIM>
{
private:
    friend class TestMonodomainPde;
    
public:

    //Constructor
    MonodomainPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            :  AbstractCardiacPde<SPACE_DIM>(pCellFactory)
    {}
    
    
    //The following are hidden from the coverage test while it is waiting
    //for a re-factor. (Ticket #157)
#define COVERAGE_IGNORE
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
#undef COVERAGE_IGNORE
    
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(Point<SPACE_DIM> )
    {
        return this->mIntracellularConductivityTensor;
    }
    
    
    double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double )
    {
        unsigned index = node.GetIndex();
        return  -(this->mSurfaceAreaToVolumeRatio)*(this->mIionicCacheReplicated[index])
                - this->mIntracellularStimulusCacheReplicated[index];
    }
    
    
    double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> )
    {
        return (this->mSurfaceAreaToVolumeRatio)*(this->mCapacitance);
    }
};

#endif /*MONODOMAINPDE_HPP_*/
