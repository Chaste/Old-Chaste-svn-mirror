#ifndef MONODOMAINPDE_HPP_
#define MONODOMAINPDE_HPP_
#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "MatrixDouble.hpp"
#include "AbstractCardiacPde.hpp"

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
 * Monodomain equation is of the form:
 * c(x) du/dt = a/(2*Rm) *Grad.(Grad(u))  +  LinearSourceTerm(x)  +  NonlinearSourceTerm(x, u)
 * 
 */
template <int SPACE_DIM>
class MonodomainPde : public AbstractCardiacPde<SPACE_DIM>
{
private:
    friend class TestMonodomainPde;

    double mDiffusionCoefficient;

public:
    
    //Constructor     
    MonodomainPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, double tStart, double pdeTimeStep) 
       :  AbstractCardiacPde<SPACE_DIM>(pCellFactory, tStart, pdeTimeStep)          
    {
        // Initialise the diffusion coefficient
        mDiffusionCoefficient = 0.0005;
    }

    
    /**
     * Set the diffusion coefficient
     */
    void SetDiffusionCoefficient(const double& rDiffusionCoefficient)
    {
        mDiffusionCoefficient = rDiffusionCoefficient;
    }

    
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
 
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(Point<SPACE_DIM> )
    {
        identity_matrix<double> id(SPACE_DIM);
        return  mDiffusionCoefficient * id;
    }
    
    /** ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double voltage)
     * 
     *  Main function is this class:
     *  computeNonlinearSourceTerm first checks to see if the ode set of equations have been
     *  solved for in this timestep. If not, it integrates the odes over the timestep, and uses
     *  the new results for the gating variables, together with the OLD voltage, to calculate and
     *  return the ionic current.
     */
    double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double )
    {
        int index = node.GetIndex();
//        return this->mSolutionCacheReplicated[index];
        return -this->mIionicCacheReplicated[index] - this->mIntracellularStimulusCacheReplicated[index];

    }
    
    
    double ComputeLinearSourceTermAtNode(const Node<SPACE_DIM>& )
    {   
        return 0;
    }
    
    // Capacitance = 1
    double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> )
    {
        return 1;
    }
    
   
};

#endif /*MONODOMAINPDE_HPP_*/
