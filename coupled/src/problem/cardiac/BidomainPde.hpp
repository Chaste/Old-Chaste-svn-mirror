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
public:    
    //Constructor     
    BidomainPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, double tStart, double pdeTimeStep) 
       :  AbstractCardiacPde<SPACE_DIM>(pCellFactory, tStart, pdeTimeStep)          
    {
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
    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> )
    {
        assert(0);
        return 0 * MatrixDouble::Identity(SPACE_DIM);
    }
   
};



#endif /*BIDOMAINPDE_HPP_*/
