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
    double mSurfaceAreaToVolumeRatio;
    double mCapacitance;

    c_matrix<double, SPACE_DIM, SPACE_DIM> mIntracellularConductivityTensor;
    c_matrix<double, SPACE_DIM, SPACE_DIM> mExtracellularConductivityTensor;


public:    
    //Constructor     
    BidomainPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, double tStart, double pdeTimeStep) 
       :  AbstractCardiacPde<SPACE_DIM>(pCellFactory, tStart, pdeTimeStep)          
    {
        //\todo: pick good default values;
        mSurfaceAreaToVolumeRatio = 1;
        mCapacitance = 1; 
        
        double const_intra_conductivity = 0.0005;
        double const_extra_conductivity = 0.0005;
        
        //\todo: this loop needed? if so, is there a Zero() type call for ublas
        for(int i=0;i<SPACE_DIM;i++)
        {
            for(int j=0; j<SPACE_DIM; j++)
            {
                mIntracellularConductivityTensor(i,j) = 0;
                mExtracellularConductivityTensor(i,j) = 0;
            }
        }

        for(int i=0;i<SPACE_DIM;i++)
        {
            mIntracellularConductivityTensor(i,i) = const_intra_conductivity;
            mExtracellularConductivityTensor(i,i) = const_extra_conductivity;
        }
    }

    ~BidomainPde()
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

    void SetExtracellularConductivityTensor(c_matrix<double, SPACE_DIM, SPACE_DIM> extracellularConductivity)
    {
        mExtracellularConductivityTensor = extracellularConductivity;
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

    c_matrix<double, SPACE_DIM, SPACE_DIM> GetExtracellularConductivityTensor()
    {
        return mExtracellularConductivityTensor;
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
