#ifndef EXPONENTIALMATERIALLAW_HPP_
#define EXPONENTIALMATERIALLAW_HPP_

#include "AbstractIsotropicIncompressibleMaterialLaw.hpp"
#include "Exception.hpp"

/**
 *  ExponentialMaterialLaw
 *
 *  An exponential isotropic incompressible hyperelastic material law for finite
 *  elasticity
 *
 *  The law is given by a strain energy function
 *      W(I_1,I_2,I_3) = a exp( b(I_1-3) ) - p/2 C^{-1}
 *  in 3d, or
 *      W(I_1,I_2,I_3) = a exp( b(I_1-2) ) - p/2 C^{-1}
 *  in 2d.
 *
 *  Here I_i are the principal invariants of C, the Lagrangian deformation tensor.
 *  (I1=trace(C), I2=trace(C)^2-trace(C^2), I3=det(C)).

 *  Note: only dimension equals 2 or 3 is permitted.
 */

template<unsigned DIM>
class ExponentialMaterialLaw : public AbstractIsotropicIncompressibleMaterialLaw<DIM>
{
private :
    double mA;
    double mB;
    
public :
    double Get_dW_dI1(double I1, double I2)
    {
        return mA * mB * exp(mB*(I1-DIM));
    }
    double Get_dW_dI2(double I1, double I2)
    {
        assert(DIM==3);
        return 0.0;
    }
    double Get_d2W_dI1(double I1, double I2)
    {
        return mA * mB * mB * exp(mB*(I1-DIM));
    }
    double Get_d2W_dI2(double I1, double I2)
    {
        assert(DIM==3);
        return 0.0;
    }
    double Get_d2W_dI1I2(double I1, double I2)
    {
        assert(DIM==3);
        return 0.0;
    }
    
    double GetA()
    {
        return mA;
    }
    double GetB()
    {
        return mB;
    }
    
public :
    /**
     *  Constructor, Taking in the parameters a and b. a must be positive.
     */
    ExponentialMaterialLaw(double a, double b)
    {
        if (DIM!=2 && DIM !=3)
        {
            EXCEPTION("Can only have 2 or 3d incompressible Mooney-Rivlin laws");
        }
        if (a<0.0)
        {
            EXCEPTION("a must be positive");
        }
        
        mA = a;
        mB = b;
    }
};

#endif /*EXPONENTIALMATERIALLAW_HPP_*/
