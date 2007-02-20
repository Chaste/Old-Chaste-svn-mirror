#ifndef ABSTRACTISOTROPICINCOMPRESSIBLEMATERIALLAW_HPP_
#define ABSTRACTISOTROPICINCOMPRESSIBLEMATERIALLAW_HPP_

#include "AbstractIncompressibleMaterialLaw.hpp"

/** 
 *  AbstractIsotropicIncompressibleMaterialLaw
 * 
 *  An isotropic incompressible hyperelastic material law for finite elastiticy
 *  
 *  The law is given by a strain energy function W(I1,I2,I3), where I_i are the principal
 *  invariants of C, the Lagrangian deformation tensor. (I1=trace(C), I2=trace(C)^2-trace(C^2),
 *  I3=det(C)). Since it is incompressible, the full strain energy has the form
 *  W^{full} = W(I_1,I_2) - p/2 C^{-1}
 * 
 *  Note: only dimension equals 2 or 3 should be permitted. 
 */
template<int DIM>
class AbstractIsotropicIncompressibleMaterialLaw : public AbstractIncompressibleMaterialLaw<DIM>
{
protected :
    virtual double Get_dW_dI1(double I1, double I2)=0;
    virtual double Get_dW_dI2(double I1, double I2)=0;
    virtual double Get_d2W_dI1(double I1, double I2)=0;
    virtual double Get_d2W_dI2(double I1, double I2)=0;
    virtual double Get_d2W_dI1I2(double I1, double I2)=0;

public :
    /**
     *  Compute the (2nd Piola Kirchoff) stress T and the stress derivative dT/dE 
     *  for a given strain.
     *  
     *  NOTE: the strain E is not expected to be passed in, instead the Lagrangian
     *  deformation tensor C is required (recall, E = 0.5(C-I)
     * 
     *  dT/dE is a fourth-order tensor, where dT/dE[M][N][P][Q] = dT^{MN}/dE_{PQ}
     * 
     *  @param C The Lagrangian deformation tensor (F^T F)
     *  @param invC The inverse of C. Should be computed by the user. (Change this?)
     *  @param pressure the current pressure
     *  @param T the stress will be returned in this parameter
     *  @param dTdE the stress derivative will be returned in this parameter, assuming
     *    the final parameter is true
     *  @param computeDTdE a boolean flag saying whether the stress derivative is 
     *    required or not.
     * 
     *  This is the implemtation for an isotropic material law, so the stress etc is
     *  computed by calling methods returning dW/dI1, dW/dI2 etc.
     */  
    void ComputeStressAndStressDerivative(Tensor<2,DIM>&          C, 
                                          Tensor<2,DIM>&          invC,
                                          double                  pressure,
                                          SymmetricTensor<2,DIM>& T,
                                          double                  dTdE[DIM][DIM][DIM][DIM],
                                          bool                    computeDTdE)
    {
        assert((DIM==2) || (DIM==3)); 
        
        static Tensor<2,DIM> identity;
        static bool first=true;
        if(first)
        {
            for(unsigned i=0; i<DIM; i++)
            {
                for(unsigned j=0; j<DIM; j++)
                {
                    identity[i][j] = i==j ? 1.0 : 0.0;
                }
            }
        }
        first = false;

        
        double I1 = trace(C);

        double I2 = 0.0;
        if(DIM==3)
        {
            I2 =  C[0][1]*C[1][0] + C[1][2]*C[2][1] + C[2][0]*C[0][2]
                - C[0][0]*C[1][1] - C[1][1]*C[2][2] - C[2][2]*C[0][0];
        }        
        
        double  dW_dI1 = Get_dW_dI1(I1,I2);
        double  dW_dI2; // only computed if DIM==3
        
        double  d2W_dI1;
        double  d2W_dI2;
        double  d2W_dI1I2;

        // Compute stress:
        //
        //  T = dW_dE
        //    = 2 * dI1_dC_MN * dI1_dC_MN   +   2 * dI1_dC_MN * dI1_dC_MN  -  p * invC
        //    = 2 * dI1_dC_MN * delta_MN    +   2 * dI1_dC_MN * (I1 delta_MN - C_MN)  -  p * invC         
        T = 2*dW_dI1*identity - pressure*invC;
        if(DIM==3)
        {
            dW_dI2 = Get_dW_dI2(I1,I2);
            T += 2*dW_dI2*(I1*identity - C);
        }
        
        
        // Compute stress derivative if required:
        // 
        // The stress derivative dT_{MN}/dE_{PQ} can be expanded to be seen to be
        //
        //  dT_dE =    4 * true_d2WdI1 * dI1_dC_MN * dI1_dC_PQ
        //           + 4 * true_dWdI1  * d2I1_dC2
        //           + 4 * true_d2WdI2 * dI2_dC_MN * dI2_dC_PQ
        //           + 4 * true_dWdI2  * d2I2_dC2
        //           + 4 * true_d2WdI1I2 * (dI1_dC_MN*dI2_dC_PQ + dI1_dC_PQ*dI2_dC_MN)
        //          - 2 * pressure * d_invC_dC;
        //
        // where
        //   dI1_dC_MN = (M==N); // ie delta_{MN} 
        //   dI1_dC_PQ = (P==Q); 
        //   d2I1_dC2  = 0;
        //
        //   dI2_dC_MN = I1*(M==N)-C[M][N]; 
        //   dI2_dC_PQ = I1*(P==Q)-C[P][Q];                          
        //   d2I2_dC2  = (M==N)*(P==Q)-(M==P)*(N==Q);
        //                
        //   d_invC_dC = -invC[M][P]*invC[Q][N];
        if(computeDTdE)
        {
            d2W_dI1 = Get_d2W_dI1(I1,I2);
            
            if(DIM==3)
            {
                d2W_dI2   = Get_d2W_dI2(I1,I2);
                d2W_dI1I2 = Get_d2W_dI1I2(I1,I2);
            }

            for(unsigned M=0;M<DIM;M++)
            {
                for(unsigned N=0;N<DIM;N++)
                {
                    for(unsigned P=0;P<DIM;P++)
                    {
                        for(unsigned Q=0;Q<DIM;Q++)
                        {
                            dTdE[M][N][P][Q]  =    4 * d2W_dI1  * (M==N) * (P==Q)
                                                 + 2 * pressure * invC[M][P] * invC[Q][N];

                            if(DIM==3)
                            {
                                dTdE[M][N][P][Q] +=    4 * d2W_dI2   * (I1*(M==N)-C[M][N]) * (I1*(P==Q)-C[P][Q])
                                                     + 4 * dW_dI2    * ((M==N)*(P==Q)-(M==P)*(N==Q))
                                                     + 4 * d2W_dI1I2 * ((M==N)*(I1*(P==Q)-C[P][Q]) + (P==Q)*(I1*(M==N)-C[M][N]));
                            }
                        }
                    }
                }
            }
        }
    }

    virtual ~AbstractIsotropicIncompressibleMaterialLaw()
    {
    }

};



#endif /*ABSTRACTISOTROPICINCOMPRESSIBLEMATERIALLAW_HPP_*/
