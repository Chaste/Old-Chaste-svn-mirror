#ifndef ABSTRACTINCOMPRESSIBLEMATERIALLAW_HPP_
#define ABSTRACTINCOMPRESSIBLEMATERIALLAW_HPP_


#include <numerics/vectors.h>
#include <numerics/matrices.h>

#include <base/tensor.h>
#include <base/symmetric_tensor.h>


/** 
 *  AbstractIncompressibleMaterialLaw
 * 
 *  An incompressible hyperelastic material law for finite elastiticy
 *  
 *  The law is given by a strain energy function W(E), where E is the strain, such
 *  that the stress T = dW/dE
 */
template<int DIM>
class AbstractIncompressibleMaterialLaw
{
public :
    /**
     *  Compute the stress T and the stress derivative dT/dE for a given strain.
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
     */  
    virtual void GetStressAndStressDerivative(Tensor<2,DIM>           C, 
                                              Tensor<2,DIM>           invC,
                                              double                  pressure,
                                              SymmetricTensor<2,DIM>& T,
                                              double                  dTdE[DIM][DIM][DIM][DIM],
                                              bool                    computeDTdE)=0;
                                              
    virtual ~AbstractIncompressibleMaterialLaw()
    {
    }
};


#endif /*ABSTRACTINCOMPRESSIBLEMATERIALLAW_HPP_*/
