/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ABSTRACTINCOMPRESSIBLEMATERIALLAW_HPP_
#define ABSTRACTINCOMPRESSIBLEMATERIALLAW_HPP_


#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <base/tensor.h>
#include <base/symmetric_tensor.h>
#include "FourthOrderTensor.hpp"

/**
 *  AbstractIncompressibleMaterialLaw
 *
 *  An incompressible hyperelastic material law for finite elastiticy
 *
 *  The law is given by a strain energy function W(E), where E is the strain, such
 *  that the stress T = dW/dE
 */
template<unsigned DIM>
class AbstractIncompressibleMaterialLaw
{
public :
    /**
     *  Compute the (2nd Piola Kirchoff) stress T and the stress derivative dT/dE for 
     *  a given strain.
     *  
     *  NOTE: the strain E is not expected to be passed in, instead the Lagrangian
     *  deformation tensor C is required (recall, E = 0.5(C-I))
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
    virtual void ComputeStressAndStressDerivative(Tensor<2,DIM>&          C,
                                                  Tensor<2,DIM>&          invC,
                                                  double                  pressure,
                                                  SymmetricTensor<2,DIM>& T,
                                                  FourthOrderTensor<DIM>& dTdE,
                                                  bool                    computeDTdE)=0;
                                                  
                                                  
    /**
     *  Compute the Cauchy stress (the true stress), given the deformation gradient
     *  F and the pressure. The Cauchy stress is given by
     *  
     *  sigma^{ij} = (1/detF) F^i_M T^{MN} F^j_N
     * 
     *  where T is the 2nd Piola Kirchoff stress, dW/dE
     * 
     *  @param F the deformation gradient
     *  @param pressure the pressure
     *  @sigma sigma an empty matrix, which will be filled in with the Cauchy stress
     * 
     *  Note: the compute the material part of the stress (the pressure-independent
     *  part), just pass in pressure=0.0
     */
    void ComputeCauchyStress(Tensor<2,DIM>& F, double pressure, Tensor<2,DIM>& sigma)
    {
        double detF = determinant(F);
        
        Tensor<2,DIM> C = transpose(F) * F;
        Tensor<2,DIM> invC = invert(C);
        
        SymmetricTensor<2,DIM> T;
        
        static FourthOrderTensor<DIM> dTdE; // not filled in, made static for efficiency
        
        ComputeStressAndStressDerivative(C,invC,pressure,T,dTdE,false);
        
        // looping it probably more eficient then doing sigma = (1/detF)F*T*transpose(F)
        // which doesn't seem to compile anyway, as F is a Tensor<2,DIM> and T is a
        // SymmetricTensor<2,DIM>
        for (unsigned i=0; i<DIM; i++)
        {
            for (unsigned j=0; j<DIM; j++)
            {
                sigma[i][j] = 0.0;
                for (unsigned M=0; M<DIM; M++)
                {
                    for (unsigned N=0; N<DIM; N++)
                    {
                        sigma[i][j] += F[i][M]*T[M][N]*F[j][N];
                    }
                }
                sigma[i][j] /= detF;
            }
        }
    }
    
    
    
    /**
     *  Compute the 1st Piola Kirchoff stress, given the deformation gradient F
     *  and the pressure. The 1st Piola Kirchoff stress given by
     *  
     *  S^{Mi} = T^{MN} F^i_M, 
     * 
     *  where T is the 2nd PK stress, dW/dE. 
     * 
     *  Note that this stress is not symmetric and the least useful of the three
     *  stresses. 
     * 
     *  @param F the deformation gradient
     *  @param pressure the pressure
     *  @sigma S an empty matrix, which will be filled in with the stress
     * 
     *  Note: the compute the material part of the stress (the pressure-independent
     *  part), just pass in pressure=0.0
     */
    void Compute1stPiolaKirchoffStress(Tensor<2,DIM>& F, double pressure, Tensor<2,DIM>& S)
    {
        Tensor<2,DIM> C = transpose(F) * F;
        Tensor<2,DIM> invC = invert(C);
        
        SymmetricTensor<2,DIM> T;
        
        static FourthOrderTensor<DIM> dTdE; // not filled in, made static for efficiency
        
        ComputeStressAndStressDerivative(C,invC,pressure,T,dTdE,false);
        
        
        // looping it probably more eficient then doing S = T*transpose(F)
        // which doesn't seem to compile anyway, as F is a Tensor<2,DIM> and T is a
        // SymmetricTensor<2,DIM>
        for (unsigned M=0; M<DIM; M++)
        {
            for (unsigned i=0; i<DIM; i++)
            {
                S[M][i] = 0.0;
                for (unsigned N=0; N<DIM; N++)
                {
                    S[M][i] += T[M][N] * F[i][N];
                }
            }
        }
    }
    
    
    
    /**
     *  Compute the 2nd Piola Kirchoff stress, given the deformation tensor C
     *  and the pressure. The 2nd Piola Kirchoff stress given by
     *  
     *  T^{MN} = dW/dE_{MN} = 2dW/dC_{MN} 
     * 
     *  @param C the Lagrange deformation tensor (C=F^T F), *not* F, and *not* E
     *  @param pressure the pressure
     *  @sigma T an empty matrix, which will be filled in with the stress
     * 
     *  Note: to compute the material part of the stress (the pressure-independent
     *  part), just pass in pressure=0.0
     */
    void Compute2ndPiolaKirchoffStress(Tensor<2,DIM>& C, double pressure, SymmetricTensor<2,DIM>& T)
    {
        Tensor<2,DIM> invC = invert(C);
        
        static FourthOrderTensor<DIM> dTdE; // not filled in, made static for efficiency
        
        ComputeStressAndStressDerivative(C,invC,pressure,T,dTdE,false);
    }
    
    /**
     *  Get the pressure corresponding to E=0, ie C=identity
     */
    virtual double GetZeroStrainPressure()=0;
    
    virtual ~AbstractIncompressibleMaterialLaw()
    {}
    
    /** 
     *  Set a scale factor by which (dimensional) material parameters are scaled. This method
     *  can be optionally implemented in the child class; if no implementation is made an
     *  exception is thrown. A scale factor may be used/needed to improve GMRES convergence.
     *  Note that is a material law is scaled like this any dimensionally equivalent terms
     *  (eg gravity, tractions, active tensions) must also be scaled. Also, computed pressure
     *  will come out scaled.
     */
    virtual void ScaleMaterialParameters(double scaleFactor)
    {
        EXCEPTION("[the material law you are using]::ScaleMaterialParameters() has not be implemented\n");
    }
};


#endif /*ABSTRACTINCOMPRESSIBLEMATERIALLAW_HPP_*/
