#ifndef _SIMPLEDG0PARABOLICASSEMBLER_HPP_
#define _SIMPLEDG0PARABOLICASSEMBLER_HPP_

#include <vector>
#include <iostream>
#include <petscvec.h>

#include "LinearSystem.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "AbstractLinearAssembler.hpp"
#include "AbstractDynamicAssemblerMixin.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"

/**
 *  SimpleDg0ParabolicAssembler
 *
 *  Assembler for solving AbstractLinearParabolicPdes
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SimpleDg0ParabolicAssembler
    : public AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 1>,
      public AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, 1>
{
private:
    AbstractLinearParabolicPde<SPACE_DIM>* mpParabolicPde;
    
protected:
    /**
     *  The term to be added to the element stiffness matrix: 
     *  
     *   grad_phi[row] \dot ( pde_diffusion_term * grad_phi[col]) + 
     *  (1.0/mDt) * pPde->ComputeDuDtCoefficientFunction(rX) * rPhi[row] * rPhi[col]
     */
    virtual c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,SPACE_DIM> &rGradU /* not used */ )
    {
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> pde_diffusion_term = mpParabolicPde->ComputeDiffusionTerm(rX);
        
        return    prod( trans(rGradPhi), c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) )
                  + this->mDtInverse * mpParabolicPde->ComputeDuDtCoefficientFunction(rX) * outer_prod(rPhi, rPhi);
    }
    
    /**
     *  The term to be added to the element stiffness vector: 
     */
    virtual c_vector<double,1*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double, 1, SPACE_DIM> &rGradU /* not used */ )
            
    {
        return (mpParabolicPde->ComputeNonlinearSourceTerm(rX, u(0)) + mpParabolicPde->ComputeLinearSourceTerm(rX)
                + this->mDtInverse * mpParabolicPde->ComputeDuDtCoefficientFunction(rX) * u(0)) * rPhi;
    }
    
    
    /**
     *  The term arising from boundary conditions to be added to the element
     *  stiffness vector
     */
    virtual c_vector<double, ELEMENT_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
        c_vector<double, ELEMENT_DIM> &rPhi,
        ChastePoint<SPACE_DIM> &rX )
    {
        // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
        double D_times_gradu_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX);
        return rPhi * D_times_gradu_dot_n;
    }
    
    
public:
    /**
     * Constructor stores the mesh, pde and boundary conditons, and calls base constructor.
     */
    SimpleDg0ParabolicAssembler(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                AbstractLinearParabolicPde<SPACE_DIM>* pPde,
                                BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
                                unsigned numQuadPoints = 2,
                                double linearSolverRelativeTolerance=1e-6) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>(),
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,1>(numQuadPoints,linearSolverRelativeTolerance),
            AbstractDynamicAssemblerMixin<ELEMENT_DIM,SPACE_DIM,1>()
    {
        // note - we don't check any of these are NULL here (that is done in Solve() instead),
        // to allow the user or a subclass to set any of these later
        mpParabolicPde = pPde;
        this->SetMesh(pMesh);
        this->SetBoundaryConditionsContainer(pBoundaryConditions);
        
        this->SetMatrixIsConstant();
    }
    
    /**
     * Called by AbstractDynamicAssemblerMixin at the beginning of Solve() 
     */
    virtual void PrepareForSolve()
    {
        AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,1>::PrepareForSolve();
        assert(mpParabolicPde != NULL);
    }
    
    Vec Solve(Vec currentSolutionOrGuess=NULL, double currentTime=0.0)
    {
        return AbstractDynamicAssemblerMixin<ELEMENT_DIM,SPACE_DIM,1>::Solve(currentSolutionOrGuess,currentTime);
    }
};


#endif //_SIMPLEDG0PARABOLICASSEMBLER_HPP_
