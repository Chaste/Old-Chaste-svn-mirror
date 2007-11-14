#ifndef _SIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLELINEARELLIPTICASSEMBLER_HPP_


#include <vector>
#include <petscvec.h>

#include "LinearSystem.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AbstractLinearAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"


/**
 *  SimpleLinearEllipticAssembler
 *
 *  Assembler for solving AbstractLinearEllipticPdes
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SimpleLinearEllipticAssembler : public AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 1>
{
    friend class TestSimpleLinearEllipticAssembler;
    
private:
    AbstractLinearEllipticPde<SPACE_DIM>* mpEllipticPde;
    
protected:
    /**
     *  The term to be added to the element stiffness matrix: 
     *  
     *   grad_phi[row] \dot ( pde_diffusion_term * grad_phi[col]) 
     */
    virtual c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,SPACE_DIM> &rGradU)
    {
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> pde_diffusion_term = mpEllipticPde->ComputeDiffusionTerm(rX);
        
        // if statement just saves computing phi*phi^T if it is to be multiplied by zero
        if(mpEllipticPde->ComputeLinearInUCoeffInSourceTerm(rX)!=0)
        {
            return   prod( trans(rGradPhi), c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) )
                   - mpEllipticPde->ComputeLinearInUCoeffInSourceTerm(rX)*outer_prod(rPhi,rPhi);
        }
        else
        {
            return   prod( trans(rGradPhi), c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) );
        }
    }
    
    /**
     *  The term arising from boundary conditions to be added to the element
     *  stiffness vector
     */
    virtual c_vector<double,1*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,SPACE_DIM> &rGradU)
    {
        return mpEllipticPde->ComputeConstantInUSourceTerm(rX) * rPhi;
    }
    
    
    
    virtual c_vector<double, ELEMENT_DIM> ComputeVectorSurfaceTerm(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
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
    SimpleLinearEllipticAssembler(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                  AbstractLinearEllipticPde<SPACE_DIM>* pPde,
                                  BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
                                  unsigned numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>(),
            AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,1>(numQuadPoints)
    {
        // note - we don't check any of these are NULL here (that is done in Solve() instead),
        // to allow the user or a subclass to set any of these later
        mpEllipticPde = pPde;
        this->SetMesh(pMesh);
        this->SetBoundaryConditionsContainer(pBoundaryConditions);
    }
    
    /**
     * This method is called at the beginning of Solve() in AbstractLinearStaticProblemAssembler
     */
    void PrepareForSolve()
    {
        AbstractLinearAssembler<ELEMENT_DIM,SPACE_DIM,1>::PrepareForSolve();
        assert(mpEllipticPde != NULL);
    }
};


#endif //_SIMPLELINEARELLIPTICASSEMBLER_HPP_
