#ifndef _SIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLELINEARELLIPTICASSEMBLER_HPP_


#include <vector>
#include <petscvec.h>

#include "LinearSystem.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AbstractLinearStaticProblemAssembler.hpp"
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
class SimpleLinearEllipticAssembler : public AbstractLinearStaticProblemAssembler<ELEMENT_DIM, SPACE_DIM, 1>
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
        Point<SPACE_DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,SPACE_DIM> &rGradU)
    {
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> pde_diffusion_term = mpEllipticPde->ComputeDiffusionTerm(rX);
        
        return prod( trans(rGradPhi), c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) );
    }
    
    /**
     *  The term arising from boundary conditions to be added to the element
     *  stiffness vector
     */
    virtual c_vector<double,1*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        Point<SPACE_DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,SPACE_DIM> &rGradU)
    {
        return mpEllipticPde->ComputeLinearSourceTerm(rX) * rPhi;
    }
    
    
    
    virtual c_vector<double, ELEMENT_DIM> ComputeVectorSurfaceTerm(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
            c_vector<double, ELEMENT_DIM> &rPhi,
            Point<SPACE_DIM> &rX )
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
            AbstractLinearStaticProblemAssembler<ELEMENT_DIM,SPACE_DIM,1>(numQuadPoints)
    {
        // note - we don't check any of these are NULL here (that is done in Solve() instead),
        // to allow the user or a subclass to set any of these later
        mpEllipticPde = pPde;
        this->mpMesh = pMesh;
        this->mpBoundaryConditions = pBoundaryConditions;
    }
    
    /**
     * This method is called at the beginning of Solve() in AbstractLinearStaticProblemAssembler
     */
    void PrepareForSolve()
    {
        AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>::PrepareForSolve();
        assert(mpEllipticPde != NULL);
        assert(this->mpMesh != NULL);
        assert(this->mpBoundaryConditions != NULL);
    }
};


#endif //_SIMPLELINEARELLIPTICASSEMBLER_HPP_
