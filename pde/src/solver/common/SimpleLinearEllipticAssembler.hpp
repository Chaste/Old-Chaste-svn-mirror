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
#include "AbstractBasisFunction.hpp"


/**
 *  SimpleLinearEllipticAssembler
 *
 *  Assembler for solving AbstractLinearEllipticPdes
 */
template<int ELEMENT_DIM, int SPACE_DIM>
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
    virtual c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)> ComputeLhsTerm(
        const c_vector<double,ELEMENT_DIM+1> &rPhi,
        const c_matrix<double,ELEMENT_DIM,ELEMENT_DIM+1> &rGradPhi,
        const Point<SPACE_DIM> &rX,
        const c_vector<double,1> &u)
    {
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> pde_diffusion_term = mpEllipticPde->ComputeDiffusionTerm(rX);
        
        return prod( trans(rGradPhi), c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) );
    }
    
    /**
     *  The term arising from boundary conditions to be added to the element
     *  stiffness vector
     */
    virtual c_vector<double,1*(ELEMENT_DIM+1)> ComputeRhsTerm(const c_vector<double, ELEMENT_DIM+1> &rPhi,
                                                              const Point<SPACE_DIM> &rX,
                                                              const c_vector<double,1> &u)
    {
        return mpEllipticPde->ComputeLinearSourceTerm(rX) * rPhi;
    }
    
    
    virtual c_vector<double, ELEMENT_DIM> ComputeSurfaceRhsTerm(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
                                                                const c_vector<double, ELEMENT_DIM> &phi,
                                                                const Point<SPACE_DIM> &x )
    {
        // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
        double D_times_gradu_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, x);
        return phi * D_times_gradu_dot_n;
    }
    
    
    
public:
    /**
     * Constructor stores the mesh, pde and boundary conditons, and calls base constructor.
     */
    SimpleLinearEllipticAssembler(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                  AbstractLinearEllipticPde<SPACE_DIM>* pPde,
                                  BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
                                  int numQuadPoints = 2) :
            AbstractLinearStaticProblemAssembler<ELEMENT_DIM,SPACE_DIM,1>(numQuadPoints)
    {
        // note - we don't check any of these are NULL here (that is done in Solve() instead),
        // to allow the user or a subclass to set any of these later
        mpEllipticPde = pPde;
        this->mpMesh = pMesh;
        this->mpBoundaryConditions = pBoundaryConditions;
    }
    
    /**
     * Constructor which also takes in basis functions
     */
    SimpleLinearEllipticAssembler(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                  AbstractLinearEllipticPde<SPACE_DIM>* pPde,
                                  BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
                                  AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                  AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                  int numQuadPoints = 2) :
            AbstractLinearStaticProblemAssembler<ELEMENT_DIM,SPACE_DIM,1>(pBasisFunction, pSurfaceBasisFunction, numQuadPoints)
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
        assert(mpEllipticPde != NULL);
        assert(this->mpMesh != NULL);
        assert(this->mpBoundaryConditions != NULL);
    }
};


#endif //_SIMPLELINEARELLIPTICASSEMBLER_HPP_
