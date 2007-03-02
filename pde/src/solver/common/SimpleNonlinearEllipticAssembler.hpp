#ifndef _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_


#include <vector>
#include <petscsnes.h>
#include <petscvec.h>
#include <petscmat.h>

#include "AbstractNonlinearStaticAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractBasisFunction.hpp"
#include "AbstractNonlinearEllipticPde.hpp"

                                    
                                    
/**
 * Concrete simple class that assembles and solves the nonlinear system
 * for a nonlinear elliptic PDE. 
 * 
 * USAGE: call the constructor with the mesh, pde and boundary conditions,
 * then call Solve() with the initial guess.
 *
 * ///\todo [old todo, maybe not true anymore after refactor(?)]
 * This class could do with some tidying. More (3D) tests are also needed.
 * It probably needs re-writing to take advantage of parallel machines.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SimpleNonlinearEllipticAssembler : public AbstractNonlinearStaticAssembler<ELEMENT_DIM, SPACE_DIM, 1>
{
    // Allow tests to access private members, in order to test computation of
    // residual & jacobian directly.
    friend class TestSimpleNonlinearEllipticAssembler;
    
private:
    /*< The pde to be solved */
    AbstractNonlinearEllipticPde<SPACE_DIM> *mpNonlinearEllipticPde;      
 
 
    /**
     *  This method returns the matrix to be added to element stiffness matrix
     *  for a given gauss point. The arguments are the bases, bases gradients, 
     *  x and current solution computed at the Gauss point. The returned matrix
     *  will be multiplied by the gauss weight and jacobian determinent and 
     *  added to the element stiffness matrix (see AssembleOnElement()).
     */
    virtual c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        Point<SPACE_DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,SPACE_DIM> &rGradU) 
    {
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> ret;

        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> f_of_u = mpNonlinearEllipticPde->ComputeDiffusionTerm(rX,u(0));
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> f_of_u_prime = mpNonlinearEllipticPde->ComputeDiffusionTermPrime(rX,u(0));
        
        //LinearSourceTerm(x)   not needed as it is a constant wrt u
        double forcing_term_prime = mpNonlinearEllipticPde->ComputeNonlinearSourceTermPrime(rX, u(0));

        // note rGradU is a 1 by SPACE_DIM matrix, the 1 representing the dimension of
        // u (ie in this problem the unknown is a scalar). rGradU0 is rGradU as a vector
        matrix_row< c_matrix<double, 1, SPACE_DIM> > rGradU0( rGradU, 0);
        c_vector<double, ELEMENT_DIM> temp1 = prod(f_of_u_prime,rGradU0);
        c_vector<double, ELEMENT_DIM+1> temp1a = prod(temp1, rGradPhi);
        
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> integrand_values1 = outer_prod(temp1a, rPhi);
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> temp2 = prod(f_of_u, rGradPhi);
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> integrand_values2 = prod(trans(rGradPhi), temp2);
        c_vector<double, ELEMENT_DIM+1> integrand_values3 = forcing_term_prime * rPhi;
        
        ret = integrand_values1 + integrand_values2 - outer_prod( scalar_vector<double>(ELEMENT_DIM+1), integrand_values3);
        
        return ret;
    }
       

    /**
     *  This method returns the vector to be added to element stiffness vector
     *  for a given gauss point. The arguments are the bases, 
     *  x and current solution computed at the Gauss point. The returned vector
     *  will be multiplied by the gauss weight and jacobian determinent and 
     *  added to the element stiffness matrix (see AssembleOnElement()).
     */
    virtual c_vector<double,1*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        Point<SPACE_DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,SPACE_DIM> &rGradU) 
    {
        c_vector<double, 1*(ELEMENT_DIM+1)> ret;
        
       //c_vector<double, SPACE_DIM> gradU = prod(grad_phi, Ui);
        
        // For solving NonlinearEllipticEquation
        // which should be defined in/by NonlinearEllipticEquation.hpp:
        // d/dx [f(U,x) du/dx ] = -g
        // where g(x,U) is the forcing term
        double ForcingTerm = mpNonlinearEllipticPde->ComputeLinearSourceTerm(rX);
        ForcingTerm += mpNonlinearEllipticPde->ComputeNonlinearSourceTerm(rX, u(0));
        //make RHS general: consists of linear and nonlinear source terms
        
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> FOfU = mpNonlinearEllipticPde->ComputeDiffusionTerm(rX,u(0));
        
        // note rGradU is a 1 by SPACE_DIM matrix, the 1 representing the dimension of
        // u (ie in this problem the unknown is a scalar). rGradU0 is rGradU as a vector.
        matrix_row< c_matrix<double, 1, SPACE_DIM> > rGradU0( rGradU, 0);
        c_vector<double, ELEMENT_DIM+1> integrand_values1 =
            prod(c_vector<double, ELEMENT_DIM>(prod(rGradU0, FOfU)), rGradPhi);
           
        ret = integrand_values1 - (ForcingTerm * rPhi);
        return ret;
    }        
        
        
    /**
     *  This method returns the vector to be added to element stiffness vector
     *  for a given gauss point in BoundaryElement. The arguments are the bases, 
     *  x and current solution computed at the Gauss point. The returned vector
     *  will be multiplied by the gauss weight and jacobian determinent and 
     *  added to the element stiffness matrix (see AssembleOnElement()).
     */
    virtual c_vector<double, 1*ELEMENT_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
        c_vector<double, ELEMENT_DIM> &rPhi,
        Point<SPACE_DIM> &rX )
    {
        double Dgradu_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX);
        
        // I'm not sure why we want -phi, but it seems to work:)
        return  (-Dgradu_dot_n)* rPhi ;
    }
           

public :

    /** 
     * Constructor - takes in the mesh, pde and boundary conditions container to be solved. Can
     * also define the number of quad points (in each dimension), the default value of which is 2
     */
    SimpleNonlinearEllipticAssembler( ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                      AbstractNonlinearEllipticPde<SPACE_DIM>* pPde,
                                      BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>* pBoundaryConditions,
                                      unsigned numQuadPoints = 2) :  
            AbstractNonlinearStaticAssembler<ELEMENT_DIM,SPACE_DIM,1>(numQuadPoints)
    {
        // Store data structures
        assert(pMesh!=NULL);
        assert(pPde!=NULL);
        assert(pBoundaryConditions!=NULL);

        this->mpMesh = pMesh;
        mpNonlinearEllipticPde = pPde;
        this->mpBoundaryConditions = pBoundaryConditions;
    }

    /**
     * Alternative constructor - takes in the mesh, pde and boundary conditions container to be solved, and 
     * basis functions to use. Can also define the number of quad points (in each dimension), the default 
     * value of which is 2
     */
    SimpleNonlinearEllipticAssembler( ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                      AbstractNonlinearEllipticPde<SPACE_DIM>* pPde,
                                      BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>* pBoundaryConditions,
                                      AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                      AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                      unsigned numQuadPoints = 2) :
            AbstractNonlinearStaticAssembler<ELEMENT_DIM,SPACE_DIM,1>(pBasisFunction, pSurfaceBasisFunction, numQuadPoints)
    {
        // Store data structures
        assert(pMesh!=NULL);
        assert(pPde!=NULL);
        assert(pBoundaryConditions!=NULL);

        this->mpMesh = pMesh;
        mpNonlinearEllipticPde = pPde;
        this->mpBoundaryConditions = pBoundaryConditions;
    }  
};




#endif  // _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
