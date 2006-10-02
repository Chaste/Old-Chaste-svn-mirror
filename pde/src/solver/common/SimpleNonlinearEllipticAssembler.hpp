#ifndef _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_


#include <vector>
#include <petscsnes.h>
#include <petscvec.h>
#include <petscmat.h>

#include "AbstractAssembler.hpp"
#include "AbstractNonlinearEllipticAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractBasisFunction.hpp"
#include "AbstractNonlinearEllipticPde.hpp"
#include "ReplicatableVector.hpp"


                                    
                                    
/**
 * Concrete simple class that assembles and solves the nonlinear system
 * for a nonlinear elliptic PDE.
 *
 * \todo This class could do with some tidying. More (3D) tests are also needed.
 * It probably needs re-writing to take advantage of parallel machines.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
class SimpleNonlinearEllipticAssembler : public AbstractNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>
{
    // Allow tests to access private members, in order to test computation of
    // residual & jacobian directly.
    friend class TestSimpleNonlinearEllipticAssembler;
    
private:
    AbstractNonlinearEllipticPde<SPACE_DIM> *mpPde;
    
    
    void AssembleResidualOnElement( const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                                   c_vector<double, ELEMENT_DIM+1> &rBElem,
                                   c_vector<double, ELEMENT_DIM+1> Ui
                                  );


    void AssembleResidualOnSurfaceElement( const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
                                          c_vector<double, ELEMENT_DIM> &rBsubElem,                                          
                                          vector<double> Ui
                                         );


    void AssembleJacobianOnElement( const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                                   c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> &rAElem,
                                   vector<double> Ui
                                  );
        
public :

    /**
     * Constructors just call the base class versions.
     */
    SimpleNonlinearEllipticAssembler( ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                      AbstractNonlinearEllipticPde<SPACE_DIM>* pPde,
                                      BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>* pBoundaryConditions,
                                      int numQuadPoints = 2) :  
            AbstractNonlinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM>(numQuadPoints)
    {
        // Store data structures
        assert(pMesh!=NULL);
        assert(pPde!=NULL);
        assert(pBoundaryConditions!=NULL);

        this->mpMesh = pMesh;
        mpPde = pPde;
        this->mpBoundaryConditions = pBoundaryConditions;
    }

    SimpleNonlinearEllipticAssembler( ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                      AbstractNonlinearEllipticPde<SPACE_DIM>* pPde,
                                      BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>* pBoundaryConditions,
                                      AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                      AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                      int numQuadPoints = 2) :
            AbstractNonlinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, numQuadPoints)
    {
        // Store data structures
        assert(pMesh!=NULL);
        assert(pPde!=NULL);
        assert(pBoundaryConditions!=NULL);

        this->mpMesh = pMesh;
        mpPde = pPde;
        this->mpBoundaryConditions = pBoundaryConditions;
    }
            
    

};


/*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * *****************************************************************************
 *                     Computation of Residual Vector
 * *****************************************************************************
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */


/**
 * Compute Residual on Surface Elements, applying Neumann boundary conditions.
 *
 */
template<int ELEMENT_DIM, int SPACE_DIM>
void SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::AssembleResidualOnSurfaceElement(
    const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
    c_vector<double, ELEMENT_DIM> &rBsubElem,
    vector<double> Ui)
{
    AbstractBasisFunction<ELEMENT_DIM-1> &rBasisFunction =
        *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>::mpSurfaceBasisFunction);
    GaussianQuadratureRule<ELEMENT_DIM-1> &quad_rule =
        *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>::mpSurfaceQuadRule);
        
    double jacobian_determinant = rSurfaceElement.GetJacobianDeterminant();
    
    //const int num_nodes = rSurfaceElement.GetNumNodes();
    
    for (int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
    {
        Point<ELEMENT_DIM-1> quad_point=quad_rule.GetQuadPoint(quad_index);
        
        c_vector<double, ELEMENT_DIM>  phi = rBasisFunction.ComputeBasisFunctions(quad_point);
        
        // location of the gauss point in the original element will be stored in x
        Point<SPACE_DIM> x(0,0,0);
        
        for (int i=0; i<rSurfaceElement.GetNumNodes(); i++)
        {
        
            for (int j=0; j<SPACE_DIM; j++)
            {
                x.SetCoordinate(j, x[j] + phi[i]*rSurfaceElement.GetNodeLocation(i,j));
            }
        }
        
        /**
         * \todo Neumann BC value depends on u?
         */
        //double U = inner_prod(phi,Ui);
        double Dgradu_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, x);
        
        // I'm not sure why we want -phi, but it seems to work:)
        
        noalias(rBsubElem) += (Dgradu_dot_n * jacobian_determinant * quad_rule.GetWeight(quad_index) * -1) * phi ;
    }
}


/**
 * Compute the residual vector on a single Element
 *
 */
template<int ELEMENT_DIM, int SPACE_DIM>
void SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::AssembleResidualOnElement(
    const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
    c_vector<double, ELEMENT_DIM+1> &rBElem,
    c_vector<double, ELEMENT_DIM+1> Ui)
{
    AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
        *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>::mpBasisFunction);
    GaussianQuadratureRule<ELEMENT_DIM> *pQuadRule =
        AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>::mpQuadRule;
        
    const c_matrix<double, SPACE_DIM, SPACE_DIM> *inverseJacobian = rElement.GetInverseJacobian();
    double jacobian_determinant = rElement.GetJacobianDeterminant();
    
    const int num_nodes = rElement.GetNumNodes();
    
    for (int quad_index=0; quad_index<pQuadRule->GetNumQuadPoints(); quad_index++)
    {
        Point<ELEMENT_DIM> quad_point=pQuadRule->GetQuadPoint(quad_index);
        
        c_vector<double, ELEMENT_DIM+1> phi     = rBasisFunction.ComputeBasisFunctions(quad_point);
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1 > grad_phi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
                                                                 (quad_point, *inverseJacobian);
                                                                 
        Point<SPACE_DIM> x(0,0,0);
        for (int i=0; i<num_nodes; i++)
        {
            for (int j=0; j<SPACE_DIM; j++)
            {
                x.SetCoordinate(j, x[j] + phi[i]*rElement.GetNodeLocation(i,j));
                
            }
        }
        
        // Need to compute add U as double and gradU as vector double
        // U = sum(Ui phi_i)
        double U = inner_prod(phi, Ui);
        c_vector<double, SPACE_DIM> gradU = prod(grad_phi, Ui);
        
        // For solving NonlinearEllipticEquation
        // which should be defined in/by NonlinearEllipticEquation.hpp:
        // d/dx [f(U,x) du/dx ] = -g
        // where g(x,U) is the forcing term
        double ForcingTerm = mpPde->ComputeLinearSourceTerm(x);
        ForcingTerm += mpPde->ComputeNonlinearSourceTerm(x, U);
        //make RHS general: consists of linear and nonlinear source terms
        
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> FOfU = mpPde->ComputeDiffusionTerm(x,U);
        
        c_vector<double, ELEMENT_DIM+1> integrand_values1 =
            prod(c_vector<double, ELEMENT_DIM>(prod(gradU, FOfU)), grad_phi);
            
        noalias(rBElem) += (jacobian_determinant*pQuadRule->GetWeight(quad_index))
                           * (integrand_values1-(ForcingTerm * phi)) ;
    }
}






/*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * *****************************************************************************
 *                     Computation of Jacobian Matrix
 * *****************************************************************************
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */







template<int ELEMENT_DIM, int SPACE_DIM>
void SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::AssembleJacobianOnElement(
    const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> &rAElem,
    vector<double> Ui)
{
    AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
        *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>::mpBasisFunction);
    GaussianQuadratureRule<ELEMENT_DIM> *pQuadRule =
        AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>::mpQuadRule;
        
    const c_matrix<double, SPACE_DIM, SPACE_DIM> *inverseJacobian = rElement.GetInverseJacobian();
    double jacobian_determinant = rElement.GetJacobianDeterminant();
    
    // Initialise element contributions to zero
    const int num_nodes = rElement.GetNumNodes();
    
    for (int quad_index=0; quad_index<pQuadRule->GetNumQuadPoints(); quad_index++)
    {
        Point<ELEMENT_DIM> quad_point=pQuadRule->GetQuadPoint(quad_index);
        
        c_vector<double, ELEMENT_DIM+1> phi     = rBasisFunction.ComputeBasisFunctions(quad_point);
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1 > grad_phi = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
                                                                 (quad_point, *inverseJacobian);
                                                                 
        Point<SPACE_DIM> x(0,0,0);
        //Need to compute add U as double and gradU as vector double
        // get U =sum(Ui phi_i)
        double U = inner_prod(phi, Ui);
        c_vector<double, SPACE_DIM> gradU=prod(grad_phi, Ui);
        
        for (int i=0; i<num_nodes; i++)
        {
            for (int j=0; j<SPACE_DIM; j++)
            {
                x.SetCoordinate(j, x[j] + phi[i]*rElement.GetNodeLocation(i,j));
                
            }
        }
        
        // For solving NonlinearEllipticEquation
        // which should be defined in/by NonlinearEllipticEquation.hpp:
        // d/dx [f(U,x) du/dx ] = -g
        // where g(x,U) is the forcing term
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> f_of_u = mpPde->ComputeDiffusionTerm(x,U);
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> f_of_u_prime = mpPde->ComputeDiffusionTermPrime(x,U);
        
        //LinearSourceTerm(x)   not needed as it is a constant wrt U_i
        double forcing_term_prime = mpPde->ComputeNonlinearSourceTermPrime(x, U);
        c_vector<double, ELEMENT_DIM> temp1 = prod(f_of_u_prime,gradU);
        c_vector<double, ELEMENT_DIM+1> temp1a = prod(temp1, grad_phi);
        
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> integrand_values1 = outer_prod(temp1a, phi);
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> temp2 = prod(f_of_u, grad_phi);
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> integrand_values2 = prod(trans(grad_phi), temp2);
        c_vector<double, ELEMENT_DIM+1> integrand_values3 = forcing_term_prime * phi;
        
        rAElem +=  jacobian_determinant*pQuadRule->GetWeight(quad_index)*
                   (integrand_values1 + integrand_values2 -
                    outer_prod( scalar_vector<double>(ELEMENT_DIM+1), integrand_values3));
                    
    }
}






#endif  // _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
