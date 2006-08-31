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

/*
 * Since we need to pass function pointers to the PETSc SNES routines, we can't
 * make these functions below methods. This is a pain, since it also means we
 * need to pass round a pointer to our assembler object as the void *pContext,
 * and cast it within the function to access data members.
 *
 * All the functions are defined as stubs which call methods on *pContext.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeResidualPetsc(SNES snes, Vec currentSolution, Vec residualVector,
                                    void *pContext);
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeJacobianPetsc(SNES snes, Vec currentSolution,
                                    Mat *pGlobalJacobian, Mat *pPreconditioner,
                                    MatStructure *pMatStructure, void *pContext);
                                    
                                    
                                    
/**
 * Concrete simple class that assembles and solves the nonlinear system
 * for a nonlinear elliptic PDE.
 *
 * \todo This class could do with some tidying. More (3D) tests are also needed.
 * It probably needs re-writing to take advantage of parallel machines.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
class SimpleNonlinearEllipticAssembler: public AbstractNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>
{
    // Allow tests to access private members, in order to test computation of
    // residual & jacobian directly.
    friend class TestSimpleNonlinearEllipticAssembler;
    
private:
    ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *mpMesh;
    AbstractNonlinearEllipticPde<SPACE_DIM> *mpPde;
    BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1> *mpBoundaryConditions;
    AbstractNonlinearSolver *mpSolver;
    bool mUseAnalyticalJacobian;
    
public:
    // These methods are called indirectly by the PETSc iterative solvers,
    // so must be public.
    PetscErrorCode ComputeResidual(const Vec currentSolution, Vec residualVector);
    PetscErrorCode ComputeJacobianAnalytically(const Vec currentSolution, Mat *pGlobalJacobian);
    PetscErrorCode ComputeJacobianNumerically(const Vec input, Mat *pJacobian);
    PetscErrorCode ComputeJacobian(const Vec input, Mat *pJacobian);
    
private:
    void ComputeResidualOnElement(
        const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
        c_vector<double, ELEMENT_DIM+1> &rBElem,
        AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
        c_vector<double, ELEMENT_DIM+1> Ui);
    void ComputeResidualOnSurfaceElement(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
        c_vector<double, ELEMENT_DIM> &rBsubElem,
        AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM, 1> &rBoundaryConditions,
        vector<double> Ui);
    void ComputeJacobianOnElement(
        const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> &rAElem,
        AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
        vector<double> Ui);
        
public:
    /**
     * Constructors just call the base class versions.
     */
    SimpleNonlinearEllipticAssembler(int numPoints = 2) :
            AbstractNonlinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM>(numPoints)
    {}
    SimpleNonlinearEllipticAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                     AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                     int numPoints = 2) :
            AbstractNonlinearEllipticAssembler<ELEMENT_DIM,SPACE_DIM>(pBasisFunction, pSurfaceBasisFunction, numPoints)
    {}
    
    Vec AssembleSystem(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh,
                       AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
                       BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1> *pBoundaryConditions,
                       AbstractNonlinearSolver *pSolver,
                       Vec initialGuess,
                       bool UseAnalyticalJacobian = false);
                       
};

/**
 * Assemble and solve the system for a nonlinear elliptic PDE.
 *
 * @param pMesh Pointer to the mesh to solve on
 * @param pPde Pointer to the object specifying the equation to solve
 * @param pBoundaryConditions Pointer to the container object for our boundary conditions
 * @param pSolver Pointer to the nonlinear solver object
 * @param pBasisFunction Pointer to object for computing basis functions
 * @param pGaussianQuadratureRule Pointer to database object for Gaussian quadrature
 * @param initialGuess An initial guess for the iterative solver
 * @param UseAnalyticalJacobian Set to true to use an analytically calculated
 *     jacobian matrix rather than a numerically approximated one.
 * @return A PETSc vector giving the solution at each mesh node.
 */
template <int ELEMENT_DIM, int SPACE_DIM>
Vec SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::AssembleSystem(
    ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh,
    AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
    BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1> *pBoundaryConditions,
    AbstractNonlinearSolver *pSolver,
    Vec initialGuess,
    bool UseAnalyticalJacobian)
{
    // Store data structures as private members
    mpMesh = pMesh;
    mpPde = pPde;
    mpBoundaryConditions = pBoundaryConditions;
    
    mUseAnalyticalJacobian = UseAnalyticalJacobian;
    
    Vec residual;
    VecDuplicate(initialGuess, &residual);
    
    Vec answer = pSolver->Solve(&ComputeResidualPetsc<ELEMENT_DIM, SPACE_DIM>,
                                &ComputeJacobianPetsc<ELEMENT_DIM, SPACE_DIM>, residual, initialGuess, this);
                                
    VecDestroy(residual);
    
    return answer;
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
void SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeResidualOnSurfaceElement(
    const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
    c_vector<double, ELEMENT_DIM> &rBsubElem,
    AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM, 1> &rBoundaryConditions,
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
        double Dgradu_dot_n = rBoundaryConditions.GetNeumannBCValue(&rSurfaceElement, x);
        
        // I'm not sure why we want -phi, but it seems to work:)
        
        noalias(rBsubElem) += (Dgradu_dot_n * jacobian_determinant * quad_rule.GetWeight(quad_index) * -1) * phi ;
    }
}


/**
 * Compute the residual vector on a single Element
 *
 */
template<int ELEMENT_DIM, int SPACE_DIM>
void SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeResidualOnElement(
    const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
    c_vector<double, ELEMENT_DIM+1> &rBElem,
    AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
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
        double ForcingTerm = pPde->ComputeLinearSourceTerm(x);
        ForcingTerm += pPde->ComputeNonlinearSourceTerm(x, U);
        //make RHS general: consists of linear and nonlinear source terms
        
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> FOfU = pPde->ComputeDiffusionTerm(x,U);
        
        c_vector<double, ELEMENT_DIM+1> integrand_values1 =
            prod(c_vector<double, ELEMENT_DIM>(prod(gradU, FOfU)), grad_phi);
            
        noalias(rBElem) += (jacobian_determinant*pQuadRule->GetWeight(quad_index))
                           * (integrand_values1-(ForcingTerm * phi)) ;
    }
}

/**
 * Function called by PETSc to compute the residual vector given the current solution guess.
 * Calls a method on a SimpleNonlinearEllipticAssembler object to do the calculation.
 *
 * @param snes This is not used by us, but required by PETSc.
 * @param currentSolution The solution guess for the current iteration.
 * @param residualVector We fill this with the residual vector.
 * @param pContext Pointer to a SimpleNonlinearEllipticAssembler object.
 * @return An error code if any PETSc routines fail.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeResidualPetsc(SNES snes, Vec currentSolution,
                                    Vec residualVector, void *pContext)
{
    SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM> *pAssembler =
        (SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>*) pContext;
        
    return pAssembler->ComputeResidual(currentSolution, residualVector);
}


/**
 * Compute the residual vector given the current solution guess.
 *
 * @param currentSolution The solution guess for the current iteration.
 * @param residualVector We fill this with the residual vector.
 * @return An error code if any PETSc routines fail.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeResidual(
    const Vec currentSolution,
    Vec residualVector)
{
    // Set residual vector to zero
    PetscScalar zero = 0.0;
    
#if (PETSC_VERSION_MINOR == 2) //Old API
    PETSCEXCEPT( VecSet(&zero, residualVector) );
#else
    PETSCEXCEPT( VecSet(residualVector, zero) );
#endif
    
    // Replicate the currentSolution data
    ReplicatableVector current_solution_replicated_array;
    current_solution_replicated_array.ReplicatePetscVector(currentSolution);
    
    // Get our ownership range
    int lo, hi;
    VecGetOwnershipRange(currentSolution, &lo, &hi);
    
    // Get an iterator over the elements of the mesh
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter
    = mpMesh->GetElementIteratorBegin();
    
    // Assume all elements have the same number of nodes...
    const int num_nodes = (*iter)->GetNumNodes();
    
    // Iterate over all elements, summing the contribution of each to the residual
    while (iter != mpMesh->GetElementIteratorEnd())
    {
        const Element<ELEMENT_DIM, SPACE_DIM> &element = **iter;
        
        
        // Ui contains the values of the current solution at the nodes of this element
        c_vector<double, ELEMENT_DIM+1> Ui;
        for (int i=0; i<num_nodes; i++)
        {
            int node_index = element.GetNodeGlobalIndex(i);
            Ui(i) = current_solution_replicated_array[node_index];
        }
        
        // Will contain the contribution of a single element to the residual
        c_vector<double, ELEMENT_DIM+1> b_elem;
        b_elem.clear();
        ComputeResidualOnElement(element, b_elem, mpPde, Ui);
        
        // Update the residual vector for this element
        for (int i=0; i<num_nodes; i++)
        {
            int node_index = element.GetNodeGlobalIndex(i);
            //Make sure it's only done once
            if (lo<=node_index && node_index<hi)
            {
                PetscScalar value = b_elem(i);
                PETSCEXCEPT( VecSetValue(residualVector,node_index,value,ADD_VALUES) );
            }
        }
        iter++;
    }
    
    /*
     * 
     * BOUNDARY CONDITIONS
     * 
     */
    // Add the integrals associated with Neumann boundary conditions to the residual
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator surf_iter
    = mpMesh->GetBoundaryElementIteratorBegin();
    
    if (surf_iter != mpMesh->GetBoundaryElementIteratorEnd())
    {
        const int num_surf_nodes = (*surf_iter)->GetNumNodes();
        
        
        while (surf_iter != mpMesh->GetBoundaryElementIteratorEnd())
        {
            const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& surf_element = **surf_iter;
            
            // UiSurf contains the values of the current solution at the nodes of this surface element
            c_vector<double,ELEMENT_DIM> UiSurf;
            for (int i=0; i<num_surf_nodes; i++)
            {
                int node = surf_element.GetNodeGlobalIndex(i);
                UiSurf(i) = current_solution_replicated_array[node];
            }
            
            /**
             * \todo
             * Check surf_element is in the Neumann surface in an efficient manner.
             */
            if (mpBoundaryConditions->HasNeumannBoundaryCondition(&surf_element))
            {
                c_vector<double, ELEMENT_DIM> b_surf_elem;
                b_surf_elem.clear();
                ComputeResidualOnSurfaceElement(surf_element, b_surf_elem, mpPde, *mpBoundaryConditions, UiSurf);
                
                for (int i=0; i<num_surf_nodes; i++)
                {
                    int node = surf_element.GetNodeGlobalIndex(i);
                    PetscScalar value = b_surf_elem(i);
                    if (lo<=node && node<hi)
                    {
                        PETSCEXCEPT( VecSetValue(residualVector, node, value, ADD_VALUES) );
                    }
                }
            }
            surf_iter++;
        }
    }
    
    // Apply Dirichlet boundary conditions for nonlinear problem
    mpBoundaryConditions->ApplyDirichletToNonlinearResidual(currentSolution, residualVector);
    
    
    VecAssemblyBegin(residualVector);
    VecAssemblyEnd(residualVector);
    
    return 0; // No error
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

/**
 * Function called by PETSc to compute the jacobian matrix given the current solution guess.
 * Calls a method on a SimpleNonlinearEllipticAssembler object to do the calculation.
 *
 * @param snes This is not used by us, but required by PETSc.
 * @param currentSolution The solution guess for the current iteration.
 * @param pGlobalJacobian Pointer to object to fill with the jacobian matrix.
 * @param pPreconditioner This is not used by us, but required by PETSc.
 * @param pMatStructure This is not used by us, but required by PETSc.
 * @param pContext Pointer to a SimpleNonlinearEllipticAssembler object.
 * @return An error code if any PETSc routines fail.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode ComputeJacobianPetsc(SNES snes, Vec currentSolution,
                                    Mat *pGlobalJacobian, Mat *pPreconditioner,
                                    MatStructure *pMatStructure, void *pContext)
{
    // Extract an assembler from the void*
    SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM> *pAssembler =
        (SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>*) pContext;
        
    return pAssembler->ComputeJacobian(currentSolution, pGlobalJacobian);
}

/**
 * Compute the Jacobian matrix given a current guess at the solution.
 * Choose whether to use a numerical or analytical method based on a flag
 * provided by the user.
 *
 * @param currentSolution The solution guess for the current iteration.
 * @param pGlobalJacobian Pointer to object to fill with the jacobian matrix.
 * @return An error code if any PETSc routines fail.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeJacobian(
    const Vec currentSolution, Mat *pGlobalJacobian)
{
    if (mUseAnalyticalJacobian)
    {
        return ComputeJacobianAnalytically(currentSolution, pGlobalJacobian);
    }
    else
    {
        return ComputeJacobianNumerically(currentSolution, pGlobalJacobian);
    }
}




template<int ELEMENT_DIM, int SPACE_DIM>
void SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeJacobianOnElement(
    const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> &rAElem,
    AbstractNonlinearEllipticPde<SPACE_DIM> *pPde,
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
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> f_of_u = pPde->ComputeDiffusionTerm(x,U);
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> f_of_u_prime = pPde->ComputeDiffusionTermPrime(x,U);
        
        //LinearSourceTerm(x)   not needed as it is a constant wrt U_i
        double forcing_term_prime = pPde->ComputeNonlinearSourceTermPrime(x, U);
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

/**
 * Compute the jacobian matrix given the current solution guess.
 *
 * @param currentSolution The solution guess for the current iteration.
 * @param pGlobalJacobian Pointer to object to fill with the jacobian matrix.
 * @return An error code if any PETSc routines fail.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeJacobianAnalytically(
    const Vec currentSolution, Mat *pGlobalJacobian)
{
    // Set all entries of jacobian to 0
    MatZeroEntries(*pGlobalJacobian);
    // Replicate the currentSolution data
    ReplicatableVector current_solution_replicated_array;
    current_solution_replicated_array.ReplicatePetscVector(currentSolution);
    // Get our ownership range
    int lo, hi;
    VecGetOwnershipRange(currentSolution, &lo, &hi);
    // Get an iterator over the elements of the mesh
    typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter =
        mpMesh->GetElementIteratorBegin();
    // Assume all elements have the same number of nodes...
    const int num_nodes = (*iter)->GetNumNodes();
    // Will hold the contribution to the global jacobian from a single element
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> a_elem;
    while (iter != mpMesh->GetElementIteratorEnd())
    {
        const Element<ELEMENT_DIM, SPACE_DIM> &element = **iter;
        a_elem.clear();
        
        // Ui contains the values of the current solution at the nodes of this element
        vector<double> Ui(num_nodes);
        for (int i=0; i<num_nodes; i++)
        {
            int node = element.GetNodeGlobalIndex(i);
            Ui(i) = current_solution_replicated_array[node];
        }
        ComputeJacobianOnElement(element, a_elem, mpPde, Ui);
        // Update global jacobian matrix with this element's contribution
        for (int i=0; i<num_nodes; i++)
        {
            PetscInt node_i = element.GetNodeGlobalIndex(i);
            if (lo<=node_i && node_i<hi)
            {
                for (int j=0; j<num_nodes; j++)
                {
                    PetscInt node_j = element.GetNodeGlobalIndex(j);
                    PetscScalar value = a_elem(i,j);
                    MatSetValue(*pGlobalJacobian, node_i, node_j, value, ADD_VALUES);
                }
            }
        }
        iter++;
    }
    MatAssemblyBegin(*pGlobalJacobian, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*pGlobalJacobian, MAT_FLUSH_ASSEMBLY);
    // Apply dirichlet boundary conditions
    mpBoundaryConditions->ApplyDirichletToNonlinearJacobian(*pGlobalJacobian);
    MatAssemblyBegin(*pGlobalJacobian, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*pGlobalJacobian, MAT_FINAL_ASSEMBLY);
    
    return 0; // No error
}



/**
 * Computes the Jacobian numerically i.e. an approximation, using partial derivatives.
 *
 * @param input Indepedent variable, u in f(u), for example
 * @param *pJacobian A pointer to the Jacobian matrix
 *
 * @return An error code if any PETSc routines fail.
 */
template <int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode SimpleNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>::ComputeJacobianNumerically(
    const Vec input, Mat *pJacobian)
{

    int num_nodes = mpMesh->GetNumNodes();
    
    // Set up working vectors
    Vec residual;
    Vec perturbed_residual;
    Vec result;
    
    VecCreate(PETSC_COMM_WORLD, &residual);
    VecCreate(PETSC_COMM_WORLD, &result);
    VecCreate(PETSC_COMM_WORLD, &perturbed_residual);
    
    VecSetSizes(residual,          PETSC_DECIDE,num_nodes);
    VecSetSizes(result,            PETSC_DECIDE,num_nodes);
    VecSetSizes(perturbed_residual,PETSC_DECIDE,num_nodes);
    
    VecSetFromOptions(residual);
    VecSetFromOptions(result);
    VecSetFromOptions(perturbed_residual);
    
    // Copy the input vector; we perturb the copy
    Vec inputcopy;
    PETSCEXCEPT( VecDuplicate(input, &inputcopy) );
    PETSCEXCEPT( VecCopy(input, inputcopy) );
    
    // Compute the current residual
    ComputeResidual(input, residual);
    
    // Amount to perturb each input element by
    double h = 0.00001;
    PetscScalar subtract = -1;
    PetscScalar one_over_h = 1.0/h;
    
    int lo, hi;
    VecGetOwnershipRange(inputcopy, &lo, &hi);
    // Iterate over entries in the input vector.
    for (int global_index_outer = 0; global_index_outer < num_nodes; global_index_outer++)
    {
        //Only perturb if we own it
        if (lo<=global_index_outer && global_index_outer<hi)
        {
            PETSCEXCEPT( VecSetValue(inputcopy, global_index_outer,h, ADD_VALUES) );
        }
        ComputeResidual(inputcopy, perturbed_residual);
        
        // result = (perturbed_residual - residual) / h
#if (PETSC_VERSION_MINOR == 2) //Old API
        PETSCEXCEPT( VecWAXPY(&subtract, residual, perturbed_residual, result) );
        PETSCEXCEPT( VecScale(&one_over_h, result) );
#else
        PETSCEXCEPT( VecWAXPY(result, subtract, residual, perturbed_residual) );
        PETSCEXCEPT( VecScale(result, one_over_h) );
#endif
        
        double *p_result;
        PETSCEXCEPT( VecGetArray(result, &p_result) );
        for (int global_index=lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;
            PETSCEXCEPT( MatSetValue(*pJacobian, global_index, global_index_outer,
                                     p_result[local_index], INSERT_VALUES) );
        }
        PETSCEXCEPT( VecRestoreArray(result, &p_result) );
        
        if (lo<=global_index_outer && global_index_outer<hi)
        {
            PETSCEXCEPT( VecSetValue(inputcopy, global_index_outer, -h, ADD_VALUES) );
        }
    }
    
    VecDestroy(residual);
    VecDestroy(perturbed_residual);
    VecDestroy(result);
    VecDestroy(inputcopy);
    
    MatAssemblyBegin(*pJacobian, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*pJacobian, MAT_FINAL_ASSEMBLY);
    
    return 0; // No error
}

#endif  // _SIMPLENONLINEARELLIPTICASSEMBLER_HPP_
