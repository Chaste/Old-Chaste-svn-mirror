#ifndef _ABSTRACTNONLINEARELLIPTICASSEMBLER_HPP_
#define _ABSTRACTNONLINEARELLIPTICASSEMBLER_HPP_

/**
 * Abstract superclass for classes that assemble and solve the nonlinear system
 * for a nonlinear elliptic PDE.
 */

#include <petscsnes.h>
#include <petscvec.h>
#include <vector>

#include "AbstractAssembler.hpp"
#include "AbstractNonlinearEllipticPde.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "GaussianQuadratureRule.hpp"
#include "SimpleNonlinearSolver.hpp"
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
PetscErrorCode AssembleResidualPetsc(SNES snes, Vec currentGuess, Vec residualVector,
                                     void *pContext);
                                    
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode AssembleJacobianPetsc(SNES snes, Vec currentGuess,
                                     Mat *pGlobalJacobian, Mat *pPreconditioner,
                                     MatStructure *pMatStructure, void *pContext);
                                    



template<int ELEMENT_DIM, int SPACE_DIM>
class AbstractNonlinearEllipticAssembler : public AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>
{
protected:
    //AbstractNonlinearEllipticPde<SPACE_DIM>* mpPde;
    AbstractNonlinearSolver* mpSolver;
    bool mUseAnalyticalJacobian;
    
    //====================================================================================
    //====================================================================================
    //
    //            pure methods which need to be implemented in the derived class
    //
    //====================================================================================
    //====================================================================================
        
    virtual void AssembleResidualOnElement( const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                                            c_vector<double, ELEMENT_DIM+1> &rBElem,
                                            c_vector<double, ELEMENT_DIM+1> Ui
                                          )=0;


    virtual void AssembleResidualOnSurfaceElement( const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
                                                   c_vector<double, ELEMENT_DIM> &rBsubElem,                                          
                                                   vector<double> Ui
                                                 )=0;


    virtual void AssembleJacobianOnElement( const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                                            c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> &rAElem,
                                            vector<double> Ui
                                          )=0;
        
    


    //====================================================================================
    //====================================================================================
    //
    //                              Compute residual 
    //
    //====================================================================================
    //====================================================================================
public :
    /**
     * Compute the residual vector given the current solution guess.
     *
     * @param currentSolution The solution guess for the current iteration.
     * @param residualVector We fill this with the residual vector.
     * @return An error code if any PETSc routines fail.
     * 
     * This method is called indirectly by the PETSc iterative solvers, so must be 
     * public.
     * 
     */
    PetscErrorCode AssembleResidual(const Vec currentSolution, Vec residualVector)
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
          = this->mpMesh->GetElementIteratorBegin();
        
        // Assume all elements have the same number of nodes...
        const int num_nodes = (*iter)->GetNumNodes();
        
        // Iterate over all elements, summing the contribution of each to the residual
        while (iter != this->mpMesh->GetElementIteratorEnd())
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
            AssembleResidualOnElement(element, b_elem, Ui);
            
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
        = this->mpMesh->GetBoundaryElementIteratorBegin();
        
        if (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
        {
            const int num_surf_nodes = (*surf_iter)->GetNumNodes();
            
            
            while (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
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
                if (this->mpBoundaryConditions->HasNeumannBoundaryCondition(&surf_element))
                {
                    c_vector<double, ELEMENT_DIM> b_surf_elem;
                    b_surf_elem.clear();
                    AssembleResidualOnSurfaceElement(surf_element, b_surf_elem, UiSurf);
                    
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
        this->mpBoundaryConditions->ApplyDirichletToNonlinearResidual(currentSolution, residualVector);
        
        
        VecAssemblyBegin(residualVector);
        VecAssemblyEnd(residualVector);
        
        return 0; // No error
    }
    
    
    
    
    
    //====================================================================================
    //====================================================================================
    //
    //                           Compute Jacobian methods 
    //   
    // AssembleJacobian() is called (via the global function AssembleJacobianPetc by 
    // Petsc), which then calls either AssembleJacobianNumerically() or 
    // AssembleJacobianAnalytically(). If the latter is called AssembleJacobianOnElement()
    // will be called
    //
    //====================================================================================
    //====================================================================================
    
    
    /**
     * Compute the Jacobian matrix given a current guess at the solution.
     * Choose whether to use a numerical or analytical method based on a flag
     * provided by the user.
     *
     * @param currentSolution The solution guess for the current iteration.
     * @param pGlobalJacobian Pointer to object to fill with the jacobian matrix.
     * @return An error code if any PETSc routines fail.
     * 
     * This method is called indirectly by the PETSc iterative solvers, so must be 
     * public.
     * 
     */
    PetscErrorCode AssembleJacobian(const Vec currentSolution, Mat *pGlobalJacobian)
    {
        if (this->mUseAnalyticalJacobian)
        {
            return AssembleJacobianAnalytically(currentSolution, pGlobalJacobian);
        }
        else
        {
            return AssembleJacobianNumerically(currentSolution, pGlobalJacobian);
        }
    }
    
protected:  
    /**
     * Computes the Jacobian numerically i.e. an approximation, using partial derivatives.
     *
     * @param currentGuess Indepedent variable, u in f(u), for example
     * @param pJacobian A pointer to the Jacobian matrix
     *
     * @return An error code if any PETSc routines fail.
     */
    PetscErrorCode AssembleJacobianNumerically(const Vec currentGuess, Mat *pJacobian)
    {
    
        int num_nodes = this->mpMesh->GetNumNodes();
        
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
        
        // Copy the currentGuess vector; we perturb the copy
        Vec current_guess_copy;
        PETSCEXCEPT( VecDuplicate(currentGuess, &current_guess_copy) );
        PETSCEXCEPT( VecCopy(currentGuess, current_guess_copy) );
        
        // Compute the current residual
        AssembleResidual(currentGuess, residual);
        
        // Amount to perturb each input element by
        double h = 0.00001;
        PetscScalar subtract = -1;
        PetscScalar one_over_h = 1.0/h;
        
        int lo, hi;
        VecGetOwnershipRange(current_guess_copy, &lo, &hi);
        // Iterate over entries in the input vector.
        for (int global_index_outer = 0; global_index_outer < num_nodes; global_index_outer++)
        {
            //Only perturb if we own it
            if (lo<=global_index_outer && global_index_outer<hi)
            {
                PETSCEXCEPT( VecSetValue(current_guess_copy, global_index_outer,h, ADD_VALUES) );
            }
            AssembleResidual(current_guess_copy, perturbed_residual);
            
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
                PETSCEXCEPT( VecSetValue(current_guess_copy, global_index_outer, -h, ADD_VALUES) );
            }
        }
        
        VecDestroy(residual);
        VecDestroy(perturbed_residual);
        VecDestroy(result);
        VecDestroy(current_guess_copy);
        
        MatAssemblyBegin(*pJacobian, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(*pJacobian, MAT_FINAL_ASSEMBLY);
        
        return 0; // No error
    }


    
    
    /**
     * Compute the jacobian matrix given the current solution guess.
     *
     * @param currentGuess The solution guess for the current iteration.
     * @param pGlobalJacobian Pointer to object to fill with the jacobian matrix.
     * @return An error code if any PETSc routines fail.
     */
    PetscErrorCode AssembleJacobianAnalytically(const Vec currentGuess, Mat *pGlobalJacobian)
    {
        // Set all entries of jacobian to 0
        MatZeroEntries(*pGlobalJacobian);
        // Replicate the currentGuess data
        ReplicatableVector current_solution_replicated_array;
        current_solution_replicated_array.ReplicatePetscVector(currentGuess);
        // Get our ownership range
        int lo, hi;
        VecGetOwnershipRange(currentGuess, &lo, &hi);
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter =
            this->mpMesh->GetElementIteratorBegin();
        // Assume all elements have the same number of nodes...
        const int num_nodes = (*iter)->GetNumNodes();
        // Will hold the contribution to the global jacobian from a single element
        c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> a_elem;
        while (iter != this->mpMesh->GetElementIteratorEnd())
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
            AssembleJacobianOnElement(element, a_elem, Ui);
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
        this->mpBoundaryConditions->ApplyDirichletToNonlinearJacobian(*pGlobalJacobian);
        MatAssemblyBegin(*pGlobalJacobian, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(*pGlobalJacobian, MAT_FINAL_ASSEMBLY);
        
        return 0; // No error
    }
    

    
public:


    /**
    * Constructors just call the base class versions.
    */
    AbstractNonlinearEllipticAssembler(int numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>(numQuadPoints)
    {
        mpSolver = new SimpleNonlinearSolver;
    }
    
    AbstractNonlinearEllipticAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                       AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                       int numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>(pBasisFunction, pSurfaceBasisFunction, numQuadPoints)
    {
        mpSolver = new SimpleNonlinearSolver;
    }

    ~AbstractNonlinearEllipticAssembler()
    {
        delete mpSolver;
    }

    
   /**
     * Assemble and solve the system for a nonlinear elliptic PDE.
     *
     * @param initialGuess An initial guess for the iterative solver
     * @param UseAnalyticalJacobian Set to true to use an analytically calculated
     *     jacobian matrix rather than a numerically approximated one.
     * @return A PETSc vector giving the solution at each mesh node.
     */
    Vec Solve(Vec initialGuess, bool UseAnalyticalJacobian = false)
    {
        mUseAnalyticalJacobian = UseAnalyticalJacobian;
        
        Vec residual;
        VecDuplicate(initialGuess, &residual);
        
        Vec answer = this->mpSolver->Solve( &AssembleResidualPetsc<ELEMENT_DIM, SPACE_DIM>,
                                            &AssembleJacobianPetsc<ELEMENT_DIM, SPACE_DIM>, 
                                            residual, 
                                            initialGuess, 
                                            this);
                                    
        VecDestroy(residual);
        
        return answer;
    }
                                                             
                                                                 
    void SetNonlinearSolver(AbstractNonlinearSolver* pNonlinearSolver)
    {                           
        delete mpSolver;
        mpSolver = pNonlinearSolver;
    }
                               
};














//========================================================================================//
//========================================================================================//
//
//                Global functions called by the Petsc nonlinear solver
//
//========================================================================================//
//========================================================================================//



/**
 * Function called by PETSc to compute the residual vector given the current solution guess.
 * Calls a method on a SimpleNonlinearEllipticAssembler object to do the calculation.
 *
 * @param snes This is not used by us, but required by PETSc.
 * @param currentGuess The solution guess for the current iteration.
 * @param residualVector We fill this with the residual vector.
 * @param pContext Pointer to a SimpleNonlinearEllipticAssembler object.
 * @return An error code if any PETSc routines fail.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode AssembleResidualPetsc(SNES snes, Vec currentGuess,
                                     Vec residualVector, void *pContext)
{
    AbstractNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM> *pAssembler =
        (AbstractNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>*) pContext;
        
    return pAssembler->AssembleResidual(currentGuess, residualVector);
}



/**
 * Function called by PETSc to compute the jacobian matrix given the current solution guess.
 * Calls a method on a SimpleNonlinearEllipticAssembler object to do the calculation.
 *
 * @param snes This is not used by us, but required by PETSc.
 * @param currentGuess The solution guess for the current iteration.
 * @param pGlobalJacobian Pointer to object to fill with the jacobian matrix.
 * @param pPreconditioner This is not used by us, but required by PETSc.
 * @param pMatStructure This is not used by us, but required by PETSc.
 * @param pContext Pointer to a SimpleNonlinearEllipticAssembler object.
 * @return An error code if any PETSc routines fail.
 */
template<int ELEMENT_DIM, int SPACE_DIM>
PetscErrorCode AssembleJacobianPetsc(SNES snes, Vec currentGuess,
                                    Mat *pGlobalJacobian, Mat *pPreconditioner,
                                    MatStructure *pMatStructure, void *pContext)
{
    // Extract an assembler from the void*
    AbstractNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM> *pAssembler =
        (AbstractNonlinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM>*) pContext;
        
    return pAssembler->AssembleJacobian(currentGuess, pGlobalJacobian);
}


#endif //_ABSTRACTNONLINEARELLIPTICASSEMBLER_HPP_
