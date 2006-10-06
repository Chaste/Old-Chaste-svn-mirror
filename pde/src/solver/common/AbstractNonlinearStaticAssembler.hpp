#ifndef _ABSTRACTNONLINEARSTATICASSEMBLER_HPP_
#define _ABSTRACTNONLINEARSTATICASSEMBLER_HPP_

/**
 * Abstract superclass for classes that assemble and solve the nonlinear system
 * for a nonlinear elliptic PDE.
 */

#include <petscsnes.h>
#include <petscvec.h>
#include <vector>

#include "AbstractLinearAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "SimpleNonlinearSolver.hpp"



/*
 * Since we need to pass function pointers to the PETSc SNES routines, we can't
 * make these functions below methods. This is a pain, since it also means we
 * need to pass round a pointer to our assembler object as the void *pContext,
 * and cast it within the function to access data members.
 *
 * All the functions are defined as stubs which call methods on *pContext.
 * 
 * [The implementations are at the bottom of this file]
 */
template<int ELEMENT_DIM, int SPACE_DIM, int PROBLEM_DIM>
PetscErrorCode AssembleResidualPetsc(SNES snes, Vec currentGuess, Vec residualVector,
                                     void *pContext);
                                    
template<int ELEMENT_DIM, int SPACE_DIM, int PROBLEM_DIM>
PetscErrorCode AssembleJacobianPetsc(SNES snes, Vec currentGuess,
                                     Mat *pGlobalJacobian, Mat *pPreconditioner,
                                     MatStructure *pMatStructure, void *pContext);
                                    


/**
 *  AbstractNonlinearStaticAssembler
 * 
 *  to become AbstractNonlinearStaticAssembler ?
 */
template<int ELEMENT_DIM, int SPACE_DIM, int PROBLEM_DIM>
class AbstractNonlinearStaticAssembler : public AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
protected:
    AbstractNonlinearSolver* mpSolver;
    bool mUseAnalyticalJacobian;
    
public :
    /**
     * Compute the residual vector given the current solution guess.
     *
     * @param currentGuess The solution guess for the current iteration.
     * @param residualVector We fill this with the residual vector.
     * @return An error code if any PETSc routines fail.
     * 
     * NOTE: this method is called indirectly by the PETSc iterative 
     * solvers, so must be public.
     */
    PetscErrorCode AssembleResidual(const Vec currentGuess, Vec residualVector)
    {
        // call assemble system with the current guess and the residual vector
        // to be assembled
        this->AssembleSystem(currentGuess, 0.0, residualVector, NULL);
        return 0;
    }
    
    
    /**
     * Compute the Jacobian matrix given a current guess at the solution.
     * Choose whether to use a numerical or analytical method based on a flag
     * provided by the user (in Solve()).
     *
     * @param currentGuess The solution guess for the current iteration.
     * @param pGlobalJacobian Pointer to object to fill with the jacobian matrix.
     * @return An error code if any PETSc routines fail.
     * 
     * NOTE: this method is called indirectly by the PETSc iterative 
     * solvers, so must be public.
     */
    PetscErrorCode AssembleJacobian(const Vec currentGuess, Mat *pGlobalJacobian)
    {
        if (mUseAnalyticalJacobian)
        {
            // call assemble system with the current guess and the jacobian to
            // be assembled
            this->AssembleSystem(currentGuess,0.0,NULL,pGlobalJacobian);
            return 0; // no error
        }
        else
        {
            return AssembleJacobianNumerically(currentGuess, pGlobalJacobian);
        }
    }
    
protected:  
    /**
     * Computes the Jacobian numerically i.e. an approximation, using numerical
     * partial derivatives.
     *
     * @param currentGuess Independent variable, u in f(u), for example
     * @param pJacobian A pointer to the Jacobian matrix
     *
     * @return An error code if any PETSc routines fail.
     */
    PetscErrorCode AssembleJacobianNumerically(const Vec currentGuess, Mat *pJacobian)
    {
        int num_nodes = PROBLEM_DIM*this->mpMesh->GetNumNodes();
        
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


public:


    /**
    * Constructors just call the base class versions.
    */
    AbstractNonlinearStaticAssembler(int numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(numQuadPoints)
    {
        mpSolver = new SimpleNonlinearSolver;  
        this->mProblemIsLinear = false;
        
        mUseAnalyticalJacobian = true;
    }
    
    AbstractNonlinearStaticAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                                       AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                                       int numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(pBasisFunction, pSurfaceBasisFunction, numQuadPoints)
    {
        mpSolver = new SimpleNonlinearSolver;
        this->mpProblemIsLinear = false;

        mUseAnalyticalJacobian = true;
    }

    ~AbstractNonlinearStaticAssembler()
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
    virtual Vec Solve(Vec initialGuess, bool UseAnalyticalJacobian = false)
    {
        assert(this->mpMesh!=NULL);
        assert(this->mpBoundaryConditions!=NULL);

        // check size of initial guess is correct
        int size_of_init_guess;
        VecGetSize(initialGuess, &size_of_init_guess);
        if(size_of_init_guess!=PROBLEM_DIM * this->mpMesh->GetNumNodes())
        {
            std::stringstream error_message;
            error_message << "Size of initial guess vector, " << size_of_init_guess 
                          << ", does not match size of problem, " 
                          << PROBLEM_DIM * this->mpMesh->GetNumNodes();
            
            EXCEPTION(error_message.str());
        }
        
        
        mUseAnalyticalJacobian = UseAnalyticalJacobian;
        
        Vec residual;
        VecDuplicate(initialGuess, &residual);
        
        Vec answer = this->mpSolver->Solve( &AssembleResidualPetsc<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>,
                                            &AssembleJacobianPetsc<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>, 
                                            residual, 
                                            initialGuess, 
                                            this);
                                    
        VecDestroy(residual);
        
        return answer;
    }
                                                             

    /** 
     *  SetNonlinearSolver - by default a SimpleNonlinearSolver is created
     *  and used in this class, this method can be called to use a different
     *  AbstractNonlinearSolver
     */
    void SetNonlinearSolver(AbstractNonlinearSolver* pNonlinearSolver)
    {                           
        delete mpSolver;
        mpSolver = pNonlinearSolver;
    }                        
    
    
    
    /** 
     *  A helpful method for creating an initial guess vector
     */
    Vec CreateConstantInitialGuess(double value)
    {
        assert(this->mpMesh!=NULL);
        
        int size = PROBLEM_DIM * this->mpMesh->GetNumNodes();

        Vec initial_guess;
        VecCreate(PETSC_COMM_WORLD, &initial_guess);
        VecSetSizes(initial_guess, PETSC_DECIDE, size);
        VecSetFromOptions(initial_guess);
        
#if (PETSC_VERSION_MINOR == 2) //Old API
        VecSet(&value, initial_guess);
#else
        VecSet(initial_guess, value);
#endif
        
        VecAssemblyBegin(initial_guess);
        VecAssemblyEnd(initial_guess);
        return initial_guess;
    }
    
}; //end of class AbstractNonlinearStaticAssembler






//========================================================================================//
//========================================================================================//
//                                                                                        //
//                Global functions called by the Petsc nonlinear solver                   //
//                                                                                        //
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
template<int ELEMENT_DIM, int SPACE_DIM, int PROBLEM_DIM>
PetscErrorCode AssembleResidualPetsc(SNES snes, Vec currentGuess,
                                     Vec residualVector, void *pContext)
{
    // Extract an assembler from the void*
    AbstractNonlinearStaticAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> *pAssembler =
        (AbstractNonlinearStaticAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*) pContext;
        
    double two_norm;
    VecNorm(residualVector, NORM_2, &two_norm);
    std::cout << "||residual|| = " << two_norm << "\n";
        
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
template<int ELEMENT_DIM, int SPACE_DIM, int PROBLEM_DIM>
PetscErrorCode AssembleJacobianPetsc(SNES snes, Vec currentGuess,
                                    Mat *pGlobalJacobian, Mat *pPreconditioner,
                                    MatStructure *pMatStructure, void *pContext)
{
    std::cout << "begin assemble jacobian\n";
    // Extract an assembler from the void*
    AbstractNonlinearStaticAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> *pAssembler =
        (AbstractNonlinearStaticAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*) pContext;

    PetscErrorCode ierr = pAssembler->AssembleJacobian(currentGuess, pGlobalJacobian);
    std::cout << "end assemble jacobian\n";
        
    return ierr;
}


#endif //_ABSTRACTNONLINEARSTATICASSEMBLER_HPP_
