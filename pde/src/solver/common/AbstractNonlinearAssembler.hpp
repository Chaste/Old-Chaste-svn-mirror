#ifndef _ABSTRACTNONLINEARASSEMBLER_HPP_
#define _ABSTRACTNONLINEARASSEMBLER_HPP_

/**
 * Abstract superclass for classes that assemble and solve the nonlinear system
 * for a nonlinear elliptic PDE.
 */

#include <petscsnes.h>
#include <petscvec.h>
#include <vector>

#include "AbstractStaticAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "SimplePetscNonlinearSolver.hpp"
#include "RandomNumberGenerator.hpp"
#include "PetscTools.hpp"



/*
 * Since we need to pass function pointers to the PETSc SNES routines, we can't
 * make these functions below methods. This is a pain, since it also means we
 * need to pass round a pointer to our assembler object as the void *pContext,
 * and cast it within the function to access data members.
 *
 * All the functions are defined as stubs which call methods on *pContext.
 *
 * Note: these are global functions, hence the need for long names to avoid
 * potential conflicting names later
 *
 * [The implementations are at the bottom of this file]
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PetscErrorCode AbstractNonlinearAssembler_AssembleResidual(SNES snes,
        Vec currentGuess,
        Vec residualVector,
        void *pContext);
        
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PetscErrorCode AbstractNonlinearAssembler_AssembleJacobian(SNES snes,
        Vec currentGuess,
        Mat *pGlobalJacobian,
        Mat *pPreconditioner,
        MatStructure *pMatStructure,
        void *pContext);
        
        
        
/**
 *  AbstractNonlinearAssembler
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractNonlinearAssembler : public AbstractStaticAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
private:
    Vec mInitialGuess;
protected:
    AbstractNonlinearSolver* mpSolver;
    bool mWeAllocatedSolverMemory;
    
    bool mUseAnalyticalJacobian;
    
    /**
     * Apply Dirichlet boundary conditions to either the residual or jacobian.
     */
    void ApplyDirichletConditions(Vec currentGuess, bool applyToMatrix)
    {
        Vec& residual = this->mpLinearSystem->rGetRhsVector();
        Mat& jacobian = this->mpLinearSystem->rGetLhsMatrix();
        assert((jacobian && applyToMatrix) || (!jacobian && !applyToMatrix));
        
        if (residual)
        {
            this->mpBoundaryConditions->ApplyDirichletToNonlinearResidual(currentGuess, residual);
        }
        if (jacobian)
        {
            this->mpBoundaryConditions->ApplyDirichletToNonlinearJacobian(jacobian);
        }
    }
    
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
        this->PrepareForSolve();
        delete this->mpLinearSystem;
        this->mpLinearSystem = new LinearSystem(residualVector, NULL);
        this->AssembleSystem(true, false, currentGuess, 0.0);
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
            delete this->mpLinearSystem;
            this->mpLinearSystem = new LinearSystem(NULL, *pGlobalJacobian);
            this->AssembleSystem(false, true, currentGuess, 0.0);
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
        unsigned num_nodes = PROBLEM_DIM*this->mpMesh->GetNumNodes();
        
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
        
        PetscInt ilo, ihi;
        VecGetOwnershipRange(current_guess_copy, &ilo, &ihi);
        unsigned lo=ilo;
        unsigned hi=ihi;
        
        // Iterate over entries in the input vector.
        for (unsigned global_index_outer = 0; global_index_outer < num_nodes; global_index_outer++)
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
            for (unsigned global_index=lo; global_index < hi; global_index++)
            {
                unsigned local_index = global_index - lo;
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
    
    virtual void AssembleSystem(bool assembleVector, bool assembleMatrix,
                                Vec currentGuess=NULL, double currentTime=0.0)
    {
        // If the problem is nonlinear the currentGuess MUST be specifed
        assert( currentGuess );
        AbstractStaticAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AssembleSystem(
            assembleVector, assembleMatrix, currentGuess, currentTime);
    }
    
    bool ProblemIsNonlinear()
    {
        return true;
    }

    /**
     * No separate initialisation is needed in the nonlinear case; PrepareForSolve does
     * everything.  We just check the size of the initial guess.
     */
    void InitialiseForSolve(Vec initialGuess)
    {
        // Check size of initial guess is correct
        PetscInt size_of_init_guess;
        VecGetSize(initialGuess, &size_of_init_guess);
        PetscInt problem_size = PROBLEM_DIM * this->mpMesh->GetNumNodes();
        if (size_of_init_guess != problem_size)
        {
            std::stringstream error_message;
            error_message << "Size of initial guess vector, " << size_of_init_guess
                          << ", does not match size of problem, " << problem_size;
            EXCEPTION(error_message.str());
        }
    }
    
    /**
     * Perform the work of a single solve, but without any initialisation.
     * 
     * @param currentSolutionOrGuess  either the current solution (dynamic problem) or
     *     initial guess (static problem).  MUST be provided. 
     * @param currentTime  for a dynamic problem, the current time
     * @return the solution vector
     */
    Vec StaticSolve(Vec currentSolutionOrGuess=NULL,
                    double currentTime=0.0,
                    bool assembleMatrix=true)
    {
        assert(this->mpBoundaryConditions!=NULL); // no flagged mesh to worry about here!
        assert(currentSolutionOrGuess!=NULL);
        assert(assembleMatrix); ///\todo do something sensible if assembleMatrix is false.
        
        // run the solver, telling it which global functions to call in order to assemble
        // the residual or jacobian
        Vec answer = this->mpSolver->Solve(&AbstractNonlinearAssembler_AssembleResidual<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>,
                                           &AbstractNonlinearAssembler_AssembleJacobian<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>,
                                           currentSolutionOrGuess,
                                           this);
        return answer;
    }
    
public:
    /**
     * Constructors just call the base class versions.
     */
    AbstractNonlinearAssembler(unsigned numQuadPoints = 2) :
            AbstractStaticAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(numQuadPoints)
    {
        mpSolver = new SimplePetscNonlinearSolver;
        mWeAllocatedSolverMemory = true;
        
        mUseAnalyticalJacobian = false;
        
        mInitialGuess = NULL;
    }
    
    ~AbstractNonlinearAssembler()
    {
        if (mWeAllocatedSolverMemory)
        {
            delete mpSolver;
        }
    }
        
    /**
     * Set whether to use an analytical jacobian.  This is provided for use when
     * solving dynamic nonlinear problems; when solving static problems there is
     * an argument to the Solve method which specifies this property, and overrides
     * any user call to this method.
     * 
     * If this method is not called the default is false i.e. numerical jacobian.
     */
    void SetUseAnalyticalJacobian(bool useAnalyticalJacobian)
    {
        mUseAnalyticalJacobian = useAnalyticalJacobian;
    }
    
    /**
      * Assemble and solve the system for a nonlinear elliptic PDE.
      *
      * @param initialGuess An initial guess for the iterative solver
      * @param useAnalyticalJacobian Set to true to use an analytically calculated
      *     jacobian matrix rather than a numerically approximated one.
      * @return A PETSc vector giving the solution at each mesh node.
      */
    virtual Vec Solve(Vec initialGuess, bool useAnalyticalJacobian=false)
    {
        assert(initialGuess!=NULL);
        SetUseAnalyticalJacobian(useAnalyticalJacobian);
        
        this->PrepareForSolve();
        InitialiseForSolve(initialGuess);
        
        return StaticSolve(initialGuess);
    }
    
    /**
     *  SetNonlinearSolver - by default a SimplePetscNonlinearSolver is created
     *  and used in this class, this method can be called to use a different
     *  AbstractNonlinearSolver
     */
    void SetNonlinearSolver(AbstractNonlinearSolver* pNonlinearSolver)
    {
        if (mWeAllocatedSolverMemory)
        {
            delete mpSolver;
        }
        mpSolver = pNonlinearSolver;
        mWeAllocatedSolverMemory = false;
    }
    
    /**
     *  A helpful method for creating an initial guess vector
     */
    Vec CreateConstantInitialGuess(double value)
    {
        assert(this->mpMesh!=NULL);
        unsigned size = PROBLEM_DIM * this->mpMesh->GetNumNodes();
        return PetscTools::CreateVec(size, value);
    }
    
    /**
     *  VerifyJacobian
     * 
     *  A helper method for use when writing concrete assemblers. Once the user has calculated
     *  (on paper) the weak form and the form of the ComputeMatrixTerm method, they can check 
     *  whether the analytic Jacobian matches the numerical Jacobian (which only calls 
     *  ComputeVectorTerm and ComputeVectorSurfaceTerm) to verify the correctness of the code.
     * 
     *  @param tol A tolerance which defaults to 1e-5
     *  @param print  Whether to print the different matrix J_numerical-J_analytical
     *  @return true if the componentwise difference between the matrices is less than 
     *    the tolerance, false otherwise.
     * 
     *  This method should NOT be run in simulations - it is only to verify the correctness 
     *  of the concrete assembler code.
     */
    bool VerifyJacobian(double tol=1e-4, bool print=false)
    {
        unsigned size = PROBLEM_DIM * this->mpMesh->GetNumNodes();
        
        Vec initial_guess;
        VecCreate(PETSC_COMM_WORLD, &initial_guess);
        VecSetSizes(initial_guess, PETSC_DECIDE, size);
        VecSetFromOptions(initial_guess);
        
        
        for (unsigned i=0; i<size ; i++)
        {
            VecSetValue(initial_guess, i, 0.0, INSERT_VALUES);
        }
        
        VecAssemblyBegin(initial_guess);
        VecAssemblyEnd(initial_guess);
        
        
        Mat analytic_jacobian; //Jacobian Matrix
        Mat numerical_jacobian; //Jacobian Matrix
        
#if (PETSC_VERSION_MINOR == 2) //Old API
        MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,size,size,&analytic_jacobian);
        MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,size,size,&numerical_jacobian);
#else
        MatCreate(PETSC_COMM_WORLD,&analytic_jacobian);
        MatSetSizes(analytic_jacobian,PETSC_DECIDE,PETSC_DECIDE,size,size);
        MatCreate(PETSC_COMM_WORLD,&numerical_jacobian);
        MatSetSizes(numerical_jacobian,PETSC_DECIDE,PETSC_DECIDE,size,size);
#endif
        MatSetFromOptions(analytic_jacobian);
        MatSetFromOptions(numerical_jacobian);
        
        mUseAnalyticalJacobian = true;
        AssembleJacobian(initial_guess, &analytic_jacobian);
        
        mUseAnalyticalJacobian = false;
        AssembleJacobian(initial_guess, &numerical_jacobian);
        
        bool all_less_than_tol = true;
        
        if (print)
        {
            std::cout << "Difference between numerical and analyical Jacobians:\n\n";
        }
        for (unsigned i=0; i<size; i++)
        {
            for (unsigned j=0; j<size; j++)
            {
                double val_a[1];
                double val_n[1];
                PetscInt row[1];
                PetscInt col[1];
                row[0] = i;
                col[0] = j;
                
                MatGetValues(numerical_jacobian,1,row,1,col,val_n);
                MatGetValues(analytic_jacobian,1,row,1,col,val_a);
                
                if (print)
                {
                    std::cout << val_n[0] - val_a[0]<< " ";
                }
                
                if (fabs(val_n[0]-val_a[0]) > tol)
                {
#define COVERAGE_IGNORE // would have to write a bad concrete assembler class just to cover this line
                    all_less_than_tol = false;
#undef COVERAGE_IGNORE
                }
            }
            if (print)
            {
                std::cout << "\n" <<std::flush;
            }
        }
        std::cout << std::flush;
        
        MatDestroy(numerical_jacobian);
        MatDestroy(analytic_jacobian);
        VecDestroy(initial_guess);
        
        return all_less_than_tol;
    }

}; //end of class AbstractNonlinearAssembler






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
 *
 * Note: this is a global function, hence the need for a long name to avoid
 * potential conflicting names later
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PetscErrorCode AbstractNonlinearAssembler_AssembleResidual(SNES snes, Vec currentGuess,
        Vec residualVector, void *pContext)
{
    // Extract an assembler from the void*
    AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> *pAssembler =
        (AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*) pContext;
        
    PetscErrorCode ierr = pAssembler->AssembleResidual(currentGuess, residualVector);
    
    //double two_norm;
    //VecNorm(residualVector, NORM_2, &two_norm);
    //std::cout << "||residual|| = " << two_norm << "\n";
    
    return ierr;
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
 *
 * Note: this is a global function, hence the need a long name to avoid
 * potential conflicting names later
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PetscErrorCode AbstractNonlinearAssembler_AssembleJacobian(SNES snes, Vec currentGuess,
        Mat *pGlobalJacobian, Mat *pPreconditioner,
        MatStructure *pMatStructure, void *pContext)
{
    //std::cout << "begin assemble jacobian\n";
    
    // Extract an assembler from the void*
    AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> *pAssembler =
        (AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*) pContext;
        
    PetscErrorCode ierr = pAssembler->AssembleJacobian(currentGuess, pGlobalJacobian);
    //std::cout << "end assemble jacobian\n";
    
    return ierr;
}



#endif //_ABSTRACTNONLINEARASSEMBLER_HPP_
