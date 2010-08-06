/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef _ABSTRACTNONLINEARASSEMBLER_HPP_
#define _ABSTRACTNONLINEARASSEMBLER_HPP_

/***  THIS CLASS WILL BE DELETED NEXT WEEK ***/
#define COVERAGE_IGNORE



/**
 * Abstract superclass for classes that assemble and solve the nonlinear system
 * for a nonlinear elliptic PDE.
 */



#include <vector>

#include "AbstractStaticAssembler.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "SimplePetscNonlinearSolver.hpp"
#include "RandomNumberGenerator.hpp"
#include "PetscTools.hpp"


/*
 * Since we need to pass function pointers to the PETSc SNES routines, we can't
 * make these functions below methods. This is a pain, since it also means we
 * need to pass round a pointer to our assembler object as the void* pContext,
 * and cast it within the function to access data members.
 *
 * All the functions are defined as stubs which call methods on* pContext.
 *
 * Note: these are global functions, hence the need for long names to avoid
 * potential conflicting names later
 *
 * [The implementations are at the bottom of this file]
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
PetscErrorCode AbstractNonlinearAssembler_AssembleResidual(SNES snes,
        Vec currentGuess,
        Vec residualVector,
        void* pContext);

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
PetscErrorCode AbstractNonlinearAssembler_AssembleJacobian(SNES snes,
        Vec currentGuess,
        Mat* pGlobalJacobian,
        Mat* pPreconditioner,
        MatStructure* pMatStructure,
        void* pContext);


/**
 *  AbstractNonlinearAssembler
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
class AbstractNonlinearAssembler : public AbstractStaticAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, true, CONCRETE>
{
    typedef AbstractStaticAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, true, CONCRETE> BaseClassType; /**< Base class type (to save typing). */

private:

    /** PETSc vector storing an initial guess for the solution. */
    Vec mInitialGuess;
    unsigned mFillSize; /**< Calculated in InitialiseForSolve and used to communicate connectivity information from the mesh to the nonlinear solver*/

protected:

    /** The nonlinear solver. */
    AbstractNonlinearSolver* mpSolver;

    /** Whether memory has been allocated for the solver. */
    bool mWeAllocatedSolverMemory;

    /** Whether to use an analytical expression for the Jacobian. */
    bool mUseAnalyticalJacobian;

    /**
     * Apply Dirichlet boundary conditions to either the residual or Jacobian.
     *
     * @param currentGuess  the solution guess for the current iteration
     * @param applyToMatrix  whether to apply the boundary conditions to the Jacobian matrix
     */
    void ApplyDirichletConditions(Vec currentGuess, bool applyToMatrix);

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
    PetscErrorCode AssembleResidual(const Vec currentGuess, Vec residualVector);

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
    PetscErrorCode AssembleJacobian(const Vec currentGuess, Mat* pGlobalJacobian);

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
    PetscErrorCode AssembleJacobianNumerically(const Vec currentGuess, Mat* pJacobian);

    /**
     *  AssembleSystem - the major method for all assemblers
     *
     *  Assemble the linear system for a linear PDE, or the residual or Jacobian for
     *  nonlinear PDEs. Loops over each element (and each each surface element if
     *  there are non-zero Neumann boundary conditions), calls AssembleOnElement()
     *  and adds the contribution to the linear system.
     *
     *  Takes in current solution and time if necessary but only used if the problem
     *  is a dynamic one. This method uses PROBLEM_DIM and can assemble linear systems
     *  for any number of unknown variables.
     *
     *  @param assembleVector  Whether to assemble the RHS vector of the linear system
     *     (i.e. the residual vector for nonlinear problems).
     *  @param assembleMatrix  Whether to assemble the LHS matrix of the linear system
     *     (i.e. the jacobian matrix for nonlinear problems).
     *  @param currentGuess The current solution in a linear dynamic problem,
     *     or the current guess in a nonlinear problem. Should be NULL for linear static
     *     problems. Defaults to NULL.
     *  @param currentTime The current time for dynamic problems. Not used in static
     *     problems. Defaults to 0.0.
     */
    virtual void AssembleSystem(bool assembleVector, bool assembleMatrix,
                                Vec currentGuess=NULL, double currentTime=0.0);

    /**
     * Whether grad_u should be calculated
     */
    bool ProblemIsNonlinear();

    /**
     * No separate initialisation is needed in the nonlinear case; PrepareForSolve does
     * everything.  We just check the size of the initial guess.
     *
     * @param initialGuess an initial guess
     */
    void InitialiseForSolve(Vec initialGuess);

    /**
     * Perform the work of a single solve, but without any initialisation.
     *
     * @param currentSolutionOrGuess  either the current solution (dynamic problem) or
     *     initial guess (static problem).  MUST be provided.
     * @param currentTime  for a dynamic problem, the current time
     * @param assembleMatrix  Whether to assemble the LHS matrix of the linear system
     *     (i.e. the jacobian matrix for nonlinear problems).
     * @return the solution vector
     */
    Vec StaticSolve(Vec currentSolutionOrGuess=NULL,
                    double currentTime=0.0,
                    bool assembleMatrix=true);

public:

    /**
     * Constructors just call the base class versions.
     *
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    AbstractNonlinearAssembler(unsigned numQuadPoints = 2);

    /**
     * Destructor.
     */
    ~AbstractNonlinearAssembler();

    /**
     * Set whether to use an analytical jacobian.  This is provided for use when
     * solving dynamic nonlinear problems; when solving static problems there is
     * an argument to the Solve method which specifies this property, and overrides
     * any user call to this method.
     *
     * If this method is not called the default is false i.e. numerical jacobian.
     *
     * @param useAnalyticalJacobian Set to true to use an analytically calculated
      *     jacobian matrix rather than a numerically approximated one.
     */
    void SetUseAnalyticalJacobian(bool useAnalyticalJacobian);

    /**
      * Assemble and solve the system for a nonlinear elliptic PDE.
      *
      * @param initialGuess An initial guess for the iterative solver
      * @param useAnalyticalJacobian Set to true to use an analytically calculated
      *     jacobian matrix rather than a numerically approximated one.
      * @return A PETSc vector giving the solution at each mesh node.
      */
    virtual Vec Solve(Vec initialGuess, bool useAnalyticalJacobian=false);

    /**
     * SetNonlinearSolver - by default a SimplePetscNonlinearSolver is created
     * and used in this class, this method can be called to use a different
     * AbstractNonlinearSolver.
     *
     * @param pNonlinearSolver a nonlinear solver
     */
    void SetNonlinearSolver(AbstractNonlinearSolver* pNonlinearSolver);

    /**
     * A helpful method for creating an initial guess vector.
     *
     * @param value constant value of the initial guess
     */
    Vec CreateConstantInitialGuess(double value);

    /**
     *  VerifyJacobian
     *
     *  A helper method for use when writing concrete assemblers. Once the user has calculated
     *  (on paper) the weak form and the form of the ComputeMatrixTerm method, they can check
     *  whether the analytic Jacobian matches the numerical Jacobian (which only calls
     *  ComputeVectorTerm and ComputeVectorSurfaceTerm) to verify the correctness of the code.
     *
     *  @param tol A tolerance which defaults to 1e-5
     *  @return true if the componentwise difference between the matrices is less than
     *    the tolerance, false otherwise.
     *
     *  This method should NOT be run in simulations - it is only to verify the correctness
     *  of the concrete assembler code.
     */
    bool VerifyJacobian(double tol=1e-4);

}; //end of class AbstractNonlinearAssembler


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
void AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::ApplyDirichletConditions(Vec currentGuess, bool applyToMatrix)
{
    Vec& residual = this->mpLinearSystem->rGetRhsVector();
    Mat& jacobian = this->mpLinearSystem->rGetLhsMatrix();
    assert((jacobian && applyToMatrix) || (!jacobian && !applyToMatrix));

    if (residual)
    {
        this->mpBoundaryConditions->ApplyDirichletToNonlinearResidual(
            currentGuess, residual, *(this->mpMesh->GetDistributedVectorFactory()));
    }
    if (jacobian)
    {
        this->mpBoundaryConditions->ApplyDirichletToNonlinearJacobian(jacobian);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
PetscErrorCode AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::AssembleResidual(const Vec currentGuess, Vec residualVector)
{
    // call assemble system with the current guess and the residual vector
    // to be assembled
    this->PrepareForSolve();
    delete this->mpLinearSystem;
    this->mpLinearSystem = new LinearSystem(residualVector, (Mat) NULL);
    this->AssembleSystem(true, false, currentGuess, 0.0);
    return 0;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
PetscErrorCode AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::AssembleJacobian(const Vec currentGuess, Mat* pGlobalJacobian)
{
    if (mUseAnalyticalJacobian)
    {
        // call assemble system with the current guess and the jacobian to
        // be assembled
        delete this->mpLinearSystem;
        this->mpLinearSystem = new LinearSystem(NULL,* pGlobalJacobian);
        this->AssembleSystem(false, true, currentGuess, 0.0);
        return 0; // no error
    }
    else
    {
        return AssembleJacobianNumerically(currentGuess, pGlobalJacobian);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
PetscErrorCode AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::AssembleJacobianNumerically(const Vec currentGuess, Mat* pJacobian)
{
    unsigned num_unknowns = PROBLEM_DIM*this->mpMesh->GetNumNodes();

    // Set up working vectors
    Vec residual;
    Vec perturbed_residual;
    Vec result;
    residual=PetscTools::CreateVec(num_unknowns);
    result=PetscTools::CreateVec(num_unknowns);
    perturbed_residual=PetscTools::CreateVec(num_unknowns);
    
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
    unsigned lo = ilo;
    unsigned hi = ihi;

    // Iterate over entries in the input vector.
    for (unsigned global_index_outer = 0; global_index_outer < num_unknowns; global_index_outer++)
    {
        //Only perturb if we own it
        if (lo<=global_index_outer && global_index_outer<hi)
        {
            PETSCEXCEPT( VecSetValue(current_guess_copy, global_index_outer,h, ADD_VALUES) );
        }
        AssembleResidual(current_guess_copy, perturbed_residual);

        // result = (perturbed_residual - residual) / h
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
        PETSCEXCEPT( VecWAXPY(&subtract, residual, perturbed_residual, result) );
        PETSCEXCEPT( VecScale(&one_over_h, result) );
#else
        PETSCEXCEPT( VecWAXPY(result, subtract, residual, perturbed_residual) );
        PETSCEXCEPT( VecScale(result, one_over_h) );
#endif

        double* p_result;
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


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
void AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::AssembleSystem(bool assembleVector, bool assembleMatrix,
                                Vec currentGuess, double currentTime)
{
    // If the problem is nonlinear the currentGuess MUST be specifed
    assert( currentGuess );
    BaseClassType::AssembleSystem(assembleVector, assembleMatrix, currentGuess, currentTime);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
bool AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::ProblemIsNonlinear()
{
    return true;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
void AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::InitialiseForSolve(Vec initialGuess)
{
    // Check size of initial guess is correct
    PetscInt size_of_init_guess;
    VecGetSize(initialGuess, &size_of_init_guess);
    PetscInt problem_size = PROBLEM_DIM * this->mpMesh->GetNumNodes();
    mFillSize = PROBLEM_DIM * this->mpMesh->CalculateMaximumNodeConnectivityPerProcess();
    if (size_of_init_guess != problem_size)
    {
        std::stringstream error_message;
        error_message << "Size of initial guess vector, " << size_of_init_guess
                      << ", does not match size of problem, " << problem_size;
        EXCEPTION(error_message.str());
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
Vec AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::StaticSolve(Vec currentSolutionOrGuess,
                    double currentTime,
                    bool assembleMatrix)
{
    assert(this->mpBoundaryConditions!=NULL);
    assert(currentSolutionOrGuess!=NULL);
    assert(assembleMatrix); ///\todo do something sensible if assembleMatrix is false.

    // run the solver, telling it which global functions to call in order to assemble
    // the residual or jacobian
    Vec answer = this->mpSolver->Solve(&AbstractNonlinearAssembler_AssembleResidual<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>,
                                       &AbstractNonlinearAssembler_AssembleJacobian<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>,
                                       currentSolutionOrGuess,
                                       mFillSize, 
                                       this);
    return answer;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::AbstractNonlinearAssembler(unsigned numQuadPoints)
    : BaseClassType(numQuadPoints),
      mInitialGuess(NULL),
      mWeAllocatedSolverMemory(true),
      mUseAnalyticalJacobian(false)
{
    mpSolver = new SimplePetscNonlinearSolver;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::~AbstractNonlinearAssembler()
{
    if (mWeAllocatedSolverMemory)
    {
        delete mpSolver;
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
void AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::SetUseAnalyticalJacobian(bool useAnalyticalJacobian)
{
    mUseAnalyticalJacobian = useAnalyticalJacobian;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
Vec AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::Solve(Vec initialGuess, bool useAnalyticalJacobian)
{
    assert(initialGuess!=NULL);
    SetUseAnalyticalJacobian(useAnalyticalJacobian);

    this->PrepareForSolve();
    InitialiseForSolve(initialGuess);

    return StaticSolve(initialGuess);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
void AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::SetNonlinearSolver(AbstractNonlinearSolver* pNonlinearSolver)
{
    if (mWeAllocatedSolverMemory)
    {
        delete mpSolver;
    }
    mpSolver = pNonlinearSolver;
    mWeAllocatedSolverMemory = false;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
Vec AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::CreateConstantInitialGuess(double value)
{
    assert(this->mpMesh!=NULL);
    unsigned size = PROBLEM_DIM * this->mpMesh->GetNumNodes();
    return PetscTools::CreateAndSetVec(size, value);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
bool AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>::VerifyJacobian(double tol)
{
    unsigned size = PROBLEM_DIM * this->mpMesh->GetNumNodes();

    Vec initial_guess=PetscTools::CreateAndSetVec(size, 0.0);
    
    Mat analytic_jacobian; //Jacobian Matrix
    Mat numerical_jacobian; //Jacobian Matrix

    PetscTools::SetupMat(analytic_jacobian, size, size, size);
    PetscTools::SetupMat(numerical_jacobian, size, size, size);

    mUseAnalyticalJacobian = true;
    AssembleJacobian(initial_guess, &analytic_jacobian);

    mUseAnalyticalJacobian = false;
    AssembleJacobian(initial_guess, &numerical_jacobian);
    
    double minus_one = -1.0;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    // MatAYPX(*a, X, Y) does  Y = X + a*Y. 
    MatAYPX(&minus_one, analytic_jacobian, numerical_jacobian);
#else
    // MatAYPX( Y, a, X, structure) does Y = a*Y + X. 
    MatAYPX(numerical_jacobian, minus_one,analytic_jacobian,DIFFERENT_NONZERO_PATTERN);
#endif
    double norm;
    MatNorm(numerical_jacobian,NORM_INFINITY,&norm);
    
    MatDestroy(numerical_jacobian);
    MatDestroy(analytic_jacobian);
    VecDestroy(initial_guess);

    return (norm<tol);
}


//========================================================================================//
//========================================================================================//
//                                                                                        //
//                Global functions called by the Petsc nonlinear solver                   //
//                                                                                        //
//========================================================================================////========================================================================================//



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
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
PetscErrorCode AbstractNonlinearAssembler_AssembleResidual(SNES snes, Vec currentGuess,
        Vec residualVector, void* pContext)
{
    // Extract an assembler from the void*
    AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>* pAssembler =
        (AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>*) pContext;

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
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, class CONCRETE>
PetscErrorCode AbstractNonlinearAssembler_AssembleJacobian(SNES snes, Vec currentGuess,
        Mat* pGlobalJacobian, Mat* pPreconditioner,
        MatStructure* pMatStructure, void* pContext)
{
    //std::cout << "begin assemble jacobian\n";

    // Extract an assembler from the void*
    AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>* pAssembler =
        (AbstractNonlinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CONCRETE>*) pContext;

    PetscErrorCode ierr = pAssembler->AssembleJacobian(currentGuess, pGlobalJacobian);
    //std::cout << "end assemble jacobian\n";

    return ierr;
}

#undef COVERAGE_IGNORE


#endif //_ABSTRACTNONLINEARASSEMBLER_HPP_
