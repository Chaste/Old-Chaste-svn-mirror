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

#ifndef ABSTRACTNONLINEARASSEMBLERSOLVERHYBRID_HPP_
#define ABSTRACTNONLINEARASSEMBLERSOLVERHYBRID_HPP_



#include "BoundaryConditionsContainer.hpp"
#include "AbstractNonlinearEllipticPde.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "PetscTools.hpp"
#include "AbstractFeObjectAssembler.hpp"
#include "LinearSystem.hpp"
#include "SimplePetscNonlinearSolver.hpp"


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////
////
////    Definitions of GLOBAL functions used by Petsc nonlinear solver
////    (implementations at the bottom of this file)
////
////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


/**
 * Function called by PETSc to compute the residual vector given the current solution guess.
 * Calls a method on a AbstractNonlinearAssemblerSolverHybrid object to do the calculation.
 *
 * @param snes This is not used by us, but required by PETSc.
 * @param currentGuess The solution guess for the current iteration.
 * @param residualVector We fill this with the residual vector.
 * @param pContext Pointer to a AbstractNonlinearAssemblerSolverHybrid object.
 * @return Petsc expects this function to return a PetscErrorCode. We always return 0;
 *   exceptions are thrown if there is an error.
 *
 * Note: this is a global function, hence the need for a long name to avoid
 * potential conflicting names
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PetscErrorCode AbstractNonlinearAssemblerSolverHybrid_ComputeResidual(SNES snes,
                                                                      Vec currentGuess,
                                                                      Vec residualVector,
                                                                      void* pContext);



/**
 * Function called by PETSc to compute the jacobian matrix given the current solution guess.
 * Calls a method on a AbstractNonlinearAssemblerSolverHybrid object to do the calculation.
 *
 * @param snes This is not used by us, but required by PETSc.
 * @param currentGuess The solution guess for the current iteration.
 * @param pGlobalJacobian Pointer to object to fill with the jacobian matrix.
 * @param pPreconditioner This is not used by us, but required by PETSc.
 * @param pMatStructure This is not used by us, but required by PETSc.
 * @param pContext Pointer to a AbstractNonlinearAssemblerSolverHybrid object.
 * @return Petsc expects this function to return a PetscErrorCode. We always return 0;
 *   exceptions are thrown if there is an error.
 *
 * Note: this is a global function, hence the need a long name to avoid
 * potential conflicting names
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PetscErrorCode AbstractNonlinearAssemblerSolverHybrid_ComputeJacobian(SNES snes,
                                                                      Vec currentGuess,
                                                                      Mat* pGlobalJacobian,
                                                                      Mat* pPreconditioner,
                                                                      MatStructure* pMatStructure,
                                                                      void* pContext);








////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////
////
////         The ASSEMBLER-SOLVER class
////
////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////



/**
 * 
 *  Abstract solver for solving nonlinear elliptic PDEs. 
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractNonlinearAssemblerSolverHybrid : public AbstractFeObjectAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM,true,true,NONLINEAR>
{
protected:
    /** Mesh class */
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpMesh;

    /** Boundary conditions container */
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* mpBoundaryConditions;

    /** The nonlinear solver. */
    AbstractNonlinearSolver* mpSolver;

    /** Whether memory has been allocated for the solver. */
    bool mWeAllocatedSolverMemory;

    /** Whether to use an analytical expression for the Jacobian. */
    bool mUseAnalyticalJacobian;



    /**
     * Apply Dirichlet boundary conditions to either the residual or Jacobian.
     *
     * @param currentGuess  the solution guess for the current nonlinear solver iteration
     * @param residual the residual to apply boundary conditions to (can be NULL if bcs to 
     *   be applied to jacobian only)
     * @param pMat the jacobian to apply boundary conditions to (can be NULL if bcs to 
     *   be applied to residual only)
     */
    void ApplyDirichletConditions(Vec currentGuess, Vec residual, Mat* pMat);

    /**
     * Computes the Jacobian numerically i.e. an approximation, using numerical
     * partial derivatives.
     *
     * @param currentGuess Independent variable, u in f(u), for example
     * @param pJacobian A pointer to the Jacobian matrix
     */
    void ComputeJacobianNumerically(const Vec currentGuess, Mat* pJacobian);

public :

    /**
     * Compute the residual vector given the current solution guess.
     *
     * @param currentGuess The solution guess for the current nonlinear solve iteration.
     * @param residualVector We fill this with the residual vector.
     *
     * NOTE: this method is called indirectly by the PETSc iterative
     * solvers, so must be public.
     */
    void ComputeResidual(const Vec currentGuess, Vec residualVector);

    /**
     * Compute the Jacobian matrix given a current guess at the solution.
     * Choose whether to use a numerical or analytical method based on a flag
     * provided by the user (in Solve()).
     *
     * @param currentGuess The solution guess for the current iteration.
     * @param pJacobian Pointer to object to fill with the jacobian matrix.
     *
     * NOTE: this method is called indirectly by the PETSc iterative
     * solvers, so must be public.
     */
    void ComputeJacobian(const Vec currentGuess, Mat* pJacobian);



    /**
     * Constructor
     *
     * @param pMesh The mesh
     * @param pBoundaryConditions The boundary conditions to apply
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    AbstractNonlinearAssemblerSolverHybrid(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                           BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pBoundaryConditions,
                                           unsigned numQuadPoints = 2);

    /**
     * Destructor.
     */
    virtual ~AbstractNonlinearAssemblerSolverHybrid();

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
     *  A helper method for use when writing concrete assemblers. Once the user has calculated
     *  (on paper) the weak form and the form of the ComputeMatrixTerm method, they can check
     *  whether the analytic Jacobian matches the numerical Jacobian to verify the 
     *  correctness of the code.
     *
     *  @param tol A tolerance which defaults to 1e-5
     *  @return true if the componentwise difference between the matrices is less than
     *    the tolerance, false otherwise.
     *
     *  This method should NOT be run in simulations - it is only to verify the correctness
     *  of the concrete assembler code. Allocates dense matrices.
     */
    bool VerifyJacobian(double tol);
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AbstractNonlinearAssemblerSolverHybrid(
            AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
            BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pBoundaryConditions,
            unsigned numQuadPoints)
    :  AbstractFeObjectAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM,true,true,NONLINEAR>(pMesh,numQuadPoints),
       mpMesh(pMesh),
       mpBoundaryConditions(pBoundaryConditions)
{
    // if this is run with SPACE_DIM != ELEMENT_DIM the class has to be checked -
    // lots of places where should be using SPACE_DIM not ELEMENT_DIM??
    assert(SPACE_DIM==ELEMENT_DIM);
    assert(pMesh!=NULL);
    assert(pBoundaryConditions!=NULL);

   // mpAssembler = new SimpleNonlinearEllipticPdeAssembler<ELEMENT_DIM,SPACE_DIM>(mpMesh,pPde,numQuadPoints);
    mpSolver = new SimplePetscNonlinearSolver;
    mWeAllocatedSolverMemory = true;

    assert(mpMesh->GetNumNodes() == mpMesh->GetDistributedVectorFactory()->GetProblemSize());
    mpMesh->SetElementOwnerships(mpMesh->GetDistributedVectorFactory()->GetLow(),
                                 mpMesh->GetDistributedVectorFactory()->GetHigh());
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::~AbstractNonlinearAssemblerSolverHybrid()
{
    if (mWeAllocatedSolverMemory)
    {
        delete mpSolver;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ApplyDirichletConditions(Vec currentGuess, Vec residual, Mat* pJacobian)
{
    if (residual)
    {
        mpBoundaryConditions->ApplyDirichletToNonlinearResidual( currentGuess, 
                                                                 residual, 
                                                                 *(mpMesh->GetDistributedVectorFactory()));
    }
    if (pJacobian)
    {
        mpBoundaryConditions->ApplyDirichletToNonlinearJacobian(*pJacobian);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeResidual(const Vec currentGuess, Vec residualVector)
{
    this->SetVectorToAssemble(residualVector,true);
    this->SetCurrentSolution(currentGuess);
    this->SetApplyNeummanBoundaryConditionsToVector(mpBoundaryConditions);
    this->AssembleVector();

    VecAssemblyBegin(residualVector);
    VecAssemblyEnd(residualVector);

    ApplyDirichletConditions(currentGuess, residualVector, NULL);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeJacobian(const Vec currentGuess, Mat* pJacobian)
{
    if (mUseAnalyticalJacobian)
    {
        this->SetMatrixToAssemble(*pJacobian);
        this->SetCurrentSolution(currentGuess);
        this->AssembleMatrix();
    
        MatAssemblyBegin(*pJacobian, MAT_FLUSH_ASSEMBLY);
        MatAssemblyEnd(*pJacobian, MAT_FLUSH_ASSEMBLY);
        
        ApplyDirichletConditions(currentGuess, NULL, pJacobian);

        MatAssemblyBegin(*pJacobian, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(*pJacobian, MAT_FINAL_ASSEMBLY);
    }
    else
    {
        ComputeJacobianNumerically(currentGuess, pJacobian);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeJacobianNumerically(const Vec currentGuess, Mat* pJacobian)
{
    unsigned num_unknowns = PROBLEM_DIM * this->mpMesh->GetNumNodes();

    // Set up working vectors
    Vec residual;
    Vec perturbed_residual;
    Vec result;
    residual = PetscTools::CreateVec(num_unknowns);
    result = PetscTools::CreateVec(num_unknowns);
    perturbed_residual = PetscTools::CreateVec(num_unknowns);
    
    // Copy the currentGuess vector; we perturb the copy
    Vec current_guess_copy;
    PETSCEXCEPT( VecDuplicate(currentGuess, &current_guess_copy) );
    PETSCEXCEPT( VecCopy(currentGuess, current_guess_copy) );

    // Compute the current residual
    ComputeResidual(currentGuess, residual);

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
        ComputeResidual(current_guess_copy, perturbed_residual);

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
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::Solve(Vec initialGuess,
                                                                                       bool useAnalyticalJacobian)
{
    assert(initialGuess!=NULL);
    mUseAnalyticalJacobian = useAnalyticalJacobian;

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

    // run the solver, telling it which global functions to call in order to assemble
    // the residual or jacobian
    return mpSolver->Solve(&AbstractNonlinearAssemblerSolverHybrid_ComputeResidual<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>,
                           &AbstractNonlinearAssemblerSolverHybrid_ComputeJacobian<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>,
                           initialGuess,
                           PROBLEM_DIM * this->mpMesh->CalculateMaximumNodeConnectivityPerProcess(), 
                           this);
}




template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetNonlinearSolver(AbstractNonlinearSolver* pNonlinearSolver)
{
    if (mWeAllocatedSolverMemory)
    {
        delete mpSolver;
    }
    mpSolver = pNonlinearSolver;
    mWeAllocatedSolverMemory = false;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::VerifyJacobian(double tol)
{
    unsigned size = PROBLEM_DIM * this->mpMesh->GetNumNodes();

    Vec initial_guess = PetscTools::CreateAndSetVec(size, 0.0);
    
    Mat analytic_jacobian; 
    Mat numerical_jacobian;

    PetscTools::SetupMat(analytic_jacobian, size, size, size);
    PetscTools::SetupMat(numerical_jacobian, size, size, size);

    mUseAnalyticalJacobian = true;
    ComputeJacobian(initial_guess, &analytic_jacobian);

    mUseAnalyticalJacobian = false;
    ComputeJacobian(initial_guess, &numerical_jacobian);
    
    double minus_one = -1.0;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    // MatAYPX(*a, X, Y) does  Y = X + a*Y. 
    MatAYPX(&minus_one, analytic_jacobian, numerical_jacobian);
#elif (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 1) //PETSc 2.3.1
    // MatAYPX( Y, a, X, structure) does Y = a*Y + X. 
    MatAYPX(numerical_jacobian, minus_one, analytic_jacobian); 
#else
    // MatAYPX( Y, a, X, structure) does Y = a*Y + X. 
    MatAYPX(numerical_jacobian, minus_one, analytic_jacobian, DIFFERENT_NONZERO_PATTERN);
#endif
    double norm;
    MatNorm(numerical_jacobian,NORM_INFINITY,&norm);
    
    MatDestroy(numerical_jacobian);
    MatDestroy(analytic_jacobian);
    VecDestroy(initial_guess);

    return (norm<tol);
}




////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////
////
////   Implementations of GLOBAL functions called by Petsc nonlinear solver
////
////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PetscErrorCode AbstractNonlinearAssemblerSolverHybrid_ComputeResidual(SNES snes, Vec currentGuess,
                                                                      Vec residualVector, void* pContext)
{
    // Extract the solver from the void*
    AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* p_solver =
        (AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*) pContext;

    p_solver->ComputeResidual(currentGuess, residualVector);

    return 0;
}




template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PetscErrorCode AbstractNonlinearAssemblerSolverHybrid_ComputeJacobian(SNES snes, Vec currentGuess,
                                                                      Mat* pGlobalJacobian, Mat* pPreconditioner,
                                                                      MatStructure* pMatStructure, void* pContext)
{
    // Extract the solver from the void*
    AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* p_solver =
        (AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*) pContext;

    p_solver->ComputeJacobian(currentGuess, pGlobalJacobian);

    return 0;
}


#endif /*ABSTRACTNONLINEARASSEMBLERSOLVERHYBRID_HPP_*/
