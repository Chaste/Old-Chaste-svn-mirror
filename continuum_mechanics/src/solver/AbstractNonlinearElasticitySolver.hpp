/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef ABSTRACTNONLINEARELASTICITYSOLVER_HPP_
#define ABSTRACTNONLINEARELASTICITYSOLVER_HPP_

#include <vector>
#include <cmath>
#include "LinearSystem.hpp"
#include "OutputFileHandler.hpp"
#include "LogFile.hpp"
#include "PetscTools.hpp"
#include "MechanicsEventHandler.hpp"
#include "ReplicatableVector.hpp"
#include "FourthOrderTensor.hpp"
#include "QuadraticMesh.hpp"
#include "GaussianQuadratureRule.hpp"
#include "CmguiDeformedSolutionsWriter.hpp"
#include "AbstractMaterialLaw.hpp"
#include "Warnings.hpp"
#include "PetscException.hpp"

#include "SolidMechanicsProblemDefinition.hpp"
#include "DeformedBoundaryElement.hpp"

// Bizarrely PETSc 2.2 has this, but doesn't put it in the petscksp.h header...
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
extern PetscErrorCode KSPInitialResidual(KSP,Vec,Vec,Vec,Vec,Vec);
#endif

//#define MECH_VERBOSE      // Print output on how nonlinear solve is progressing
//#define MECH_VERY_VERBOSE // See number of elements done whilst assembling vectors or matrices
//#define MECH_USE_HYPRE    // uses HYPRE to solve linear systems, requires Petsc to be installed with HYPRE
//#define MECH_KSP_MONITOR  // Print residual norm each iteration in linear solve (ie -ksp_monitor).

//EMTODO: allow change of max_iter, better preconditioner

#ifdef MECH_VERBOSE
#include "Timer.hpp"
#endif

typedef enum CompressibilityType_
{
    COMPRESSIBLE,
    INCOMPRESSIBLE
} CompressibilityType;

/**
 * Abstract nonlinear elasticity solver.
 */
template<unsigned DIM>
class AbstractNonlinearElasticitySolver
{
protected:

    /** Number of vertices per element. */
    static const size_t NUM_VERTICES_PER_ELEMENT = DIM+1;

    /** Number of nodes per element. */
    static const size_t NUM_NODES_PER_ELEMENT = (DIM+1)*(DIM+2)/2;      // assuming quadratic

    /** Number of nodes per boundary element. */
    static const size_t NUM_NODES_PER_BOUNDARY_ELEMENT = DIM*(DIM+1)/2; // assuming quadratic

    /**
     * Maximum absolute tolerance for Newton solve. The Newton solver uses the absolute tolerance
     * corresponding to the specified relative tolerance, but has a max and min allowable absolute
     * tolerance. (ie: if max_abs = 1e-7, min_abs = 1e-10, rel=1e-4: then if the norm of the
     * initial_residual (=a) is 1e-2, it will solve with tolerance 1e-7; if a=1e-5, it will solve
     * with tolerance 1e-9; a=1e-9, it will solve with tolerance 1e-10.
     */
    static double MAX_NEWTON_ABS_TOL;

    /** Minimum absolute tolerance for Newton solve. See documentation for MAX_NEWTON_ABS_TOL. */
    static double MIN_NEWTON_ABS_TOL;

    /** Relative tolerance for Newton solve. See documentation for MAX_NEWTON_ABS_TOL. */
    static double NEWTON_REL_TOL;

    /**
     * The mesh to be solved on. Requires 6 nodes per triangle (or 10 per tetrahedron)
     * as quadratic bases are used.
     */
    QuadraticMesh<DIM>& mrQuadMesh;

    /**
     *  This class contains all the information about the problem (except the material law):
     *  body force, surface tractions, fixed nodes, density
     */
    SolidMechanicsProblemDefinition<DIM>& mrProblemDefinition;

    /** Gaussian quadrature rule. */
    GaussianQuadratureRule<DIM>* mpQuadratureRule;

    /** Boundary Gaussian quadrature rule. */
    GaussianQuadratureRule<DIM-1>* mpBoundaryQuadratureRule;

    /**
     * Absolute tolerance for linear systems. Can be set by calling
     * SetKspAbsoluteTolerances(), but default to -1, in which case
     * a relative tolerance is used.
     */
    double mKspAbsoluteTol;

    /**
     * Number of degrees of freedom (equal to, in the incompressible case:
     * DIM*N + M if quadratic-linear bases are used, where there are N total
     * nodes and M vertices; or DIM*N in the compressible case).
     */
    unsigned mNumDofs;

    /**
     * Residual vector for the full nonlinear system, also the RHS vector in the linear
     * system used to solve the nonlinear problem using Newton's method.
     */
    Vec mResidualVector;

    /**
     * The RHS side in the linear system that is solved each Newton iteration. Since Newton's method
     * is Ju = f, where J is the Jacobian, u the (negative of the) update and f the residual, it might seem necessary
     * to store this as well as the residual. However, when applying Dirichlet boundary conditions in
     * the compressible case, we alter the rows of the matrix, but also alter the columns in order to
     * maintain symmetry. This requires making further changes to the right-hand vector, meaning that
     * it no longer properly represents the residual. Hence, we have to use two vectors.
     *
     * Overall, this can be represents as
     *  - compute residual f
     *  - compute Jacobian J
     *  - apply BCs to f.
     *  - alter the linear system from Ju=f to (J*)u=f* which enforces the dirichlet boundary conditions but enforces them symmetrically.
     *
     * mLinearSystemRhsVector represents f*.
     */
    Vec mLinearSystemRhsVector;

    /**
     * Jacobian matrix of the nonlinear system, LHS matrix for the linear system.
     */
    Mat mJacobianMatrix;

    /**
     * Helper vector (see ApplyDirichletBoundaryConditions code).
     */
    Vec mDirichletBoundaryConditionsVector;

    /**
     * Precondition matrix for the linear system.
     *
     * In the incompressible case:
     * the preconditioner is the petsc LU factorisation of
     *
     * Jp = [A B] in displacement-pressure block form,
     *      [C M]
     *
     * where the A, B and C are the matrices in the normal jacobian,
     * i.e.
     *
     * J  = [A B]
     *      [C 0]
     *
     * and M is the MASS MATRIX (ie integral phi_i phi_j dV, where phi_i are the
     * pressure basis functions).
     */
    Mat mPreconditionMatrix;

    /** Whether to write any output. */
    bool mWriteOutput;

    /** Where to write output, relative to CHASTE_TESTOUTPUT. */
    std::string mOutputDirectory;

    /** Output file handler. */
    OutputFileHandler* mpOutputFileHandler;

    /**
     * By default only the initial and final solutions are written. However, we may
     * want to write the solutions after every Newton iteration, in which case the
     * following should be set to true.
     */
    bool mWriteOutputEachNewtonIteration;

    /**
     * The current solution, in the form (assuming 2d):
     *   Incompressible problem: [u1 v1 u2 v2 ... uN vN p1 p2 .. pM]
     *   Compressible problem:   [u1 v1 u2 v2 ... uN vN]
     * where there are N total nodes and M vertices.
     */
    std::vector<double> mCurrentSolution;

    /**
     * Storage space for a 4th order tensor used in assembling the Jacobian (to avoid repeated memory allocation).
     */
    FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE;

    /** Number of Newton iterations taken in last solve. */
    unsigned mNumNewtonIterations;

    /** Deformed position: mDeformedPosition[i](j) = x_j for node i. */
    std::vector<c_vector<double,DIM> > mDeformedPosition;

    /**
     *  For a particular type of boundary condition - prescribed pressures acting on the deformed
     *  surface, we need to do surface integrals over the deformed surface. This object is a helper
     *  object for this, it takes in a base boundary element and a set of displacements, and
     *  sets up the deformed boundary element. Only used if
     *  mProblemDefinition.GetTractionBoundaryConditionType()==PRESSURE_ON_DEFORMED)
     */
    DeformedBoundaryElement<DIM-1,DIM> mDeformedBoundaryElement;

    /**
     * This solver is for static problems, however the body force or surface tractions
     * could be a function of time. The user should call SetCurrentTime() if this is
     * the case.
     */
    double mCurrentTime;

    /**
     * This is equal to either COMPRESSIBLE or INCOMPRESSIBLE (see enumeration defined at top of file)
     * and is only used in computing mNumDofs and allocating matrix memory.
     */
    CompressibilityType mCompressibilityType;

    /** Whether the boundary elements of the mesh have been checked for whether
     *  the ordering if such that the normals are outward-facing
     */
    bool mCheckedOutwardNormals;

    /**
     * Assemble the residual vector and/or Jacobian matrix (using the current solution stored
     * in mCurrentSolution, output going to mResidualVector and/or mJacobianMatrix).
     *
     * Must be overridden in concrete derived classes.
     *
     * @param assembleResidual     A bool stating whether to assemble the residual vector.
     * @param assembleLinearSystem A bool stating whether to assemble the Jacobian matrix and the RHS
     *  vector of the linear system (which is based on the residual but could be slightly different
     *  due to the way dirichlet boundary conditions are applied to the linear system - see comments in
     *  ApplyDirichletBoundaryConditions).
     */
    virtual void AssembleSystem(bool assembleResidual, bool assembleLinearSystem)=0;

    /**
     * Initialise the solver.
     */
    void Initialise();

    /**
     * Allocates memory for the Jacobian and preconditioner matrices (larger number of
     * non-zeros per row than with say linear problems)
     */
    void AllocateMatrixMemory();

    /**
     * Apply the Dirichlet boundary conditions to the linear system.
     *
     * This will always apply the Dirichlet boundary conditions to the residual vector
     * (basically, setting a component to the difference between the current value and
     * the correct value).
     *
     * If the boolean parameter is true, this will apply the boundary conditions to the
     * Jacobian and the linear system RHS vector (which should be equal to the residual
     * on entering this function). In the compressible case the boundary conditions are
     * applied by zeroing both rows and columns of the Jacobian matrix (to maintain)
     * symmetry, which means additional changes are needed for the RHS vector.
     *
     * @param applyToMatrix Whether to apply the boundary conditions to the linear system
     *     (as well as the residual).
     */
    void ApplyDirichletBoundaryConditions(bool applyToMatrix);

    /**
     * To be called at the end of AssembleSystem. Calls (Petsc) assemble methods on the
     * Vecs and Mat, and calls ApplyDirichletBoundaryConditions.
     *
     * @param assembleResidual see documentation for AssembleSystem
     * @param assembleLinearSystem see documentation for AssembleSystem
     */
    virtual void FinishAssembleSystem(bool assembleResidual, bool assembleLinearSystem);

    /**
     * Set up the residual vector (using the current solution), and get its
     * scaled norm (Calculate |r|_2 / length(r), where r is residual vector).
     *
     * @param allowException Sometimes the current solution solution will be such
     *   that the residual vector cannot be computed, as (say) the material law
     *   will throw an exception as the strains are too large. If this bool is
     *   set to true, the exception will be caught, and DBL_MAX returned as the
     *   residual norm
     */
    double ComputeResidualAndGetNorm(bool allowException);

    /** Calculate |r|_2 / length(r), where r is the current residual vector. */
    double CalculateResidualNorm();

    /**
     * Simple helper function, computes Z = X + aY, where X and Z are std::vectors and
     * Y a ReplicatableVector.
     *
     * @param rX X
     * @param rY Y (replicatable vector)
     * @param a a
     * @param rZ Z the returned vector
     */
    void VectorSum(std::vector<double>& rX, ReplicatableVector& rY, double a, std::vector<double>& rZ);

    /**
     * Print to std::cout the residual norm for this s, ie ||f(x+su)|| where f is the residual vector,
     * x the current solution and u the update vector.
     *
     * @param s s
     * @param residNorm residual norm.
     */
    void PrintLineSearchResult(double s, double residNorm);

    /**
     * Take one newton step, by solving the linear system -Ju=f, (J the jacobian, f
     * the residual, u the update), and picking s such that a_new = a_old + su (a
     * the current solution) such |f(a)| is the smallest.
     *
     * @return The current norm of the residual after the newton step.
     */
    double TakeNewtonStep();

    /**
     * Using the update vector (of Newton's method), choose s such that ||f(x+su)|| is most decreased,
     * where f is the residual vector, x the current solution (mCurrentSolution) and u the update vector.
     * This checks s=1 first (most likely to be the current solution, then 0.9, 0.8.. until ||f|| starts
     * increasing.
     *
     * @param solution The solution of the linear solve in newton's method, ie the update vector u.
     */
    double UpdateSolutionUsingLineSearch(Vec solution);

    /**
     * This function may be overloaded by subclasses. It is called after each Newton
     * iteration.
     *
     * @param counter current newton iteration number
     * @param normResidual norm of the residual
     */
    virtual void PostNewtonStep(unsigned counter, double normResidual);

    /**
     * Simple (one-line function which just calls ComputeStressAndStressDerivative() on the
     * material law given, using C,  inv(C), and p as the input and with rT and rDTdE as the
     * output. Overloaded by other assemblers (eg cardiac mechanics) which need to add extra
     * terms to the stress.
     *
     * @param pMaterialLaw The material law for this element
     * @param rC The Lagrangian deformation tensor (F^T F)
     * @param rInvC The inverse of C. Should be computed by the user.
     * @param pressure The current pressure
     * @param elementIndex Index of the current element
     * @param currentQuadPointGlobalIndex The index (assuming an outer loop over elements and an inner
     *   loop over quadrature points), of the current quadrature point.
     * @param rT The stress will be returned in this parameter
     * @param rDTdE the stress derivative will be returned in this parameter, assuming
     *   the final parameter is true
     * @param computeDTdE A boolean flag saying whether the stress derivative is
     *   required or not.
     */
    virtual void ComputeStressAndStressDerivative(AbstractMaterialLaw<DIM>* pMaterialLaw,
                                                  c_matrix<double,DIM,DIM>& rC,
                                                  c_matrix<double,DIM,DIM>& rInvC,
                                                  double pressure,
                                                  unsigned elementIndex,
                                                  unsigned currentQuadPointGlobalIndex,
                                                  c_matrix<double,DIM,DIM>& rT,
                                                  FourthOrderTensor<DIM,DIM,DIM,DIM>& rDTdE,
                                                  bool computeDTdE)
    {
        // Just call the method on the material law
        pMaterialLaw->ComputeStressAndStressDerivative(rC,rInvC,pressure,rT,rDTdE,computeDTdE);
    }

public:

    /**
     * Constructor.
     *
     * @param rQuadMesh  the quadratic mesh
     * @param rProblemDefinition an object defining in particular the body force and boundary conditions
     * @param outputDirectory output directory
     * @param compressibilityType Should be equal to COMPRESSIBLE or INCOMPRESSIBLE (see enumeration defined at top of file)
     *   (depending on which concrete class is inheriting from this) and is only used in computing mNumDofs and allocating
     *   matrix memory.
     */
    AbstractNonlinearElasticitySolver(QuadraticMesh<DIM>& rQuadMesh,
                                      SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                      std::string outputDirectory,
                                      CompressibilityType compressibilityType);

    /**
     * Destructor.
     */
    virtual ~AbstractNonlinearElasticitySolver();

    /**
     * Solve the problem.
     *
     * @param tol tolerance used in Newton solve (defaults to -1.0)
     * @param maxNumNewtonIterations (defaults to INT_MAX)
     * @param quitIfNoConvergence (defaults to true)
     */
    void Solve(double tol=-1.0,
               unsigned maxNumNewtonIterations=INT_MAX,
               bool quitIfNoConvergence=true);

    /**
     * Write the current deformation of the nodes.
     *
     * @param fileName (stem)
     * @param counterToAppend append a counter to the file name
     *
     * For example:
     * WriteCurrentDeformation("solution") --> file called "solution.nodes"
     * WriteCurrentDeformation("solution",3) --> file called "solution_3.nodes"
     */
    void WriteCurrentDeformation(std::string fileName, int counterToAppend=-1);

    /**
     * Get number of Newton iterations taken in last solve.
     */
    unsigned GetNumNewtonIterations();


    /**
     * Set whether to write any output.
     *
     * @param writeOutput (defaults to true)
     */
    void SetWriteOutput(bool writeOutput=true);

    /**
     * By default only the original and converged solutions are written. Call this
     * by get node positions written after every Newton step (mostly for
     * debugging).
     *
     * @param writeOutputEachNewtonIteration whether to write each iteration
     */
    void SetWriteOutputEachNewtonIteration(bool writeOutputEachNewtonIteration=true)
    {
        mWriteOutputEachNewtonIteration = writeOutputEachNewtonIteration;
    }

    /**
     * Set the absolute tolerance to be used when solving the linear system.
     * If this is not called a relative tolerance is used.
     *
     * @param kspAbsoluteTolerance the tolerance
     */
    void SetKspAbsoluteTolerance(double kspAbsoluteTolerance)
    {
        assert(kspAbsoluteTolerance > 0);
        mKspAbsoluteTol = kspAbsoluteTolerance;
    }

    /**
     * Get the current solution vector (advanced use only - for getting the deformed position use
     * rGetDeformedPosition().
     */
    std::vector<double>& rGetCurrentSolution()
    {
        return mCurrentSolution;
    }


    /**
     * Get the deformed position. Note returnvalue[i](j) = x_j for node i.
     */
    std::vector<c_vector<double,DIM> >& rGetDeformedPosition();

    /**
     * This solver is for static problems, however the body force or surface tractions
     * could be a function of time. This method is for setting the time.
     *
     * @param time current time
     */
    void SetCurrentTime(double time)
    {
        mCurrentTime = time;
    }

    /**
     * Convert the output to Cmgui format (placed in a folder called cmgui in the output directory).
     * Writes the original mesh as solution_0.exnode and the (current) solution as solution_1.exnode.
     */
    void CreateCmguiOutput();
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::Initialise()
{
    AllocateMatrixMemory();

    mpQuadratureRule = new GaussianQuadratureRule<DIM>(3);
    mpBoundaryQuadratureRule = new GaussianQuadratureRule<DIM-1>(3);

    mCurrentSolution.resize(mNumDofs, 0.0);
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::AllocateMatrixMemory()
{
    mResidualVector = PetscTools::CreateVec(mNumDofs);
    VecDuplicate(mResidualVector, &mLinearSystemRhsVector);
    if (mCompressibilityType == COMPRESSIBLE)
    {
        VecDuplicate(mResidualVector, &mDirichletBoundaryConditionsVector);
    }

    if (DIM==2)
    {
        // 2D: N elements around a point => 7N+3 non-zeros in that row? Assume N<=10 (structured mesh would have N_max=6) => 73.
        unsigned num_non_zeros = std::min(75u, mNumDofs);

        PetscTools::SetupMat(mJacobianMatrix, mNumDofs, mNumDofs, num_non_zeros, PETSC_DECIDE, PETSC_DECIDE);
        PetscTools::SetupMat(mPreconditionMatrix, mNumDofs, mNumDofs, num_non_zeros, PETSC_DECIDE, PETSC_DECIDE);
    }
    else
    {
        assert(DIM==3);

        // 3D: N elements around a point. nz < (3*10+6)N (lazy estimate). Better estimate is 23N+4?. Assume N<20 => 500ish

        // in 3d we get the number of containing elements for each node and use that to obtain an upper bound
        // for the number of non-zeros for each DOF associated with that node.

        int* num_non_zeros_each_row = new int[mNumDofs];
        for (unsigned i=0; i<mNumDofs; i++)
        {
            num_non_zeros_each_row[i] = 0;
        }

        for (unsigned i=0; i<mrQuadMesh.GetNumNodes(); i++)
        {
            // this upper bound neglects the fact that two containing elements will share the same nodes..
            // 4 = max num dofs associated with this node
            // 30 = 3*9+3 = 3 dimensions x 9 other nodes on this element   +  3 vertices with a pressure unknown
            unsigned num_non_zeros_upper_bound = 4 + 30*mrQuadMesh.GetNode(i)->GetNumContainingElements();

            num_non_zeros_upper_bound = std::min(num_non_zeros_upper_bound, mNumDofs);

            num_non_zeros_each_row[DIM*i + 0] = num_non_zeros_upper_bound;
            num_non_zeros_each_row[DIM*i + 1] = num_non_zeros_upper_bound;
            num_non_zeros_each_row[DIM*i + 2] = num_non_zeros_upper_bound;

            if (mCompressibilityType==INCOMPRESSIBLE)
            {
                //Could do !mrQuadMesh.GetNode(i)->IsInternal()
                if (i<mrQuadMesh.GetNumVertices()) // then this is a vertex
                {
                    num_non_zeros_each_row[DIM*mrQuadMesh.GetNumNodes() + i] = num_non_zeros_upper_bound;
                }
            }
        }

        // NOTE: PetscTools::SetupMat() or the below creates a MATAIJ matrix, which means the matrix will
        // be of type MATSEQAIJ if num_procs=1 and MATMPIAIJ otherwise. In the former case
        // MatSeqAIJSetPreallocation MUST be called [MatMPIAIJSetPreallocation will have
        // no effect (silently)], and vice versa in the latter case

        /// We want to allocate different numbers of non-zeros per row, which means
        /// PetscTools::SetupMat isn't that useful. We could call
        //PetscTools::SetupMat(mJacobianMatrix, mNumDofs, mNumDofs, 0, PETSC_DECIDE, PETSC_DECIDE);
        //PetscTools::SetupMat(mPreconditionMatrix, mNumDofs, mNumDofs, 0, PETSC_DECIDE, PETSC_DECIDE);
        /// but we would get warnings due to the lack allocation

        // possible todo: create a PetscTools::SetupMatNoAllocation()

        // In the future, when parallelising, remember to think about MAT_IGNORE_OFF_PROC_ENTRIES (see #1682)

#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
        MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,mNumDofs,mNumDofs,&mJacobianMatrix);
        MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,mNumDofs,mNumDofs,&mPreconditionMatrix);
#else //New API
        MatCreate(PETSC_COMM_WORLD,&mJacobianMatrix);
        MatCreate(PETSC_COMM_WORLD,&mPreconditionMatrix);
        MatSetSizes(mJacobianMatrix,PETSC_DECIDE,PETSC_DECIDE,mNumDofs,mNumDofs);
        MatSetSizes(mPreconditionMatrix,PETSC_DECIDE,PETSC_DECIDE,mNumDofs,mNumDofs);
#endif

        if (PetscTools::IsSequential())
        {
            MatSetType(mJacobianMatrix, MATSEQAIJ);
            MatSetType(mPreconditionMatrix, MATSEQAIJ);
            MatSeqAIJSetPreallocation(mJacobianMatrix,     PETSC_NULL, num_non_zeros_each_row);
            MatSeqAIJSetPreallocation(mPreconditionMatrix, PETSC_NULL, num_non_zeros_each_row);
        }
        else
        {
            PetscInt lo, hi;
            VecGetOwnershipRange(mResidualVector, &lo, &hi);

            int* num_non_zeros_each_row_this_proc = new int[hi-lo];
            int* zero = new int[hi-lo];
            for (unsigned i=0; i<unsigned(hi-lo); i++)
            {
                num_non_zeros_each_row_this_proc[i] = num_non_zeros_each_row[lo+i];
                zero[i] = 0;
            }

            MatSetType(mJacobianMatrix, MATMPIAIJ);
            MatSetType(mPreconditionMatrix, MATMPIAIJ);
            MatMPIAIJSetPreallocation(mJacobianMatrix,     PETSC_NULL, num_non_zeros_each_row_this_proc, PETSC_NULL, num_non_zeros_each_row_this_proc);
            MatMPIAIJSetPreallocation(mPreconditionMatrix, PETSC_NULL, num_non_zeros_each_row_this_proc, PETSC_NULL, num_non_zeros_each_row_this_proc);
        }

        MatSetFromOptions(mJacobianMatrix);
        MatSetFromOptions(mPreconditionMatrix);

        //unsigned total_non_zeros = 0;
        //for (unsigned i=0; i<mNumDofs; i++)
        //{
        //   total_non_zeros += num_non_zeros_each_row[i];
        //}
        //std::cout << total_non_zeros << " versus " << 500*mNumDofs << "\n" << std::flush;

        delete [] num_non_zeros_each_row;
    }
}

/*
 * This method applies the appropriate BCs to the residual vector, and possible also the linear system.
 *
 * For the latter, in the compressible case, the BCs are imposed in such a way as to ensure that a
 * symmetric linear system remains symmetric. For each row with boundary condition applied, both the
 * row and column are zero'd and the RHS vector modified to take into account the zero'd column.
 *
 * Suppose we have a matrix
 * [a b c] [x] = [ b1 ]
 * [d e f] [y]   [ b2 ]
 * [g h i] [z]   [ b3 ]
 * and we want to apply the boundary condition x=v without losing symmetry if the matrix is
 * symmetric. We apply the boundary condition
 * [1 0 0] [x] = [ v  ]
 * [d e f] [y]   [ b2 ]
 * [g h i] [z]   [ b3 ]
 * and then zero the column as well, adding a term to the RHS to take account for the
 * zero-matrix components
 * [1 0 0] [x] = [ v  ] - v[ 0 ]
 * [0 e f] [y]   [ b2 ]    [ d ]
 * [0 h i] [z]   [ b3 ]    [ g ]
 * Note the last term is the first column of the matrix, with one component zeroed, and
 * multiplied by the boundary condition value. This last term is then stored in
 * rLinearSystem.rGetDirichletBoundaryConditionsVector(), and in general form is the
 * SUM_{d=1..D} v_d a'_d
 * where v_d is the boundary value of boundary condition d (d an index into the matrix),
 * and a'_d is the dth-column of the matrix but with the d-th component zeroed, and where
 * there are D boundary conditions
 */
template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::ApplyDirichletBoundaryConditions(bool applyToLinearSystem)
{
    assert(mResidualVector); // BCs will be added to all the time
    if (applyToLinearSystem)
    {
        assert(mJacobianMatrix);
        assert(mLinearSystemRhsVector);
    }

    // The boundary conditions on the NONLINEAR SYSTEM are x=boundary_values
    // on the boundary nodes. However:
    // The boundary conditions on the LINEAR SYSTEM  Ju=f, where J is the
    // u the negative update vector and f is the residual is
    // u=current_soln-boundary_values on the boundary nodes

    std::vector<unsigned> rows;
    std::vector<double> values;

    rows.resize(DIM*mrProblemDefinition.rGetFixedNodes().size());
    values.resize(DIM*mrProblemDefinition.rGetFixedNodes().size());

    // Whether to apply symmetrically, ie alter columns as well as rows (see comment above)
    bool applySymmetrically = (applyToLinearSystem) && (mCompressibilityType==COMPRESSIBLE);

    if (applySymmetrically)
    {
        assert(applyToLinearSystem);
        PetscVecTools::Zero(mDirichletBoundaryConditionsVector);
        PetscMatTools::Finalise(mJacobianMatrix);
    }

    for (unsigned i=0; i<mrProblemDefinition.rGetFixedNodes().size(); i++)
    {
        unsigned node_index = mrProblemDefinition.rGetFixedNodes()[i];

        for (unsigned j=0; j<DIM; j++)
        {
            double disp = mrProblemDefinition.rGetFixedNodeDisplacements()[i](j);

            if(disp != SolidMechanicsProblemDefinition<DIM>::FREE)
            {
                unsigned dof_index = DIM*node_index+j;
                rows[DIM*i+j] = dof_index;
                values[DIM*i+j] = mCurrentSolution[dof_index] - disp;
            }
        }
    }

    if (applySymmetrically)
    {
        // Modify the matrix columns
        for (unsigned i=0; i<rows.size(); i++)
        {
            unsigned col = rows[i];
            double minus_value = -values[i];

            // Get a vector which will store the column of the matrix (column d, where d is
            // the index of the row (and column) to be altered for the boundary condition.
            // Since the matrix is symmetric when get row number "col" and treat it as a column.
            // PETSc uses compressed row format and therefore getting rows is far more efficient
            // than getting columns.
            Vec matrix_col = PetscMatTools::GetMatrixRowDistributed(mJacobianMatrix,col);

            // Zero the correct entry of the column
            PetscVecTools::SetElement(matrix_col, col, 0.0);

            // Set up the RHS Dirichlet boundary conditions vector
            // Assuming one boundary at the zeroth node (x_0 = value), this is equal to
            //   -value*[0 a_21 a_31 .. a_N1]
            // and will be added to the RHS.
            PetscVecTools::AddScaledVector(mDirichletBoundaryConditionsVector, matrix_col, minus_value);
            VecDestroy(matrix_col);
        }
    }

    if (applyToLinearSystem)
    {
        // Now zero the appropriate rows and columns of the matrix. If the matrix is symmetric we apply the
        // boundary conditions in a way the symmetry isn't lost (rows and columns). If not only the row is
        // zeroed.
        if (applySymmetrically)
        {
            PetscMatTools::ZeroRowsAndColumnsWithValueOnDiagonal(mJacobianMatrix, rows, 1.0);
            PetscMatTools::ZeroRowsAndColumnsWithValueOnDiagonal(mPreconditionMatrix, rows, 1.0);

            // Apply the RHS boundary conditions modification if required.
            PetscVecTools::AddScaledVector(mLinearSystemRhsVector, mDirichletBoundaryConditionsVector, 1.0);
        }
        else
        {
            PetscMatTools::ZeroRowsWithValueOnDiagonal(mJacobianMatrix, rows, 1.0);
            PetscMatTools::ZeroRowsWithValueOnDiagonal(mPreconditionMatrix, rows, 1.0);
        }
    }

    for (unsigned i=0; i<rows.size(); i++)
    {
        PetscVecTools::SetElement(mResidualVector, rows[i], values[i]);
    }

    if (applyToLinearSystem)
    {
        for (unsigned i=0; i<rows.size(); i++)
        {
            PetscVecTools::SetElement(mLinearSystemRhsVector, rows[i], values[i]);
        }
    }
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::FinishAssembleSystem(bool assembleResidual, bool assembleJacobian)
{
    if (assembleResidual)
    {
        PetscVecTools::Finalise(mResidualVector);
    }
    if (assembleJacobian)
    {
        PetscMatTools::SwitchWriteMode(mJacobianMatrix);
        PetscMatTools::SwitchWriteMode(mPreconditionMatrix);

        VecCopy(mResidualVector, mLinearSystemRhsVector);
    }

    // Apply Dirichlet boundary conditions
    ApplyDirichletBoundaryConditions(assembleJacobian);

    if (assembleResidual)
    {
        PetscVecTools::Finalise(mResidualVector);
    }
    if (assembleJacobian)
    {
        PetscMatTools::Finalise(mJacobianMatrix);
        PetscMatTools::Finalise(mPreconditionMatrix);
        PetscVecTools::Finalise(mLinearSystemRhsVector);
    }
}

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::ComputeResidualAndGetNorm(bool allowException)
{
    if (!allowException)
    {
        // Assemble the residual
        AssembleSystem(true, false);
    }
    else
    {
        try
        {
            // Try to assemble the residual using this current solution
            AssembleSystem(true, false);
        }
        catch(Exception& e)
        {
            // If fail (because e.g. ODEs fail to solve, or strains are too large for material law), return infinity
            return DBL_MAX;
        }
    }

    // Return the scaled norm of the residual
    return CalculateResidualNorm();
}

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::CalculateResidualNorm()
{
    double norm;
    VecNorm(mResidualVector, NORM_2, &norm);
    return norm/mNumDofs;
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::VectorSum(std::vector<double>& rX,
                                                       ReplicatableVector& rY,
                                                       double a,
                                                       std::vector<double>& rZ)
{
    assert(rX.size()==rY.GetSize());
    assert(rY.GetSize()==rZ.size());
    for (unsigned i=0; i<rX.size(); i++)
    {
        rZ[i] = rX[i] + a*rY[i];
    }
}

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::TakeNewtonStep()
{
    #ifdef MECH_VERBOSE
    Timer::Reset();
    #endif

    /////////////////////////////////////////////////////////////
    // Assemble Jacobian (and preconditioner)
    /////////////////////////////////////////////////////////////
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::ASSEMBLE);
    AssembleSystem(true, true);
    MechanicsEventHandler::EndEvent(MechanicsEventHandler::ASSEMBLE);
    #ifdef MECH_VERBOSE
    Timer::PrintAndReset("AssembleSystem");
    #endif

    ///////////////////////////////////////////////////////////////////
    //
    // Solve the linear system.
    // Three alternatives
    //   (a) Incompressible: GMRES with ILU preconditioner (or bjacobi=ILU on each process) [default]. Very poor on large problems.
    //   (b) Incompressible: GMRES with AMG preconditioner. Uncomment #define MECH_USE_HYPRE above. Requires Petsc3 with HYPRE installed.
    //   (c) Compressible: CG with ???
    ///////////////////////////////////////////////////////////////////
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::SOLVE);

    Vec solution;
    VecDuplicate(mResidualVector,&solution);

    KSP solver;
    KSPCreate(PETSC_COMM_WORLD,&solver);
    PC pc;
    KSPGetPC(solver, &pc);

    KSPSetOperators(solver, mJacobianMatrix, mPreconditionMatrix, DIFFERENT_NONZERO_PATTERN /*in precond between successive solves*/);

    if (mCompressibilityType==COMPRESSIBLE)
    {
        KSPSetType(solver,KSPCG);
        //PCSetType(pc, PCSOR);
        PCSetType(pc, PCICC); /// \todo: #1913 this only works in sequential
        PetscOptionsSetValue("-pc_factor_shift_positive_definite", "");
        
		//needed for ICC preconditioner
		//PetscOptionsSetValue("-pc_factor_shift_positive_definite", "");

		//// for debugging only
		//assert( PetscMatTools::CheckSymmetry(mJacobianMatrix) );
    }
    else
    {
        unsigned num_restarts = 100;
        KSPSetType(solver,KSPGMRES);
        KSPGMRESSetRestart(solver,num_restarts);

        #ifndef MECH_USE_HYPRE
            PCSetType(pc, PCBJACOBI); // BJACOBI = ILU on each block (block = part of matrix on each process)
        #else
            /////////////////////////////////////////////////////////////////////////////////////////////////////
            // Speed up linear solve time massively for larger simulations (in fact GMRES may stagnate without
            // this for larger problems), by using a AMG preconditioner -- needs HYPRE installed
            /////////////////////////////////////////////////////////////////////////////////////////////////////
            PetscOptionsSetValue("-pc_hypre_type", "boomeramg");
            // PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter", "1");
            // PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", "0.0");

            PCSetType(pc, PCHYPRE);
            KSPSetPreconditionerSide(solver, PC_RIGHT);

            // other possible preconditioners..
            //PCBlockDiagonalMechanics* p_custom_pc = new PCBlockDiagonalMechanics(solver, mPreconditionMatrix, mBlock1Size, mBlock2Size);
            //PCLDUFactorisationMechanics* p_custom_pc = new PCLDUFactorisationMechanics(solver, mPreconditionMatrix, mBlock1Size, mBlock2Size);
            //remember to delete memory..
            //KSPSetPreconditionerSide(solver, PC_RIGHT);
        #endif
    }




    #ifdef MECH_KSP_MONITOR
    PetscOptionsSetValue("-ksp_monitor","");
	//PetscOptionsSetValue("-ksp_norm_type","natural");
    #endif

    KSPSetFromOptions(solver);
    KSPSetUp(solver);


    ///// For printing matrix when debugging
    //std::cout << "\n\n";
    //for(unsigned i=0; i<mNumDofs; i++)
    //{
    //    for(unsigned j=0; j<mNumDofs; j++)
    //    {
    //        std::cout << PetscMatTools::GetElement(mJacobianMatrix, i, j) << " ";
    //    }
    //    std::cout << "\n";
    //}
    //std::cout << "\n\n";

	// Set the linear system absolute tolerance.
	// This is either the user provided value, or set to
	// max {1e-6 * initial_residual, 1e-12}
    if (mKspAbsoluteTol < 0)
    {
        Vec temp = PetscTools::CreateVec(mNumDofs);
        Vec temp2 = PetscTools::CreateVec(mNumDofs);
        Vec linsys_residual = PetscTools::CreateVec(mNumDofs);

        KSPInitialResidual(solver, solution, temp, temp2, linsys_residual, mLinearSystemRhsVector);
        double initial_resid_norm;
        VecNorm(linsys_residual, NORM_2, &initial_resid_norm);

        VecDestroy(temp);
        VecDestroy(temp2);
        VecDestroy(linsys_residual);
    
        double ksp_rel_tol = 1e-6;
        double absolute_tol = ksp_rel_tol * initial_resid_norm;
        if(absolute_tol < 1e-12)
        {
            absolute_tol = 1e-12;
        }
        KSPSetTolerances(solver, 1e-16, absolute_tol, PETSC_DEFAULT, PETSC_DEFAULT /* max iters */); // Note, max iters seems to be 1000 whatever we give here
    }
    else
    {
        KSPSetTolerances(solver, 1e-16, mKspAbsoluteTol, PETSC_DEFAULT, PETSC_DEFAULT /* max iters */); // Note, max iters seems to be 1000 whatever we give here
    }

    #ifdef MECH_VERBOSE
    Timer::PrintAndReset("KSP Setup");
    #endif

    KSPSolve(solver,mLinearSystemRhsVector,solution);

    /////////////////////////////////////////////
    // Error checking for linear solve
    /////////////////////////////////////////////

    // warn if ksp reports failure
    KSPConvergedReason reason;
    KSPGetConvergedReason(solver,&reason);

    if(reason != KSP_DIVERGED_ITS)
    {
        // Throw an exception if the solver failed for any reason other than DIVERGED_ITS.
        // This is not covered as would be difficult to cover - requires a bad matrix to
        // assembled, for example.
        #define COVERAGE_IGNORE
        KSPEXCEPT(reason);
        #undef COVERAGE_IGNORE
    }
    else
    {
        // DIVERGED_ITS just means it didn't converge in the given maximum number of iterations,
        // which is potentially not a problem, as the nonlinear solver may (and often will) still converge.
        // Just warn once.
        // (Very difficult to cover in normal tests, requires relative and absolute ksp tols to be very small, there
        // is no interface for setting both of these. Could be covered by setting up a problem the solver
        // finds difficult to solve, but this would be overkill.)
        #define COVERAGE_IGNORE
        WARN_ONCE_ONLY("Linear solve (with a mechanics solve) didn't converge, but this may not stop nonlinear solve converging")
        #undef COVERAGE_IGNORE
    }

    // quit if no ksp iterations were done
    int num_iters;
    KSPGetIterationNumber(solver, &num_iters);
    if (num_iters==0)
    {
        VecDestroy(solution);
        KSPDestroy(solver);
        EXCEPTION("KSP Absolute tolerance was too high, linear system wasn't solved - there will be no decrease in Newton residual. Decrease KspAbsoluteTolerance");
    }

    // warn if max ksp iterations was done
    if (num_iters==1000) // See comment on max_iters above
    {
        #define COVERAGE_IGNORE
        WARNING("Linear solver in mechanics solve may not have converged");
        #undef COVERAGE_IGNORE
    }

    #ifdef MECH_VERBOSE
    Timer::PrintAndReset("KSP Solve");
    std::cout << "[" << PetscTools::GetMyRank() << "]: Num iterations = " << num_iters << "\n" << std::flush;
    #endif

    MechanicsEventHandler::EndEvent(MechanicsEventHandler::SOLVE);

    ///////////////////////////////////////////////////////////////////////////
    // Update the solution
    //  Newton method:       sol = sol - update, where update=Jac^{-1}*residual
    //  Newton with damping: sol = sol - s*update, where s is chosen
    //   such that |residual(sol)| is minimised. Damping is important to
    //   avoid initial divergence.
    //
    // Normally, finding the best s from say 0.05,0.1,0.2,..,1.0 is cheap,
    // but this is not the case in cardiac electromechanics calculations.
    // Therefore, we initially check s=1 (expected to be the best most of the
    // time, then s=0.9. If the norm of the residual increases, we assume
    // s=1 is the best. Otherwise, check s=0.8 to see if s=0.9 is a local min.
    ///////////////////////////////////////////////////////////////////////////
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::UPDATE);
    double new_norm_resid = UpdateSolutionUsingLineSearch(solution);
    MechanicsEventHandler::EndEvent(MechanicsEventHandler::UPDATE);

    VecDestroy(solution);
    KSPDestroy(solver);

    return new_norm_resid;
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::PrintLineSearchResult(double s, double residNorm)
{
    #ifdef MECH_VERBOSE
    std::cout << "\tTesting s = " << s << ", |f| = " << residNorm << "\n" << std::flush;
    #endif
}

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::UpdateSolutionUsingLineSearch(Vec solution)
{
    double initial_norm_resid = CalculateResidualNorm();
    #ifdef MECH_VERBOSE
    std::cout << "\tInitial |f| [corresponding to s=0] is " << initial_norm_resid << "\n"  << std::flush;
    #endif

    ReplicatableVector update(solution);

    std::vector<double> old_solution = mCurrentSolution;

    std::vector<double> damping_values; // = {1.0, 0.9, .., 0.2, 0.1, 0.05} ie size 11
    for (unsigned i=10; i>=1; i--)
    {
        damping_values.push_back((double)i/10.0);
    }
    damping_values.push_back(0.05);
    assert(damping_values.size()==11);

    //// Try s=1 and see what the residual-norm is
    // let mCurrentSolution = old_solution - damping_val[0]*update; and compute residual
    unsigned index = 0;
    VectorSum(old_solution, update, -damping_values[index], mCurrentSolution);
    double current_resid_norm = ComputeResidualAndGetNorm(true);
    PrintLineSearchResult(damping_values[index], current_resid_norm);

    //// Try s=0.9 and see what the residual-norm is
    // let mCurrentSolution = old_solution - damping_val[1]*update; and compute residual
    index = 1;
    VectorSum(old_solution, update, -damping_values[index], mCurrentSolution);
    double next_resid_norm = ComputeResidualAndGetNorm(true);
    PrintLineSearchResult(damping_values[index], next_resid_norm);

    index = 2;
    // While f(s_next) < f(s_current), [f = residnorm], keep trying new damping values,
    // ie exit thus loop when next norm of the residual first increases
    while (    (next_resid_norm==DBL_MAX) // the residual is returned as infinity if the deformation is so large to cause exceptions in the material law/EM contraction model
            || ( (next_resid_norm < current_resid_norm) && (index<damping_values.size()) ) )
    {
        current_resid_norm = next_resid_norm;

        // let mCurrentSolution = old_solution - damping_val*update; and compute residual
        VectorSum(old_solution, update, -damping_values[index], mCurrentSolution);
        next_resid_norm = ComputeResidualAndGetNorm(true);
        PrintLineSearchResult(damping_values[index], next_resid_norm);

        index++;
    }

    unsigned best_index;

    if (index==damping_values.size() && (next_resid_norm < current_resid_norm))
    {
        // Difficult to come up with large forces/tractions such that it had to
        // test right down to s=0.05, but overall doesn't fail.
        // The possible damping values have been manually temporarily altered to
        // get this code to be called, it appears to work correctly. Even if it
        // didn't tests wouldn't fail, they would just be v. slightly less efficient.
        #define COVERAGE_IGNORE
        // if we exited because we got to the end of the possible damping values, the
        // best one was the last one (excl the final index++ at the end)
        current_resid_norm = next_resid_norm;
        best_index = index-1;
        #undef COVERAGE_IGNORE
    }
    else
    {
        // else the best one must have been the second last one (excl the final index++ at the end)
        // (as we would have exited when the resid norm first increased)
        best_index = index-2;
    }

    // Check out best was better than the original residual-norm
    if (initial_norm_resid < current_resid_norm)
    {
        #define COVERAGE_IGNORE
        // Have to use an assert/exit here as the following exception causes a seg fault (in cardiac mech problems?)
        // Don't know why
        std::cout << "CHASTE ERROR: (AbstractNonlinearElasticitySolver.hpp): Residual does not appear to decrease in newton direction, quitting.\n" << std::flush;
        exit(0);
        //EXCEPTION("Residual does not appear to decrease in newton direction, quitting");
        #undef COVERAGE_IGNORE
    }

    #ifdef MECH_VERBOSE
    std::cout << "\tBest s = " << damping_values[best_index] << "\n"  << std::flush;
    #endif
    VectorSum(old_solution, update, -damping_values[best_index], mCurrentSolution);

    return current_resid_norm;
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::PostNewtonStep(unsigned counter, double normResidual)
{
}

template<unsigned DIM>
AbstractNonlinearElasticitySolver<DIM>::AbstractNonlinearElasticitySolver(QuadraticMesh<DIM>& rQuadMesh,
                                                                          SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                                                          std::string outputDirectory,
                                                                          CompressibilityType compressibilityType)
    : mrQuadMesh(rQuadMesh),
      mrProblemDefinition(rProblemDefinition),
      mpQuadratureRule(NULL),
      mpBoundaryQuadratureRule(NULL),
      mKspAbsoluteTol(-1),
      mResidualVector(NULL),
      mJacobianMatrix(NULL),
      mPreconditionMatrix(NULL),
      mOutputDirectory(outputDirectory),
      mpOutputFileHandler(NULL),
      mWriteOutputEachNewtonIteration(false),
      mNumNewtonIterations(0),
      mCurrentTime(0.0),
      mCompressibilityType(compressibilityType),
      mCheckedOutwardNormals(false)
{
    assert(DIM==2 || DIM==3);

    if (mCompressibilityType==COMPRESSIBLE)
    {
        mNumDofs = DIM*mrQuadMesh.GetNumNodes();
    }
    else
    {
        mNumDofs = DIM*mrQuadMesh.GetNumNodes() + mrQuadMesh.GetNumVertices();
    }

    mWriteOutput = (mOutputDirectory != "");
    if (mWriteOutput)
    {
        mpOutputFileHandler = new OutputFileHandler(mOutputDirectory);
    }
}

template<unsigned DIM>
AbstractNonlinearElasticitySolver<DIM>::~AbstractNonlinearElasticitySolver()
{
    if (mResidualVector)
    {
        VecDestroy(mResidualVector);
        VecDestroy(mLinearSystemRhsVector);
        MatDestroy(mJacobianMatrix);
        MatDestroy(mPreconditionMatrix);
        if (mCompressibilityType==COMPRESSIBLE)
        {
            VecDestroy(mDirichletBoundaryConditionsVector);
        }
    }

    if (mpQuadratureRule)
    {
        delete mpQuadratureRule;
        delete mpBoundaryQuadratureRule;
    }
    if (mpOutputFileHandler)
    {
        delete mpOutputFileHandler;
    }
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::Solve(double tol,
                                                   unsigned maxNumNewtonIterations,
                                                   bool quitIfNoConvergence)
{
    // If the problem includes specified pressures on deformed surfaces (as opposed
    // to specified tractions), the code needs to compute normals, and they need
    // to be consistently all facing outward (or all facing inward). Check the undeformed
    // mesh boundary elements has nodes that are ordered so that all normals are
    // outward-facing
    if(mrProblemDefinition.GetTractionBoundaryConditionType()==PRESSURE_ON_DEFORMED && mCheckedOutwardNormals==false)
    {
        std::cout << "checking...\n";
        mrQuadMesh.CheckOutwardNormals();
        mCheckedOutwardNormals = true;
    }


    WriteCurrentDeformation("initial");

    if (mWriteOutputEachNewtonIteration)
    {
        WriteCurrentDeformation("newton_iteration", 0);
    }

    // Compute residual
    double norm_resid = ComputeResidualAndGetNorm(false);
    #ifdef MECH_VERBOSE
    std::cout << "\nNorm of residual is " << norm_resid << "\n";
    #endif

    mNumNewtonIterations = 0;
    unsigned iteration_number = 1;

    if (tol < 0) // i.e. if wasn't passed in as a parameter
    {
        tol = NEWTON_REL_TOL*norm_resid;

        #define COVERAGE_IGNORE // not going to have tests in cts for everything
        if (tol > MAX_NEWTON_ABS_TOL)
        {
            tol = MAX_NEWTON_ABS_TOL;
        }
        if (tol < MIN_NEWTON_ABS_TOL)
        {
            tol = MIN_NEWTON_ABS_TOL;
        }
        #undef COVERAGE_IGNORE
    }

    #ifdef MECH_VERBOSE
    std::cout << "Solving with tolerance " << tol << "\n";
    #endif

    while (norm_resid > tol && iteration_number<=maxNumNewtonIterations)
    {
        #ifdef MECH_VERBOSE
        std::cout <<  "\n-------------------\n"
                  <<   "Newton iteration " << iteration_number
                  << ":\n-------------------\n";
        #endif

        // take newton step (and get returned residual)
        norm_resid = TakeNewtonStep();

        #ifdef MECH_VERBOSE
        std::cout << "Norm of residual is " << norm_resid << "\n";
        #endif
        if (mWriteOutputEachNewtonIteration)
        {
            WriteCurrentDeformation("newton_iteration", iteration_number);
        }

        mNumNewtonIterations = iteration_number;

        PostNewtonStep(iteration_number,norm_resid);

        iteration_number++;
        if (iteration_number==20)
        {
            #define COVERAGE_IGNORE
            EXCEPTION("Not converged after 20 newton iterations, quitting");
            #undef COVERAGE_IGNORE
        }
    }

    if ((norm_resid > tol) && quitIfNoConvergence)
    {
        #define COVERAGE_IGNORE
        EXCEPTION("Failed to converge");
        #undef COVERAGE_IGNORE
    }

    // Write the final solution
    WriteCurrentDeformation("solution");
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::WriteCurrentDeformation(std::string fileName, int counterToAppend)
{
    // Only write output if the flag mWriteOutput has been set
    if (!mWriteOutput)
    {
        return;
    }

    std::stringstream file_name;
    file_name << fileName;
    if (counterToAppend >= 0)
    {
        file_name << "_" << counterToAppend;
    }
    file_name << ".nodes";

    out_stream p_file = mpOutputFileHandler->OpenOutputFile(file_name.str());

    std::vector<c_vector<double,DIM> >& r_deformed_position = rGetDeformedPosition();
    for (unsigned i=0; i<r_deformed_position.size(); i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
           * p_file << r_deformed_position[i](j) << " ";
        }
       * p_file << "\n";
    }
    p_file->close();
}

template<unsigned DIM>
unsigned AbstractNonlinearElasticitySolver<DIM>::GetNumNewtonIterations()
{
    return mNumNewtonIterations;
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::SetWriteOutput(bool writeOutput)
{
    if (writeOutput && (mOutputDirectory==""))
    {
        EXCEPTION("Can't write output if no output directory was given in constructor");
    }
    mWriteOutput = writeOutput;
}

template<unsigned DIM>
std::vector<c_vector<double,DIM> >& AbstractNonlinearElasticitySolver<DIM>::rGetDeformedPosition()
{
    mDeformedPosition.resize(mrQuadMesh.GetNumNodes(), zero_vector<double>(DIM));
    for (unsigned i=0; i<mrQuadMesh.GetNumNodes(); i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            mDeformedPosition[i](j) = mrQuadMesh.GetNode(i)->rGetLocation()[j] + mCurrentSolution[DIM*i+j];
        }
    }
    return mDeformedPosition;
}

template<unsigned DIM>
void AbstractNonlinearElasticitySolver<DIM>::CreateCmguiOutput()
{
    if (mOutputDirectory=="")
    {
        EXCEPTION("No output directory was given so no output was written, cannot convert to cmgui format");
    }

    CmguiDeformedSolutionsWriter<DIM> writer(mOutputDirectory + "/cmgui",
                                             "solution",
                                             mrQuadMesh,
                                             WRITE_QUADRATIC_MESH);

    std::vector<c_vector<double,DIM> >& r_deformed_positions = rGetDeformedPosition();
    writer.WriteInitialMesh(); // this writes solution_0.exnode and .exelem
    writer.WriteDeformationPositions(r_deformed_positions, 1); // this writes the final solution as solution_1.exnode
    writer.WriteCmguiScript(); // writes LoadSolutions.com
}

// Constant setting definitions

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::MAX_NEWTON_ABS_TOL = 1e-7;

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::MIN_NEWTON_ABS_TOL = 1e-10;

template<unsigned DIM>
double AbstractNonlinearElasticitySolver<DIM>::NEWTON_REL_TOL = 1e-4;

#endif /*ABSTRACTNONLINEARELASTICITYSOLVER_HPP_*/
