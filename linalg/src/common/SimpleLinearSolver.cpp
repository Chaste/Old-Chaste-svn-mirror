#include "SimpleLinearSolver.hpp"
#include "Exception.hpp"
#include "PetscException.hpp"
#include <sstream>

/**
 * \todo Document class + exceptional behaviour.
 *
 * This method can solve linear systems of the form Ax = b
 *
 * @param lhsMatrix A
 * @param rhsVector b
 * @return The solution Vec x.
 */

Vec SimpleLinearSolver::Solve(Mat lhsMatrix, Vec rhsVector, int size, MatNullSpace matNullSpace)
{
    Vec lhs_vector;
    VecDuplicate(rhsVector, &lhs_vector);
    
    /* The following lines are very useful for debugging
     *    MatView(lhsMatrix,    PETSC_VIEWER_STDOUT_WORLD);
     *    VecView(rhsVector,    PETSC_VIEWER_STDOUT_WORLD);
     */
    if (mLinearSystemKnown==false)
    {
        PC prec; //Type of pre-conditioner
        
        KSPCreate(PETSC_COMM_WORLD, &mSimpleSolver);
        //See
        //http://www-unix.mcs.anl.gov/petsc/petsc-2/snapshots/petsc-current/docs/manualpages/KSP/KSPSetOperators.html
        //The preconditioner flag (last argument) in the following calls says
        //how to reuse the preconditioner on subsequent iterations
        if (mMatrixIsConstant==true)
        {
            KSPSetOperators(mSimpleSolver, lhsMatrix, lhsMatrix,SAME_PRECONDITIONER);
            
        }
        else
        {
            KSPSetOperators(mSimpleSolver, lhsMatrix, lhsMatrix,SAME_NONZERO_PATTERN);
        }
        
        // Default relative tolerance appears to be 1e-5.  This ain't so great.
        KSPSetTolerances(mSimpleSolver, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
        
        //Turn off pre-conditioning if the system size is very small
        KSPGetPC(mSimpleSolver,&prec);
        if (size <= 4)
        {
            PCSetType(prec,PCNONE);
        }
        else
        {
            PCSetType(prec,PCJACOBI);
        }
        
        if (matNullSpace)
        {
            PETSCEXCEPT( KSPSetNullSpace(mSimpleSolver, matNullSpace) );
        }
        
        KSPSetFromOptions(mSimpleSolver) ;
        KSPSetUp(mSimpleSolver);
        
        mLinearSystemKnown=true;
    }
    
    PETSCEXCEPT(KSPSolve(mSimpleSolver, rhsVector, lhs_vector));
    
    // Check that solver converged and throw if not
    KSPConvergedReason reason;
    KSPGetConvergedReason(mSimpleSolver, &reason);
    
    KSPEXCEPT(reason);
    
    return lhs_vector;
}

