#include "SimpleLinearSolver.hpp"
#include "global/src/Exception.hpp"
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

Vec SimpleLinearSolver::Solve(Mat lhsMatrix, Vec rhsVector, int size)
{
    Vec lhs_vector;
	VecDuplicate(rhsVector, &lhs_vector);

    /* The following lines are very useful for debugging
     *    MatView(lhsMatrix,    PETSC_VIEWER_STDOUT_WORLD);
     *    VecView(rhsVector,    PETSC_VIEWER_STDOUT_WORLD);
     */
    if (mLinearSystemKnown==false){
        PC prec; //Type of pre-conditioner
     
        KSPCreate(PETSC_COMM_WORLD, &mSimpleSolver);
    
        KSPSetOperators(mSimpleSolver, lhsMatrix, lhsMatrix,SAME_NONZERO_PATTERN);
        
        // Default relative tolerance appears to be 1e-5.  This ain't so great.
        KSPSetTolerances(mSimpleSolver, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
        
        //Turn off pre-conditioning if the system size is very small
        KSPGetPC(mSimpleSolver,&prec);
        if (size <= 4)
        {
             PCSetType(prec,PCNONE);
        } else {
             PCSetType(prec,PCJACOBI);        
        }
        
        KSPSetFromOptions(mSimpleSolver) ;
        KSPSetUp(mSimpleSolver);
        
        mLinearSystemKnown=true;
    }    
   
    
    KSPSolve(mSimpleSolver, rhsVector, lhs_vector);
    
    // Check that solver converged and throw if not
    KSPConvergedReason reason;
    KSPGetConvergedReason(mSimpleSolver, &reason);
    if (reason<0)
    {
    	std::stringstream reason_stream;
    	reason_stream << reason;
    	VecDestroy(lhs_vector);    // Delete vec memory, since caller can't do so
    	KSPDestroy(mSimpleSolver); // Likewise for the solver
    	throw Exception("Linear Solver did not converge. Petsc reason code:"
    	                +reason_stream.str()+" .");
    }
   return lhs_vector;
}

