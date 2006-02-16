#include "SimpleLinearSolver.hpp"
#include "petscksp.h"
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
     
    KSP simple_solver;
    PC prec; //Type of pre-conditioner
     
    KSPCreate(PETSC_COMM_WORLD, &simple_solver);
    
    KSPSetOperators(simple_solver, lhsMatrix, lhsMatrix,SAME_NONZERO_PATTERN);
    // Default relative tolerance appears to be 1e-5.  This ain't so great.
    KSPSetTolerances(simple_solver, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    
    //Turn off pre-conditioning if the system size is very small
    KSPGetPC(simple_solver,&prec);
    if (size <= 4)
    {
         PCSetType(prec,PCNONE);
    } else {
         PCSetType(prec,PCJACOBI);        
    }
    
    KSPSetFromOptions(simple_solver) ;
    KSPSetUp(simple_solver);   
    
    KSPSolve(simple_solver, rhsVector, lhs_vector);
    
    // Check that solver converged and throw if not
    KSPConvergedReason reason;
    KSPGetConvergedReason(simple_solver, &reason);
    if (reason<0)
    {
    	std::stringstream reason_stream;
    	reason_stream << reason;
    	VecDestroy(lhs_vector);    // Delete vec memory, since caller can't do so
    	KSPDestroy(simple_solver); // Likewise for the solver
    	throw Exception("Linear Solver did not converge. Petsc reason code:"
    	                +reason_stream.str()+" .");
    }
    KSPDestroy(simple_solver);
    return lhs_vector;
}
