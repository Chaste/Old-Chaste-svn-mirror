#include "SimpleLinearSolver.hpp"
#include "petscksp.h"
#include "Exception.hpp"
#include <sstream>
/**
 *  TODO: Document class + exceptional behaviour
 * 
 * 
 * 
 */

#include <iostream>
Vec SimpleLinearSolver::Solve(Mat lhsMatrix, Vec rhsVector)
{
    Vec lhs_vector;
    VecCreate(PETSC_COMM_WORLD, &lhs_vector);
    int rhs_size;

    VecGetSize(rhsVector, &rhs_size);
    VecSetSizes(lhs_vector,PETSC_DECIDE,rhs_size);
    //VecSetType(lhs_vector, VECSEQ);
	//VecSetType(lhs_vector, VECMPI);
	VecSetFromOptions(lhs_vector);
    //std::cout<< "KSP about to be called \n " ;

    KSP simple_solver; 
    KSPCreate(PETSC_COMM_WORLD, &simple_solver);
    
    KSPSetOperators(simple_solver, lhsMatrix, lhsMatrix,SAME_NONZERO_PATTERN);
    //std::cout<<" Setup \n" ;
    KSPSetFromOptions(simple_solver) ;
    KSPSetUp(simple_solver);   
    
    
    KSPSolve(simple_solver,rhsVector,lhs_vector);
    
    
    // Check that solver converged and throw if not
    KSPConvergedReason reason;
    KSPGetConvergedReason(simple_solver, &reason);
    if (reason<0)
    {
    	std::stringstream reason_stream;
    	reason_stream << reason;
    	throw Exception("Linear Solver did not converge. Petsc reason code:"
    	                +reason_stream.str()+" .");
    }
    KSPDestroy(simple_solver) ;
    return lhs_vector;
    
}
