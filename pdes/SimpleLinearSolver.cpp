#include "SimpleLinearSolver.hpp"
#include "petscksp.h"

Vec SimpleLinearSolver::Solve(Mat lhsMatrix, Vec rhsVector)
{
    Vec lhs_vector;
    VecCreate(PETSC_COMM_WORLD, &lhs_vector);
    int rhs_size;

    VecGetSize(rhsVector, &rhs_size);
    VecSetSizes(lhs_vector,PETSC_DECIDE,rhs_size);
    VecSetType(lhs_vector, VECSEQ);

    KSP simple_solver; 
    KSPCreate(PETSC_COMM_WORLD, &simple_solver);
    KSPSetOperators(simple_solver, lhsMatrix, lhsMatrix,SAME_NONZERO_PATTERN);
    KSPSetUp(simple_solver);   
    KSPSolve(simple_solver,rhsVector,lhs_vector);
    return lhs_vector;
}
