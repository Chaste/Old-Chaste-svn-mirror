#ifndef ABSTRACTLINEARSOLVER_H
#define ABSTRACTLINEARSOLVER_H

#include <petscvec.h>
#include <petscmat.h>

/**
 * Abstract base class for solvers for linear systems.
 */
class AbstractLinearSolver
{

public:

	/**
	 * Solve the system lhsMatrix * x = rhsVector for x.
	 * 
	 * @param lhsMatrix The left hand side matrix
	 * @param rhsVector The right hand side vector
	 * @return The solution x.
	 */
    virtual Vec Solve(Mat lhsMatrix, Vec rhsVector, int size) = 0;
    virtual void SetMatrixIsConstant()=0;
};

#endif

