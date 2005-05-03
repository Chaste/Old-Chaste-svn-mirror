#ifndef _SIMPLELINEARSOLVER_H_
#define _SIMPLELINEARSOLVER_H_

#include "AbstractLinearSolver.hpp"
#include "petscvec.h"
#include "petscmat.h"

class SimpleLinearSolver : public AbstractLinearSolver
{

public:

    Vec Solve(Mat lhsMatrix, Vec rhsVector);

};

#endif // _SIMPLELINEARSOLVER_H_
