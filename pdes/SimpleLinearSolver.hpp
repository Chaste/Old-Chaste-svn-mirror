#ifndef SIMPLELINEARSOLVER_H
#define SIMPLELINEARSOLVER_H
#include "AbstractLinearSolver.hpp"
#include "petscvec.h"
#include "petscmat.h"

class SimpleLinearSolver : public AbstractLinearSolver
{

public:

    Vec Solve(Mat lhsMatrix, Vec rhsVector);


};

#endif
