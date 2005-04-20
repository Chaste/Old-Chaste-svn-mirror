#ifndef ABSTRACTLINEARSOLVER_H
#define ABSTRACTLINEARSOLVER_H

#include "petscvec.h"
#include "AbstractLinearSolver.hpp"
#include "petscmat.h"

class AbstractLinearSolver
{

public:

    virtual Vec Solve(Mat lhsMatrix, Vec rhsVector) = 0;


};

#endif

