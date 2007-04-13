#ifndef PARALLELVECTOR_HPP_
#define PARALLELVECTOR_HPP_

#include <vector>
#include <petscvec.h>
#include <iostream>

#include "ParallelIterator.hpp"
#include "ParallelProblem.hpp"

/***
 * Gives access to the local portion of a PETSc vector via an iterator.
 */


class ParallelVector
{
private:
    Vec mVec;
public:
    int mStride;
    double *mpVec;
    /***
     * Constructor
     * @param vec PETSc vector of which this class shall be a portion.
     */
    ParallelVector(Vec vec);
    
    double& operator[](ParallelIterator index);
    
    /**
     * Store elements that have been written to via the slice
     * back into the PETSc vector. Call after you have finished writing.
     */
    void Restore();
};



#endif /*PARALLELVECTOR_HPP_*/
