#ifndef TESTPROBLEMPORTION_HPP_
#define TESTPROBLEMPORTION_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"

#include "ParallelProblem.hpp"

#include "GlobalParallelProblem.hpp"

class TestParallelProblem : public CxxTest::TestSuite
{
public:
    
    void TestBeginEndAndSize()
    {
        gProblem.SetProblemSize(10);
        
        // create a 10 element petsc vector
        Vec vec;

        VecCreate(PETSC_COMM_WORLD, &vec);
        VecSetSizes(vec, PETSC_DECIDE, 10);
        VecSetFromOptions(vec);
        
        // calculate my range
        PetscInt petsc_lo, petsc_hi;
        VecGetOwnershipRange(vec,&petsc_lo,&petsc_hi);
        unsigned lo=(unsigned)petsc_lo;
        unsigned hi=(unsigned)petsc_hi;
        
        // check that the problem portion agrees
        TS_ASSERT_EQUALS( hi - lo, gProblem.Size() );
                          
        // check that iterator, local and global work
        
        ParallelIterator index = gProblem.Begin();
        unsigned global_index = lo;
        
        while (index != gProblem.End()
               && global_index != hi)
        {
            unsigned local_index = global_index-lo;
            TS_ASSERT_EQUALS(gProblem.Global(index), global_index);
            TS_ASSERT_EQUALS(gProblem.Local(index), local_index);
            ++index;
            global_index++;
        }
        TS_ASSERT_EQUALS(index, gProblem.End());
        TS_ASSERT_EQUALS(global_index, hi);
    } 
};

#endif /*TESTPROBLEMPORTION_HPP_*/
