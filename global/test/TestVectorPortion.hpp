#ifndef TESTVECTORPORTION_HPP_
#define TESTVECTORPORTION_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"

#include "VectorPortion.hpp"

class TestVectorPortion : public CxxTest::TestSuite
{
public:
    
    void TestReadAndRestore()
    {
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
        
        // stuff the vector with values: global index * local index
        double* p_vec;
        VecGetArray(vec, &p_vec);
        
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            p_vec[local_index] = local_index*global_index;
        }
        VecRestoreArray(vec, &p_vec);
        VecAssemblyBegin(vec);
        VecAssemblyEnd(vec);
        
        // create vec portion
        VectorPortion vec_portion(vec);
        
        // check that the portion has the correct range
        TS_ASSERT_EQUALS(vec_portion.Begin().Global,lo);
        TS_ASSERT_EQUALS(vec_portion.End().Global,hi);
        
        for (VectorPortion::Iterator index = vec_portion.Begin();
             index!= vec_portion.End();
             ++index)
        {
            TS_ASSERT_EQUALS(*index, index.Local*index.Global);
        }
        
        // now write something to the vector - local_index * global_index
        for (VectorPortion::Iterator index = vec_portion.Begin();
             index!= vec_portion.End();
             ++index)
        {
            *index =  -(double)(index.Local*index.Global);
        }
        
        TS_ASSERT( ! (vec_portion.Begin() == vec_portion.End()) );        
        
        // ask the portion to restore the main vector
        vec_portion.Restore();
        
        // check that the vector was restored
        VecGetArray(vec, &p_vec);
        
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            TS_ASSERT_EQUALS(p_vec[local_index], -(double)local_index*global_index);
        }

    }
           
        
};

#endif /*TESTVECTORPORTION_HPP_*/
