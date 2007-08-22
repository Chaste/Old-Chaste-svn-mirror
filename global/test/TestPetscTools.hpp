#ifndef TESTPETSCTOOLS_HPP_
#define TESTPETSCTOOLS_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>
#include "PetscTools.hpp"

class TestPetscTools : public CxxTest::TestSuite
{
public:
    void testPetscTools()
    {
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        bool is_sequential = (num_procs==1);
        TS_ASSERT_EQUALS( PetscTools::IsSequential(), is_sequential);

        ////////////////////////////////////////////////////
        // test CreateVec which returns a vec of constants
        ////////////////////////////////////////////////////
        Vec vec1 = PetscTools::CreateVec(10, 3.41);
        ReplicatableVector vec1_repl(vec1);
        
        TS_ASSERT_EQUALS(vec1_repl.size(), 10u);
        for(unsigned i=0; i<10; i++)
        {
            TS_ASSERT_DELTA(vec1_repl[i], 3.41, 1e-12);
        }

        ////////////////////////////////////////////////////
        // test CreateVec which uses a std::vector of data
        ////////////////////////////////////////////////////
        std::vector<double> data(10);
        for(unsigned i=0; i<10; i++)
        {
            data[i] = i+0.45;
        }
        
        Vec vec2 = PetscTools::CreateVec(data);
        ReplicatableVector vec2_repl(vec2);
        
        TS_ASSERT_EQUALS(vec2_repl.size(), 10u);
        for(unsigned i=0; i<10; i++)
        {
            TS_ASSERT_DELTA(vec2_repl[i], i+0.45, 1e-12);
        }
    }
};
#endif /*TESTPETSCTOOLS_HPP_*/
