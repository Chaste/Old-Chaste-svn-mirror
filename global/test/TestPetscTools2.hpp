#ifndef TESTPETSCTOOLS_HPP_
#define TESTPETSCTOOLS_HPP_

#include <cxxtest/TestSuite.h>

#include "PetscTools.hpp"

class TestPetscTools : public CxxTest::TestSuite
{
public:
    void TestBarrier()
    {
        // Testing the barrier method is kind of tricky, since we really want
        // to also check if it works when PETSc isn't set up.  Hence this file.
        PetscTools::Barrier();
    }
};
#endif /*TESTPETSCTOOLS_HPP_*/
