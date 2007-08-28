#ifndef TESTCONVERGENCE_HPP_
#define TESTCONVERGENCE_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "BackwardEulerLuoRudyIModel1991.hpp"
#include "SpaceConvergenceTester.hpp"

class TestSpaceConvergence : public CxxTest::TestSuite
{   
public:

    
    void Test1DSpace() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
        tester.Converge();
        TS_ASSERT(tester.mConverged);
        TS_ASSERT_EQUALS(tester.mMeshNum, 5u); 
    }

 
};

#endif /*TESTCONVERGENCE_HPP_*/
