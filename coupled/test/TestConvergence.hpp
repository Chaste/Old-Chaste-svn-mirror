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
#include "TimeConvergenceTester.hpp"
#include "SpaceConvergenceTester.hpp"
#include "KspConvergenceTester.hpp"


class TestConvergence : public CxxTest::TestSuite
{   
public:

    
    void Test1DTime() throw(Exception)
    {
        TimeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
        tester.Converge();
        TS_ASSERT(tester.mConverged);
        TS_ASSERT_DELTA(tester.mPdeTimeStep, 2.5e-3, 1e-10);
    }
    
    void Test1DSpace() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
        tester.Converge();
        TS_ASSERT(tester.mConverged);
        TS_ASSERT_EQUALS(tester.mMeshNum, 5u); 
    }
    
    void Test1DKsp() throw(Exception)
    {
        KspConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
        tester.Converge();
        TS_ASSERT(tester.mConverged);
        TS_ASSERT_DELTA(tester.mKspRtol, 1e-5, 1e-10);
    }

 
};

#endif /*TESTCONVERGENCE_HPP_*/
