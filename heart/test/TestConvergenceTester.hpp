#ifndef TESTCONVERGENCETESTER_HPP_
#define TESTCONVERGENCETESTER_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "BackwardEulerLuoRudyIModel1991.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "PdeConvergenceTester.hpp"
#include "SpaceConvergenceTester.hpp"
#include "KspConvergenceTester.hpp"
#include "OdeConvergenceTester.hpp"
#include "StimulusConvergenceTester.hpp"

class TestConvergenceTester : public CxxTest::TestSuite
{   
public:
    void Test1DOdeTime() throw(Exception)
    {
        OdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
        tester.MeshNum=1;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.OdeTimeStep, 0.0025); 
    }
    
    void Test1DPdeTime() throw(Exception)
    {
        PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1> tester;
        tester.MeshNum=1;
        tester.RelativeConvergenceCriterion=7e-4;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.01); 
    }
    void Test1DPdeTimePlane() throw(Exception)
    {
        PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1> tester;
        tester.MeshNum=1;
        tester.StimulateRegion=true;
        tester.RelativeConvergenceCriterion=5e-4;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.01); 
    }
    void Test1DSpace() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1> tester;
        tester.RelativeConvergenceCriterion=2e-2;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 2u); 
    }
    void Test2DOdeTime() throw(Exception)
    {
        OdeConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<2>, 2> tester;
        tester.MeshNum=0;
        tester.RelativeConvergenceCriterion=8e-5;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.OdeTimeStep, 0.0025); 
    }
    
};

#endif /*TESTCONVERGENCETESTER_HPP_*/
