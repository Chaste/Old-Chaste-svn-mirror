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
#include "PdeConvergenceTester.hpp"
#include "SpaceConvergenceTester.hpp"
#include "KspConvergenceTester.hpp"
#include "OdeConvergenceTester.hpp"
//#include "StimulusConvergenceTester.hpp"

class TestConvergenceTester : public CxxTest::TestSuite
{   
public:
    void Test1DOdeTime() throw(Exception)
    {
        OdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
        tester.MeshNum=1;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.OdeTimeStep, 0.0025); 
    }
    
    void Test1DPdeTime() throw(Exception)
    {
        PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.RelativeConvergenceCriterion=7e-4;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.01); 
    }
    void Test1DPdeTimePlane() throw(Exception)
    {
        PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.StimulateRegion=true;
        tester.RelativeConvergenceCriterion=5e-4;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.01); 
    }
    void Test1DSpace() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1, 1> tester;
        tester.RelativeConvergenceCriterion=2e-2;
        tester.Converge();
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.GetMeshNum(), 2); 
        TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0125, 1e-8); 
        //TS_ASSERT_EQUALS(tester.GetMeshNum(), 5); 
        //TS_ASSERT_DELTA(tester.GetSpaceStep(), 1.5625e-3, 1e-8); 
    }
    void Test2DOdeTime() throw(Exception)
    {
        OdeConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<2>, 2, 1> tester;
        tester.MeshNum=0;
        tester.RelativeConvergenceCriterion=8e-5;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.OdeTimeStep, 0.0025); 
    }
    
    void TestSpaceConvergencein1DWithRelativeTolerance() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
        tester.SetKspRelativeTolerance(1e-4);
        TS_ASSERT_DELTA(tester.GetKspRelativeTolerance(), 1e-4, 1e-10);
        TS_ASSERT_THROWS_ANYTHING(tester.GetKspAbsoluteTolerance());
        tester.RelativeConvergenceCriterion=2e-2;
        tester.Converge();
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.GetMeshNum(), 2); 
        TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0125, 1e-8);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 3.1e-3);
    }
    void TestSpaceConvergencein1DWithAbsoluteTolerance() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
        tester.SetKspAbsoluteTolerance(1e-4);
        TS_ASSERT_DELTA(tester.GetKspAbsoluteTolerance(), 1e-4, 1e-10);
        TS_ASSERT_THROWS_ANYTHING(tester.GetKspRelativeTolerance());
        tester.RelativeConvergenceCriterion=2e-2;
        tester.Converge();
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.GetMeshNum(), 2); 
        TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0125, 1e-8);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 2.9e-3);
     }
    
};

#endif /*TESTCONVERGENCETESTER_HPP_*/
