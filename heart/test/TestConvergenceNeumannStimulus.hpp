#ifndef TESTCONVERGENCENEUMANNSTIMULUS_HPP_
#define TESTCONVERGENCENEUMANNSTIMULUS_HPP_

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
//#include "StimulusConvergenceTester.hpp"

class TestConvergenceNeumannStimulus : public CxxTest::TestSuite
{   
public:

    void TestConvergenceMonodomain1d()
    {
        PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1, 1> tester;
        tester.Stimulus=NEUMANN;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.PdeTimeStep, 5.0e-3, 1e-10);
    }

    void TestConvergenceBidomain1d()
    {
        PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
        tester.Stimulus=NEUMANN;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.PdeTimeStep, 5.0e-3, 1e-10);
    }

    void TestSpaceConvergence1d()
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1, 1> tester;
        tester.Stimulus=NEUMANN;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u); 
    }
    
    void TestSpaceConvergence2d()
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<2>, 2, 1> tester;
        tester.Stimulus=NEUMANN;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u); 
    }
    
    void TestSpaceConvergence2dBidomain()
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2, 2> tester;
        tester.Stimulus=NEUMANN;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u); 
    }
    
    void xTestSpaceConvergence3d()
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<3>, 3, 2> tester;
        tester.SetKspAbsoluteTolerance(1e-3); 
        tester.Stimulus=NEUMANN;
        tester.SetMeshWidth(0.15);//cm
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 4u);
        EventHandler::Headings();
        EventHandler::Report();
    }
};

#endif /*TESTCONVERGENCENEUMANNSTIMULUS_HPP_*/
