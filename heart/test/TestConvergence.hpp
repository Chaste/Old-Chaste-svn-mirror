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
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "TimeConvergenceTester.hpp"
#include "SpaceConvergenceTester.hpp"
#include "StimulusConvergenceTester.hpp"
#include "KspConvergenceTester.hpp"
#include "OdeConvergenceTester.hpp"


class TestConvergence : public CxxTest::TestSuite
{   
public:

    
    void Test1DTime() throw(Exception)
    {
        TimeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.PdeTimeStep, 5.0e-3, 1e-10);
     }
    
    void xTest1DSpaceWithVariousKsp() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
        //tester.Converge();
        //TS_ASSERT(tester.Converged);
        //TS_ASSERT_EQUALS(tester.MeshNum, 5u); 
        tester.RelativeConvergenceCriterion=1e-7;
        
        for (int i=2; i<10; i++)
        {
            tester.KspRtol=pow(10,-i);
            tester.Converged=false;
            std::cout<<"###############Gnu new run \n#Gnu KSP = "<< tester.KspRtol<<"\n";
            tester.Converge();
            TS_ASSERT(tester.Converged);
        }
        
    }
    
    void xTest1DStimulus() throw (Exception)
    {
        StimulusConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
        int i=0;
        //for (int i=0; i<6;i++)
        {
            tester.FirstMesh=i;
            std::cout<<"###############Gnu new run \n#Gnu First mesh = "<< tester.FirstMesh<<"\n";
            tester.Converge();
            tester.PopulatedResult=false;
        }
     }
    
    void xTest1DSpace() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u); 
    }
//        
//    void Test1DKsp() throw(Exception)
//    {
//        KspConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
//        tester.Converge();
//        TS_ASSERT(tester.Converged);
//        TS_ASSERT_DELTA(tester.KspRtol, 1e-5, 1e-10);
//    }
//
//    void Test1DOdeForward() throw(Exception)
//    {
//        OdeConvergenceTester<LuoRudyIModel1991OdeSystem, BidomainProblem<1>, 1> tester;
//        tester.Converge();
//        TS_ASSERT(tester.Converged);
//        TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
//    }
//    
//    void Test1DOdeBackward() throw(Exception)
//    {
//        OdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
//        tester.Converge();
//        TS_ASSERT(tester.Converged);
//        TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
//    }
 
};

#endif /*TESTCONVERGENCE_HPP_*/
