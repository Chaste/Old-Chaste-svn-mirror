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

    void ConvergeInVarious(bool stimulateRegion)
    {
       {
            TimeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
            tester.StimulateRegion=stimulateRegion;
            tester.Converge();
            TS_ASSERT(tester.Converged);
            TS_ASSERT_DELTA(tester.PdeTimeStep, 5.0e-3, 1e-10);
        }
    
        {
            SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
            tester.StimulateRegion=stimulateRegion;
            tester.Converge();
            TS_ASSERT(tester.Converged);
            if (!stimulateRegion)
            {
                TS_ASSERT_EQUALS(tester.MeshNum, 5u); 
            }
            else
            {
                TS_ASSERT_EQUALS(tester.MeshNum, 6u);
            }
        }
            
        {
            KspConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
            tester.StimulateRegion=stimulateRegion;
            tester.Converge();
            TS_ASSERT(tester.Converged);
            if (!stimulateRegion)
            {
                TS_ASSERT_DELTA(tester.KspRtol, 1e-5, 1e-10);
            }
            else
            {
                TS_ASSERT_DELTA(tester.KspRtol, 1e-6, 1e-10);
            }
        }
    
        {
            OdeConvergenceTester<LuoRudyIModel1991OdeSystem, BidomainProblem<1>, 1> tester;
            tester.StimulateRegion=stimulateRegion;
            tester.PdeTimeStep=0.01;
            tester.Converge();
            TS_ASSERT(tester.Converged);
            TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
        }
        
        {
            OdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
            tester.StimulateRegion=stimulateRegion;
            tester.PdeTimeStep=0.01;
            tester.Converge();
            TS_ASSERT(tester.Converged);
            TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
        }
        
    }
    
    void TestStimulatePointAndRegion() throw(Exception)
    {
        ConvergeInVarious(false);
        ConvergeInVarious(true);
    }
 
};

#endif /*TESTCONVERGENCE_HPP_*/
