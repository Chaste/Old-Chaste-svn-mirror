#ifndef TESTCONVERGENCENIGHTLY_HPP_
#define TESTCONVERGENCENIGHTLY_HPP_

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

class TestConvergenceNightly : public CxxTest::TestSuite
{   
public:

    void RunConvergenceTester(AbstractUntemplatedConvergenceTester *pTester, bool stimulateRegion)
    {
            pTester->StimulateRegion=stimulateRegion;
            if ( stimulateRegion )
            {
                pTester->MeshNum = 6u;    
            }
            
            pTester->Converge();
            TS_ASSERT(pTester->Converged);
    }

    void ConvergeInVarious(bool stimulateRegion)
    {
       {
            PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
            RunConvergenceTester(&tester, stimulateRegion);           
            TS_ASSERT_DELTA(tester.PdeTimeStep, 5.0e-3, 1e-10);
        }
    
        {
            SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
            RunConvergenceTester(&tester, stimulateRegion);   
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
            RunConvergenceTester(&tester, stimulateRegion);    
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
            RunConvergenceTester(&tester, stimulateRegion);    
            TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
        }
        
        {
            OdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
            tester.PdeTimeStep=0.01;
            RunConvergenceTester(&tester, stimulateRegion);    
            TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
        }
        
    }
    
public:

    void TestStimulatePlanein1D() throw(Exception)
    {
        ConvergeInVarious(false);
    }

    void TestStimulateRegionin1D() throw(Exception)
    {
        ConvergeInVarious(true);
    }



    //Current test takes about 20 mins.
    //This is much longer (1 hour?) with default ksp
    void Test2DSpaceSymmLq() throw(Exception)
    {
        PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "bjacobi");
        PetscOptionsSetValue("-options_table", "");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2> tester;
        tester.KspRtol=5e-8;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u); 
        TS_ASSERT(tester.IsConverged());
        //TS_ASSERT_EQUALS(tester.GetMeshNum(), 4); 
        //TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0023 /*cm*/, 1e-4 /*Allowed error*/);     
    }

    void Test2DSpaceWithRegion() throw(Exception)
    {
        PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "bjacobi");
        PetscOptionsSetValue("-options_table", "");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2> tester;
        tester.StimulateRegion=true;
        tester.KspRtol=1e-8;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 6u); 
    }       
    
    //Currently takes about 3 minutes to do mesh0 and mesh1
    void Test3DSpaceWithSymmLq() throw(Exception)
    {
        PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "bjacobi");
        PetscOptionsSetValue("-options_table", "");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<3>, 3> tester;
        tester.KspRtol=1e-8;
        tester.RelativeConvergenceCriterion=4e-2;//Just to prove the thing works
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 1u); ///Just to prove the thing works
    }

    //Copied from projects/jmpf
    //Ought to take less than an hour
    void Test3DSpace10() throw(Exception)
    {
        PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "bjacobi");
        PetscOptionsSetValue("-options_table", "");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<3>, 3> tester;
        tester.KspRtol=1e-10;
        tester.OdeTimeStep /= 2.0;
        tester.PdeTimeStep /= 2.0;
        
        //tester.RelativeConvergenceCriterion=2e-11;
        tester.SetMeshWidth(0.10);//cm
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 3u);
    }   
  
};

#endif /*TESTCONVERGENCENIGHTLY_HPP_*/
