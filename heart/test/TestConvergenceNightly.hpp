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
//#include "KspConvergenceTester.hpp"
#include "OdeConvergenceTester.hpp"
//#include "StimulusConvergenceTester.hpp"

class TestConvergenceNightly : public CxxTest::TestSuite
{   
public:
    //Current test takes about an hour.
    void xTest2DSpace() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2> tester;
        tester.Converge();
        TS_ASSERT(tester.IsConverged());
        //TS_ASSERT_EQUALS(tester.GetMeshNum(), 5u);
        //TS_ASSERT_DELTA(tester.GetSpaceStep(), 1e-4 /*cm*/, 1e-7);
        TS_ASSERT_EQUALS(tester.GetMeshNum(), 5); 
        TS_ASSERT_DELTA(tester.GetSpaceStep(), 1.5625e-3 /*cm*/, 1e-8 /*Allowed error*/);     
    }
    
    //This is much briefer (20mins?)
    void xTest2DSpaceWithSymmLq() throw(Exception)
    {
        PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "bjacobi");
        PetscOptionsSetValue("-options_table", "");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2> tester;
        tester.KspRtol=5e-8;
        tester.Converge();
        //TS_ASSERT(tester.Converged);
        //TS_ASSERT_EQUALS(tester.MeshNum, 5u); 
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.GetMeshNum(), 5); 
        TS_ASSERT_DELTA(tester.GetSpaceStep(), 1.5625e-3 /*cm*/, 1e-8 /*Allowed error*/);     
    }
    
    
    //Currently takes about 3 minutes to do mesh0 and mesh1
    void Test3DSpaceWithSymmLq() throw(Exception)
    {
        PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "bjacobi");
        PetscOptionsSetValue("-options_table", "");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<3>, 3> tester;
        tester.KspRtol=5e-8;
        //tester.RelativeConvergenceCriterion=4e-2;//Just to prove the thing works
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 1u); ///Just to prove the thing works
    }


};

#endif /*TESTCONVERGENCENIGHTLY_HPP_*/
