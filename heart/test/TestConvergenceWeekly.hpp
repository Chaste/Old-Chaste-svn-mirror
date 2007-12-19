#ifndef TESTCONVERGENCEWEEKLY_HPP_
#define TESTCONVERGENCEWEEKLY_HPP_

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
#include "OdeConvergenceTester.hpp"

class TestConvergenceWeekly : public CxxTest::TestSuite
{   
public:

    void xxTest3DSpace() throw(Exception)
    {
        PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "bjacobi");
        PetscOptionsSetValue("-options_table", "");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<3>, 3> tester;
        tester.KspRtol=1e-8;
        tester.SetMeshWidth(0.15);//cm
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 4u); ///Just to prove the thing works
    }
    
    //Copied from projects/jmpf
    void Test3DSpace10() throw(Exception)
    {
        PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "bjacobi");
        PetscOptionsSetValue("-options_table", "");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<3>, 3> tester;
        tester.KspRtol=1e-10;
        tester.OdeTimeStep /= 2.0;
        tester.PdeTimeStep /= 2.0;
        tester.SetMeshWidth(0.10);//cm
        
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 3u);
    }   
  

};

#endif /*TESTCONVERGENCEWEEKLY_HPP_*/
