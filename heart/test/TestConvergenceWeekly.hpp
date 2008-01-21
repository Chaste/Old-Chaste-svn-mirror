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
        //PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "jacobi");
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
    
    //Experiments with ksp_atol follow.  Note that the way that
    //tester.KspRtol=1e-17; is set is a bit hacky, but it'll do for now.
       
    void TestSpaceConvergencein1DWithAtol() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
        tester.KspRtol=1e-17;//Unreasonably small
        
        PetscOptionsSetValue("-ksp_atol", "1e-5");
        PetscOptionsSetValue("-options_table", "");
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 1.68417e-05);
        //Has to be at least as good as the 1D with Rtol=1e-7
        //Note the final line fails with ksp_atol=1e-4
    }
    
    void TestSpaceConvergencein2DWithAtol() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2> tester;
        tester.KspRtol=1e-17;//Unreasonably small
        
        
        PetscOptionsSetValue("-ksp_atol", "1e-5");
        PetscOptionsSetValue("-options_table", "");
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 6.65582e-05);
        //Comes in at 1.17118e-5
        //Has to be at least as good as the 2D with Rtol=5e-8
        
    }
    
    //Copied from projects/jmpf since this converges on mesh4
    void Test3DSpaceRelaxWidthWithAtol() throw(Exception)
    {
        
        PetscOptionsSetValue("-ksp_atol", "1e-3");
        PetscOptionsSetValue("-options_table", "");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<3>, 3> tester;
        tester.KspRtol=1e-17;
        
        
        tester.SetMeshWidth(0.15);//cm
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 4u);
    }
  

};

#endif /*TESTCONVERGENCEWEEKLY_HPP_*/
