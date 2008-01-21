#ifndef TESTPERFORMANCE_HPP_
#define TESTPERFORMANCE_HPP_

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
#include "PerformanceTester.hpp"

class TestPerformance : public CxxTest::TestSuite
{   
public:

    
    void TestPerf() throw(Exception)
    {
        // solver and preconditioner options
        //PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "jacobi");
        PetscOptionsSetValue("-options_table", "");
        
        // write headings
        PerformanceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<3>, 3>::DisplayHeadings();
        EventHandler::Headings();

        // base line
        PerformanceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<3>, 3> tester;
        tester.MeshNum=4;
        tester.Run();
        
        EventHandler::Report();
    }
};

#endif /*TESTPERFORMANCE_HPP_*/
