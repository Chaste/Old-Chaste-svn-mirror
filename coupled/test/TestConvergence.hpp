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

#include "ConvergenceTester.hpp"

class TestConvergence : public CxxTest::TestSuite
{   
public:


//    void TestMonodomainSpaceAndTime_LuoRudyIModel1991OdeSystem_1D() throw(Exception)
//    {
//        ConvergenceTester<LuoRudyIModel1991OdeSystem, MonodomainProblem<1>, 1> tester;
//    }
    
    void TestMonodomainSpaceAndTime_LuoRudyIModel1991OdeSystem_2D() throw(Exception)
    {
        ConvergenceTester<LuoRudyIModel1991OdeSystem, MonodomainProblem<2>, 2> tester;
    }
    
};

#endif /*TESTCONVERGENCE_HPP_*/
