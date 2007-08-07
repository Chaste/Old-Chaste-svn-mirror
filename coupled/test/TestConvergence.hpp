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

    
    void Test1() throw(Exception)
    {
        ConvergenceTester<LuoRudyIModel1991OdeSystem, BidomainProblem<2>, 2> tester;
    }
    
};

#endif /*TESTCONVERGENCE_HPP_*/
