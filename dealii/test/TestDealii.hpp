#ifndef TESTDEALII_HPP_
#define TESTDEALII_HPP_

#include <cxxtest/TestSuite.h>
#include "DiffusionProblem.cpp"


class TestDealii : public CxxTest::TestSuite
{
public:
    void testDealiiOnDiffusionProblem()
    {
        DiffusionProblem diffusion_problem;
        diffusion_problem.Run();
    }
};

#endif /*TESTDEALII_HPP_*/
