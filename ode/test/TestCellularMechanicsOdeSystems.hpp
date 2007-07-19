#ifndef TESTCELLULARMECHANICSODESYSTEMS_HPP_
#define TESTCELLULARMECHANICSODESYSTEMS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include "NHSCellularMechanicsOdeSystem.hpp"
#include "EulerIvpOdeSolver.hpp"

class TestCellularMechanicsOdeSystems : public CxxTest::TestSuite
{
public :
    void TestNHSCellularMechanicsOdeSystem() throw(Exception)
    {
        NHSCellularMechanicsOdeSystem nhs_system;
        
        // Hardcoded results for two values for z when lambda1=0.
        // Note: CalculateT0(z) is a private method.
        TS_ASSERT_DELTA(nhs_system.CalculateT0(0), 0, 1e-12);
        TS_ASSERT_DELTA(nhs_system.CalculateT0(1), 58.0648, 1e-3);
        
        EulerIvpOdeSolver euler_solver;

        // lambda1=1, dlamdt = 0, so there should be no active tension
        euler_solver.SolveAndUpdateStateVariable(&nhs_system, 0, 1, 0.01);
        TS_ASSERT_DELTA(nhs_system.GetActiveTension(), 0.0, 1e-12);
        
// todo: sensible test needed..
//        nhs_system.SetLambda1DerivativeAndCalciumI(1.1, 0.1, 0.1);
//        euler_solver.SolveAndUpdateStateVariable(&nhs_system, 0, 1, 0.01);
//        TS_ASSERT_DELTA(nhs_system.GetActiveTension(), 0.0, 1e-12);
    }
};
#endif /*TESTCELLULARMECHANICSODESYSTEMS_HPP_*/
