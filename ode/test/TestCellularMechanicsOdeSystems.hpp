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
        
        // todo: verify these are correct somehow...
        nhs_system.SetLambda1DerivativeAndCalciumI(0.8, -0.1, 0.1);
        OdeSolution solution = euler_solver.Solve(&nhs_system, nhs_system.rGetStateVariables(), 0, 1, 0.001,0.001);
        solution.WriteToFile("CellularMechanics/NHS/", "lam_is_0.8", &nhs_system, "ms");

        TS_ASSERT_DELTA(nhs_system.GetActiveTension(), 0.0012, 1e-3);
        
        nhs_system.SetLambda1DerivativeAndCalciumI(0.5, 0.1, 0.1);
        solution = euler_solver.Solve(&nhs_system, nhs_system.rGetStateVariables(), 0, 1, 0.001,0.001);
        solution.WriteToFile("CellularMechanics/NHS/", "lam_is_0.5", &nhs_system, "ms", 1, false);

        TS_ASSERT_DELTA(nhs_system.GetActiveTension(), -0.1120, 1e-3);        
    }
};
#endif /*TESTCELLULARMECHANICSODESYSTEMS_HPP_*/
