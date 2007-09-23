#ifndef TESTCELLULARMECHANICSODESYSTEMS_HPP_
#define TESTCELLULARMECHANICSODESYSTEMS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include "NhsCellularMechanicsOdeSystem.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "ZeroStimulus.hpp"

class TestCellularMechanicsOdeSystems : public CxxTest::TestSuite
{
public :
    void TestNhsCellularMechanicsOdeSystem() throw(Exception)
    {
        NhsCellularMechanicsOdeSystem nhs_system;
        
        // Hardcoded results for two values for z when lambda1=0.
        // Note: CalculateT0(z) is a private method.
        TS_ASSERT_DELTA(nhs_system.CalculateT0(0), 0, 1e-12);
        TS_ASSERT_DELTA(nhs_system.CalculateT0(1), 58.0648, 1e-3);
        
        EulerIvpOdeSolver euler_solver;

        // the following is just to get a realistic Ca_I value
        ZeroStimulus zero_stimulus;
        LuoRudyIModel1991OdeSystem lr91(&euler_solver, 0.01, &zero_stimulus);
        unsigned Ca_i_index = lr91.GetStateVariableNumberByName("CaI");
        double Ca_I = lr91.rGetStateVariables()[Ca_i_index];
        
        // lambda1=1, dlamdt = 0, so there should be no active tension
        euler_solver.SolveAndUpdateStateVariable(&nhs_system, 0, 1, 0.01);
        TS_ASSERT_DELTA(nhs_system.GetActiveTension(), 0.0, 1e-12);
        
        // the following test doesn't make sense, as lambda=const, but dlam_dt > 0
        // but it is not possible to test the NHS system by itself without having a varying 
        // lambda and get non-trivial solutions. So, we'll have a non-realistic
        // test here, and TS_ASSERT against hardcoded values, just to check nothing
        // has changed. A proper test where lambda varies (which means time-looping has
        // to be done outside the solver is done in TestElectroMechanicCellularModels, 
        // where NHS is coupled to a cell model
        nhs_system.SetLambdaAndDerivative(0.5, 0.1);
        TS_ASSERT_DELTA(nhs_system.GetLambda(), 0.5, 1e-12);
        nhs_system.SetIntracellularCalciumConcentration(Ca_I);
        OdeSolution solution = euler_solver.Solve(&nhs_system, nhs_system.rGetStateVariables(), 0, 10, 0.01, 0.01);

        unsigned num_timesteps = solution.GetNumberOfTimeSteps();
        TS_ASSERT_DELTA( solution.rGetSolutions()[num_timesteps-1][0],   0.0056, 1e-2 );    
        TS_ASSERT_DELTA( solution.rGetSolutions()[num_timesteps-1][1],   0.0000, 1e-2 );    
        TS_ASSERT_DELTA( solution.rGetSolutions()[num_timesteps-1][2], -25.0359, 1e-2 );    
        TS_ASSERT_DELTA( solution.rGetSolutions()[num_timesteps-1][3],  77.2103, 1e-2 );    
        TS_ASSERT_DELTA( solution.rGetSolutions()[num_timesteps-1][4],  20.6006, 1e-2 );    
    }
};
#endif /*TESTCELLULARMECHANICSODESYSTEMS_HPP_*/
