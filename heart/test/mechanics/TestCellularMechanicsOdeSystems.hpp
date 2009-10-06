/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef TESTCELLULARMECHANICSODESYSTEMS_HPP_
#define TESTCELLULARMECHANICSODESYSTEMS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include "NhsCellularMechanicsOdeSystem.hpp"
#include "Kerchoffs2003ContractionModel.hpp"
#include "Nash2004ContractionModel.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "ZeroStimulus.hpp"
#include "TimeStepper.hpp"

class TestCellularMechanicsOdeSystems : public CxxTest::TestSuite
{
public :
    void TestNhsCellularMechanicsOdeSystem() throw(Exception)
    {
        NhsCellularMechanicsOdeSystem nhs_system;
        TS_ASSERT_EQUALS(nhs_system.IsStretchDependent(), true);
        TS_ASSERT_EQUALS(nhs_system.IsStretchRateDependent(), true);

        // Hardcoded results for two values for z when lambda1=0.
        // Note: CalculateT0(z) is a private method.
        TS_ASSERT_DELTA(nhs_system.CalculateT0(0), 0, 1e-12);
        TS_ASSERT_DELTA(nhs_system.CalculateT0(1), 58.0648, 1e-3);

        boost::shared_ptr<EulerIvpOdeSolver> p_euler_solver(new EulerIvpOdeSolver);

        // the following is just to get a realistic Ca_I value
        boost::shared_ptr<ZeroStimulus> p_zero_stimulus(new ZeroStimulus);
        LuoRudyIModel1991OdeSystem lr91(p_euler_solver, p_zero_stimulus);
        unsigned Ca_i_index = lr91.GetStateVariableNumberByName("CaI");
        double Ca_I = lr91.rGetStateVariables()[Ca_i_index];

        // lambda1=1, dlamdt = 0, so there should be no active tension
        p_euler_solver->SolveAndUpdateStateVariable(&nhs_system, 0, 1, 0.01);
        TS_ASSERT_DELTA(nhs_system.GetActiveTension(), 0.0, 1e-12);

        // the following test doesn't make sense, as lambda=const, but dlam_dt > 0
        // but it is not possible to test the NHS system by itself without having a varying
        // lambda and get non-trivial solutions. So, we'll have a non-realistic
        // test here, and TS_ASSERT against hardcoded values, just to check nothing
        // has changed. A proper test where lambda varies (which means time-looping has
        // to be done outside the solver is done in TestElectroMechanicCellularModels,
        // where NHS is coupled to a cell model
        nhs_system.SetStretchAndStretchRate(0.5, 0.1);
        TS_ASSERT_DELTA(nhs_system.GetLambda(), 0.5, 1e-12);
        nhs_system.SetIntracellularCalciumConcentration(Ca_I);
        OdeSolution solution = p_euler_solver->Solve(&nhs_system, nhs_system.rGetStateVariables(), 0, 10, 0.01, 0.01);

        unsigned num_timesteps = solution.GetNumberOfTimeSteps();
        TS_ASSERT_DELTA( solution.rGetSolutions()[num_timesteps-1][0],   0.0056, 1e-2 );
        TS_ASSERT_DELTA( solution.rGetSolutions()[num_timesteps-1][1],   0.0000, 1e-2 );
        TS_ASSERT_DELTA( solution.rGetSolutions()[num_timesteps-1][2], -25.0359, 1e-2 );
        TS_ASSERT_DELTA( solution.rGetSolutions()[num_timesteps-1][3],  77.2103, 1e-2 );
        TS_ASSERT_DELTA( solution.rGetSolutions()[num_timesteps-1][4],  20.6006, 1e-2 );
    }
    

    void TestKerchoffs2003ContractionModelConstantStrain() throw(Exception)
    {        
//todo: add tests, seconds to milliseconds, etc.

        Kerchoffs2003ContractionModel kerchoffs_model;
        TS_ASSERT_EQUALS(kerchoffs_model.IsStretchDependent(), true);
        TS_ASSERT_EQUALS(kerchoffs_model.IsStretchRateDependent(), false);

        EulerIvpOdeSolver euler_solver;
                
        ContractionModelInputParameters input_params;
        input_params.intracellularCalciumConcentration = DOUBLE_UNSET;

        TimeStepper stepper(0, 500, 1.0);  //ms
        std::cout << stepper.GetTime() << " " << kerchoffs_model.rGetStateVariables()[0] << " " << kerchoffs_model.GetActiveTension() <<  "\n";

        while(!stepper.IsTimeAtEnd())
        {
            input_params.time = stepper.GetTime();
            if( (stepper.GetTime()>100) && (stepper.GetTime()<600) )
            {
                input_params.voltage = 50;
            }
            else
            {
                input_params.voltage = -90;
            }
            kerchoffs_model.SetInputParameters(input_params);
                
            euler_solver.SolveAndUpdateStateVariable(&kerchoffs_model, stepper.GetTime(), stepper.GetNextTime(), 0.01);
            std::cout << stepper.GetTime() << " " << kerchoffs_model.rGetStateVariables()[0] << " " << kerchoffs_model.GetActiveTension() <<  "\n";
            
            stepper.AdvanceOneTimeStep();
        }
    }
    
    void dontTestNash2004ContractionModel() throw(Exception)
    {
        Nash2004ContractionModel nash_model;
        TS_ASSERT_EQUALS(nash_model.IsStretchDependent(), false);
        TS_ASSERT_EQUALS(nash_model.IsStretchRateDependent(), false);
        
        boost::shared_ptr<EulerIvpOdeSolver> p_euler_solver(new EulerIvpOdeSolver);

        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(-3.0,3.0));
        LuoRudyIModel1991OdeSystem lr91(p_euler_solver, p_stimulus);
        
        ContractionModelInputParameters input_params;
        input_params.intracellularCalciumConcentration = DOUBLE_UNSET;
        input_params.time = DOUBLE_UNSET;

        TimeStepper stepper(0, 1000, 1); 

        while(!stepper.IsTimeAtEnd())
        {
            lr91.Compute(stepper.GetTime(), stepper.GetNextTime());
            input_params.voltage = lr91.GetVoltage();
            nash_model.SetInputParameters(input_params);
            p_euler_solver->SolveAndUpdateStateVariable(&nash_model, stepper.GetTime(), stepper.GetNextTime(), 0.01);
            
            std::cout << stepper.GetTime() << " " << nash_model.GetActiveTension() <<  "\n";
            stepper.AdvanceOneTimeStep();            
        }
    }
};
#endif /*TESTCELLULARMECHANICSODESYSTEMS_HPP_*/
