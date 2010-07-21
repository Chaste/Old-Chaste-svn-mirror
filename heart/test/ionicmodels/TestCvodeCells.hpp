/*

Copyright (C) University of Oxford, 2005-2010

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

#ifndef TESTCVODECELLS_HPP_
#define TESTCVODECELLS_HPP_

#include "LuoRudy1991Cvode.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "Shannon2004.hpp"
#include "Shannon2004Cvode.hpp"

#include "HeartConfig.hpp"
#include "RegularStimulus.hpp"
#include "SimpleStimulus.hpp"
#include "ZeroStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RunAndCheckIonicModels.hpp"
#include "OdeSystemInformation.hpp"

#ifdef CHASTE_CVODE

class ExceptionalCell : public AbstractCvodeCell
{
private:
    bool mNice;
public :
    ExceptionalCell(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
                    boost::shared_ptr<AbstractStimulusFunction> pStimulus)
        : AbstractCvodeCell(pOdeSolver, 2, 0, pStimulus)
    {
        mNice = false;
        mpSystemInfo = OdeSystemInformation<ExceptionalCell>::Instance();
        Init();
        // For coverage purposes, call Init twice
        Init();
    }

    void BeNice()
    {
        mNice = true;
    }

    void EvaluateRhs(double time, N_Vector y, N_Vector ydot)
    {
        NV_Ith_S(ydot, 0) = NV_Ith_S(y, 1);
        NV_Ith_S(ydot, 1) = -NV_Ith_S(y, 0);
        if (!mNice)
        {
            EXCEPTION(DumpState("I'm feeling nasty!"));
        }
    }
};

template<>
void OdeSystemInformation<ExceptionalCell>::Initialise(void)
{
    // State variables
    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("mV");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("m");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}
#endif // CHASTE_CVODE

class TestCvodeCells : public CxxTest::TestSuite
{
public:

    void TestLuoRudyCvodeCell() throw(Exception)
    {
#ifdef CHASTE_CVODE
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double start = 50.0; // ms
        double period = 500; // ms
        boost::shared_ptr<ZeroStimulus> p_zero_stimulus(new ZeroStimulus()); // For coverage of SetStimulusFunction()
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude,
                                                                          duration,
                                                                          period,
                                                                          start));
        double start_time = 0.0;
        double end_time = 1000.0; //One second in milliseconds
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        // Make a model that uses Cvode directly:
        CellLuoRudy1991FromCellMLCvode lr91_cvode_system(p_solver, p_zero_stimulus);

        // Cover swapping to a proper stimulus
        lr91_cvode_system.SetStimulusFunction(p_stimulus);

        boost::shared_ptr<AbstractStimulusFunction> p_abs_stim = lr91_cvode_system.GetStimulusFunction();
        double period_back = boost::static_pointer_cast<RegularStimulus>(p_abs_stim)->GetPeriod();
        TS_ASSERT_DELTA(period_back,period,1e-7);
        TS_ASSERT_EQUALS(p_abs_stim,p_stimulus);

        TS_ASSERT_EQUALS(lr91_cvode_system.GetVoltageIndex(), 0u);
        TS_ASSERT_EQUALS(lr91_cvode_system.GetMaxSteps(), 0u); // 0 means 'UNSET' and Cvode uses the default.

        // 'Traditional' Chaste cell model for comparison of results:
        LuoRudyIModel1991OdeSystem lr91_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        double max_timestep = 1.0;
        double sampling_time = 1.0;

        OdeSolution solution_cvode = lr91_cvode_system.Solve(start_time, end_time, max_timestep, sampling_time);
        OdeSolution solution_chaste = lr91_ode_system.Compute(start_time, end_time);

        TS_ASSERT_DELTA(lr91_ode_system.GetIIonic(), lr91_cvode_system.GetIIonic(), 1e-3);

        unsigned step_per_row_chaste = 100u;
        bool clean_dir = false;
        solution_cvode.WriteToFile("TestCvodeCells","lr91_cvode","ms",1,clean_dir);
        solution_chaste.WriteToFile("TestCvodeCells","lr91_chaste","ms",step_per_row_chaste,clean_dir);

        double tolerance = 2e-1;
        bool voltage_only = false;
        CompareCellModelResults("lr91_cvode", "lr91_chaste",tolerance, voltage_only, "TestCvodeCells");

        // Clamping
        lr91_cvode_system.SetVoltageDerivativeToZero();
        solution_cvode = lr91_cvode_system.Solve(end_time, end_time+100.0, max_timestep, sampling_time);
        std::vector<double> voltages = solution_cvode.GetVariableAtIndex(lr91_cvode_system.GetVoltageIndex());
        for (unsigned i=0; i<voltages.size(); i++)
        {
            TS_ASSERT_EQUALS(voltages[i], lr91_cvode_system.GetVoltage());
        }

        lr91_cvode_system.SetVoltageDerivativeToZero(false);

        // Reset CVODE cell to initial conditions, and solve without sampling
        lr91_cvode_system.SetStateVariables(lr91_cvode_system.GetInitialConditions());
        lr91_cvode_system.SetMaxSteps(10000);
        TS_ASSERT_EQUALS(lr91_cvode_system.GetMaxSteps(), 10000u);
        lr91_cvode_system.Solve(start_time, end_time, max_timestep);
        TS_ASSERT_DELTA(lr91_cvode_system.GetVoltage(), lr91_ode_system.GetVoltage(), 1e-3);
        
        // Test parameter
        TS_ASSERT_EQUALS(lr91_cvode_system.GetNumberOfParameters(), 1u);
        TS_ASSERT_EQUALS(lr91_cvode_system.rGetParameterNames()[0], "fast_sodium_current_conductance");
        TS_ASSERT_EQUALS(lr91_cvode_system.rGetParameterUnits()[0], "milliS_per_cm2");
        TS_ASSERT_EQUALS(lr91_cvode_system.GetParameterIndex("fast_sodium_current_conductance"), 0u);
        TS_ASSERT_EQUALS(lr91_cvode_system.GetParameterUnits(0u), "milliS_per_cm2");
        TS_ASSERT_EQUALS(lr91_cvode_system.GetParameter(0u), 23.0);
        lr91_cvode_system.SetParameter(0u, 10.0);
        TS_ASSERT_EQUALS(lr91_cvode_system.GetParameter(0u), 10.0);
        lr91_cvode_system.SetParameter(0u, 23.0);

        // Parameter exceptions
        TS_ASSERT_THROWS_THIS(lr91_cvode_system.SetParameter(1u, -1), "The index passed in must be less than the number of parameters.");
        TS_ASSERT_THROWS_THIS(lr91_cvode_system.GetParameterIndex("b"), "No parameter named 'b'.");
        TS_ASSERT_THROWS_THIS(lr91_cvode_system.GetParameter(1u), "The index passed in must be less than the number of parameters.");
        TS_ASSERT_THROWS_THIS(lr91_cvode_system.GetParameterUnits(1u), "The index passed in must be less than the number of parameters.");

        
        // Coverage
        boost::shared_ptr<const AbstractOdeSystemInformation> p_sys_info = lr91_cvode_system.GetSystemInformation();
        TS_ASSERT(p_sys_info->rGetStateVariableNames() == lr91_cvode_system.rGetStateVariableNames());
        TS_ASSERT(p_sys_info->rGetStateVariableUnits() == lr91_cvode_system.rGetStateVariableUnits());
        TS_ASSERT_EQUALS(lr91_cvode_system.GetStateVariableIndex("membrane_voltage"),
                         lr91_cvode_system.GetVoltageIndex());
        TS_ASSERT_EQUALS(lr91_cvode_system.GetStateVariable(0),lr91_cvode_system.GetVoltage());
        TS_ASSERT_EQUALS(lr91_cvode_system.GetStateVariableUnits(0), "millivolt");

        TS_ASSERT_DELTA(lr91_cvode_system.GetRelativeTolerance(), 1e-4, 1e-10);
        TS_ASSERT_DELTA(lr91_cvode_system.GetAbsoluteTolerance(), 1e-6, 1e-10);
        TS_ASSERT_DELTA(lr91_cvode_system.GetLastStepSize(), 1.0, 1e-4);
        double old_v = lr91_cvode_system.GetVoltage();
        const double new_v = -1000.0;
        lr91_cvode_system.SetVoltage(new_v);
        TS_ASSERT_DELTA(lr91_cvode_system.GetVoltage(), new_v, 1e-12);
        lr91_cvode_system.SetVoltage(old_v);

        // Cover errors
        std::cout << "Testing Error Handling...\n";
        lr91_cvode_system.SetMaxSteps(2);
        TS_ASSERT_THROWS_CONTAINS(OdeSolution solution_cvode_fail = lr91_cvode_system.Solve(start_time, end_time, max_timestep, sampling_time),
                                  "mxstep steps taken before reaching tout");
        // Kill the cell
        boost::shared_ptr<SimpleStimulus> p_boom_stimulus(new SimpleStimulus(-50000, 2.0, 1.0));
        CellLuoRudy1991FromCellMLCvode lr91_boom(p_solver, p_boom_stimulus);
        TS_ASSERT_THROWS_CONTAINS(OdeSolution solution_boom = lr91_boom.Solve(start_time, end_time, max_timestep, sampling_time),
                                  " failed repeatedly or with |h| = hmin.");
        lr91_boom.SetStateVariables(lr91_boom.GetInitialConditions());
        lr91_boom.SetMaxSteps(10000);
        TS_ASSERT_THROWS_CONTAINS(lr91_boom.Solve(start_time, end_time, max_timestep),
                                  " failed repeatedly or with |h| = hmin.");
        // Nasty cell
        ExceptionalCell bad_cell(p_solver, p_boom_stimulus);
        TS_ASSERT_THROWS_THIS(OdeSolution solution_bad = bad_cell.Solve(start_time, end_time, max_timestep, sampling_time),
                              "CVODE Error -8 in module CVODE function CVode: At t = 0, the right-hand side routine failed in an unrecoverable manner.");
        bad_cell.SetStateVariables(bad_cell.GetInitialConditions());
        TS_ASSERT_THROWS_THIS(bad_cell.Solve(start_time, end_time, max_timestep),
                              "CVODE Error -8 in module CVODE function CVode: At t = 0, the right-hand side routine failed in an unrecoverable manner.");

        TS_ASSERT_THROWS_THIS(lr91_cvode_system.UseCellMLDefaultStimulus(),"This class has no default stimulus from CellML metadata.");
#endif // CHASTE_CVODE
    }

    void TestShannon2004() throw(Exception)
    {
#ifdef CHASTE_CVODE
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double start = 50.0; // ms
        double period = 500; // ms
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude,
                                                                          duration,
                                                                          period,
                                                                          start));

        double ode_time_step = 0.001;
        HeartConfig::Instance()->SetOdeTimeStep(ode_time_step);
        double start_time = 0.0;
        double end_time = 1000.0; //One second in milliseconds
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        // Make a model that uses Cvode directly:
        CellShannon2004FromCellMLCvode sh04_cvode_system(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(sh04_cvode_system.GetVoltageIndex(), 0u);
        TS_ASSERT_EQUALS(sh04_cvode_system.GetMaxSteps(), 0u); // 0 means 'UNSET' and Cvode uses the default.

        // 'Traditional' Chaste cell model for comparison of results:
        CellShannon2004FromCellML sh04_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        double max_timestep = 1.0;
        double sampling_time = 1.0;

        // Take a copy of the state variables, solve once, reset state variables and solve again (should get same answer)
        N_Vector state_vars = sh04_cvode_system.GetStateVariables(); // This returns the initial conditions
        sh04_cvode_system.Solve(start_time, end_time, max_timestep, sampling_time);
        N_Vector state_vars_after_solve = sh04_cvode_system.GetStateVariables(); // This returns the state variables at the end (for coverage)...

        sh04_cvode_system.SetStateVariablesUsingACopyOfThisVector(state_vars);
        OdeSolution solution_cvode = sh04_cvode_system.Solve(start_time, end_time, max_timestep, sampling_time);

        // Check we get the right answer both times
        std::vector<double> final_state_vars = solution_cvode.rGetSolutions().back();
        for (unsigned i=0; i<final_state_vars.size(); i++)
        {
            TS_ASSERT_DELTA(final_state_vars[i], NV_Ith_S(state_vars_after_solve, i), 1e-6);
        }

        // Solve using Chaste solvers for comparison.
        OdeSolution solution_chaste = sh04_ode_system.Compute(start_time, end_time);

        unsigned step_per_row_chaste = 1000u;
        bool clean_dir = false;
        solution_cvode.WriteToFile("TestCvodeCells","sh04_cvode","ms",1,clean_dir);
        solution_chaste.WriteToFile("TestCvodeCells","sh04_chaste","ms",step_per_row_chaste,clean_dir);

        double tolerance = 3e-2;
        bool voltage_only = false;
        CompareCellModelResults("sh04_cvode", "sh04_chaste",tolerance, voltage_only, "TestCvodeCells");

        // Coverage of GetIIonic method.
        TS_ASSERT_DELTA(sh04_ode_system.GetIIonic(), sh04_cvode_system.GetIIonic(), 1e-4);
        TS_ASSERT_DELTA(sh04_cvode_system.GetIIonic(), 0.0006, 1e-4);

        // Clamping
        sh04_cvode_system.SetVoltageDerivativeToZero();
        solution_cvode = sh04_cvode_system.Solve(end_time, end_time+100.0, max_timestep, sampling_time);
        std::vector<double> voltages = solution_cvode.GetVariableAtIndex(sh04_cvode_system.GetVoltageIndex());
        for (unsigned i=0; i<voltages.size(); i++)
        {
            TS_ASSERT_EQUALS(voltages[i], sh04_cvode_system.GetVoltage());
        }

        sh04_cvode_system.SetVoltageDerivativeToZero(false);

        // Coverage of mSetVoltageDerivativeToZero in non-CVODE class
        sh04_ode_system.ComputeExceptVoltage(end_time,end_time+0.01);

        // Tidy up
        state_vars->ops->nvdestroy(state_vars);
        state_vars_after_solve->ops->nvdestroy(state_vars_after_solve);
#endif // CHASTE_CVODE
    }
};


#endif /*TESTCVODECELLS_HPP_*/
