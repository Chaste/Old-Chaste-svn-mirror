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


#ifndef TESTCONTRACTIONMODELS_HPP_
#define TESTCONTRACTIONMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include "NonPhysiologicalContractionModel.hpp"
#include "NhsContractionModel.hpp"
#include "Kerchoffs2003ContractionModel.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "ZeroStimulus.hpp"
#include "TimeStepper.hpp"
#include "SimpleDataWriter.hpp"
#include "Nash2004ContractionModel.hpp"


// specify a functional form of lambda rather than get it from the mechanics.
// Use tanh so that lambda starts at 1.0 and decreases quickly to 0.8 halfway
// through the simulation
double MyLam(double t, double endTime, double minLam)
{
    double middle = (minLam + 1)/2;
    double scaled_t = 10*t/endTime - 5;
    return middle - (1-middle)*tanh(scaled_t);
}

double MyLamDeriv(double t, double endTime, double minLam)
{
    double middle = (minLam + 1)/2;
    double scaled_t = 10*t/endTime - 5;
    return -(10*(1-middle)/endTime)*(1-tanh(scaled_t)*tanh(scaled_t));
}

class TestContractionModels : public CxxTest::TestSuite
{
public :
    void TestNonPhysiologicalContractionModel() throw(Exception)
    {
        NonPhysiologicalContractionModel model1(1);
        NonPhysiologicalContractionModel model2(2);
        NonPhysiologicalContractionModel model3(3);

        ContractionModelInputParameters input_parameters;

        model1.SetInputParameters(input_parameters);
        model1.SetStretchAndStretchRate(1.0, 0.0);

        model2.SetInputParameters(input_parameters);
        model2.SetStretchAndStretchRate(1.0, 0.0);

        model3.SetInputParameters(input_parameters);
        model3.SetStretchAndStretchRate(1.0, 0.0);

        // call Run to set the end time = 0.0
        model1.RunDoNotUpdate(-0.1,0,0.1);
        model1.UpdateStateVariables(); // coverage
        model2.RunAndUpdate(-0.1,0,0.1); // does the same as RunAndUpdate
        model3.RunAndUpdate(-0.1,0,0.1);

        TS_ASSERT_DELTA(model1.GetActiveTension(), 0.0, 1e-12);
        TS_ASSERT_DELTA(model2.GetActiveTension(), 0.0, 1e-12);
        TS_ASSERT_DELTA(model3.GetActiveTension(), 0.0, 1e-12);

        model1.SetStretchAndStretchRate(0.8, 0.0);
        model2.SetStretchAndStretchRate(0.8, 0.0);
        model3.SetStretchAndStretchRate(0.8, 0.0);

        // call Run to set the end time = 1
        model1.RunDoNotUpdate(0,1,0.1);
        model2.RunDoNotUpdate(0,1,0.1);
        model3.RunDoNotUpdate(0,1,0.1);

        TS_ASSERT_DELTA(model1.GetActiveTension(), 5*sin(1.0), 1e-6);
        TS_ASSERT_DELTA(model2.GetActiveTension(), 4*sin(1.0), 1e-6);
        TS_ASSERT_DELTA(model3.GetActiveTension(), 5*exp(0.2)*sin(1.0), 1e-6);
    }


    void TestNhsContractionModelSimple() throw(Exception)
    {
        NhsContractionModel nhs_system;

        TS_ASSERT_THROWS_CONTAINS(nhs_system.GetNextActiveTension(), "If using this in an 'explicit manner'");

        TS_ASSERT_EQUALS(nhs_system.IsStretchDependent(), true);
        TS_ASSERT_EQUALS(nhs_system.IsStretchRateDependent(), true);

        // Hardcoded results for two values for z when lambda1=0.
        // Note: CalculateT0(z) is a private method.
        TS_ASSERT_DELTA(nhs_system.CalculateT0(0), 0, 1e-12);
        TS_ASSERT_DELTA(nhs_system.CalculateT0(1), 58.0648, 1e-3);

        TS_ASSERT_DELTA(nhs_system.GetCalciumTroponinValue(), 0.0, 0.01);

        boost::shared_ptr<EulerIvpOdeSolver> p_euler_solver(new EulerIvpOdeSolver);

        // the following is just to get a realistic Ca_I value
        boost::shared_ptr<ZeroStimulus> p_zero_stimulus(new ZeroStimulus);
        LuoRudyIModel1991OdeSystem lr91(p_euler_solver, p_zero_stimulus);
        unsigned Ca_i_index = lr91.GetStateVariableIndex("CaI");
        double Ca_I = lr91.rGetStateVariables()[Ca_i_index];

        // lambda1=1, dlamdt = 0, so there should be no active tension
        nhs_system.RunDoNotUpdate(0, 1, 0.01);
        nhs_system.UpdateStateVariables();

        TS_ASSERT_DELTA(nhs_system.GetActiveTension(), 0.0, 1e-12);

        // the following test doesn't make sense, as lambda=const, but dlam_dt > 0
        // but it is not possible to test the NHS system by itself without having a varying
        // lambda and get non-trivial solutions. So, we'll have a non-realistic
        // test here, and TS_ASSERT against hardcoded values, just to check nothing
        // has changed. A proper test where lambda varies (which means time-looping has
        // to be done outside the solver is done in TestElectroMechanicCellularModels,
        // where NHS is coupled to a cell model
        nhs_system.SetStretchAndStretchRate(0.5, 0.1);
        nhs_system.SetIntracellularCalciumConcentration(Ca_I);

        nhs_system.RunDoNotUpdate(0, 10, 0.01);
        nhs_system.UpdateStateVariables();

        TS_ASSERT_DELTA( nhs_system.rGetStateVariables()[0],   0.0056, 1e-2 );
        TS_ASSERT_DELTA( nhs_system.rGetStateVariables()[1],   0.0000, 1e-2 );
        TS_ASSERT_DELTA( nhs_system.rGetStateVariables()[2], -25.0574, 1e-2 );
        TS_ASSERT_DELTA( nhs_system.rGetStateVariables()[3],  77.2480, 1e-2 );
        TS_ASSERT_DELTA( nhs_system.rGetStateVariables()[4],  20.6006, 1e-2 );
    }



    /*
     *  A test which couples the NHS model with the Lr91 model, via CaI
     *
     *  There are two ways to force an active tension in the NHS model - changing
     *  the CaI concentration, or changing lam and dlam_dt. Here we keep lam=1, dlam_dt=0
     *  and stimulate an action potential in a Lr91 cell, in order to get
     *  active tension in the NHS model
     */
    void TestNhsConstantStretchVaryingCa() throw(Exception)
    {
        // setup
        double magnitude =  -25.5;
        double duration  =   2.0;  // ms
        double when      =   0.0;  // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        double end_time = 1000.0;

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        LuoRudyIModel1991OdeSystem electrophys_model(p_solver, p_stimulus);
        NhsContractionModel cellmech_model;

        // find out if electrophys model has CaTrop
        unsigned Ca_i_index = electrophys_model.GetStateVariableIndex("CaI");
        bool has_Ca_trop = false;
        unsigned Ca_trop_index=0;

        try
        {
            Ca_trop_index = electrophys_model.GetStateVariableIndex("CaTrop");
            has_Ca_trop = true;
        }
        catch(Exception& e)
        {
            has_Ca_trop = false;
        }

        std::vector<double> times;
        std::vector<double> active_tensions;
        std::vector<double> ca_trop;
        std::vector<double> z;
        std::vector<double> Ca_Is;
        //std::vector<double> Q1;   //Qi all zero as dlamdt=0
        //std::vector<double> Q2;
        //std::vector<double> Q3;


        // time loop
        for(double current_time = 0; current_time<end_time; current_time+=HeartConfig::Instance()->GetOdeTimeStep())
        {
            // solve electrophys model
            electrophys_model.Compute(current_time, current_time+HeartConfig::Instance()->GetOdeTimeStep());

            // get CaI
            double Ca_I = electrophys_model.rGetStateVariables()[Ca_i_index];
            cellmech_model.SetStretchAndStretchRate(1.0, 0.0);

            ContractionModelInputParameters input_parameters;
            input_parameters.intracellularCalciumConcentration = Ca_I;
            cellmech_model.SetInputParameters(input_parameters);

            // solve the cellular mechanics model
            cellmech_model.RunDoNotUpdate(current_time, current_time+HeartConfig::Instance()->GetOdeTimeStep(), HeartConfig::Instance()->GetOdeTimeStep());
            cellmech_model.UpdateStateVariables();


            // IF electrophy model has CaTrop, get CaTrop from
            // cellular mechanics models and send to electrophy model
            TS_ASSERT_EQUALS(has_Ca_trop, false);//NOTE - The following code is never activated....

            if(has_Ca_trop)
            {
                std::vector<double>& electrophys_model_state_vars = electrophys_model.rGetStateVariables();
                electrophys_model_state_vars[Ca_trop_index] = cellmech_model.GetCalciumTroponinValue();
            }

            times.push_back(current_time);
            active_tensions.push_back( cellmech_model.GetActiveTension() );
            ca_trop.push_back( cellmech_model.rGetStateVariables()[0] );
            z.push_back( cellmech_model.rGetStateVariables()[1] );
            Ca_Is.push_back( Ca_I );
            //Q1.push_back( cellmech_model.rGetStateVariables()[2] ); //Qi all zero as dlamdt=0
            //Q2.push_back( cellmech_model.rGetStateVariables()[3] );
            //Q3.push_back( cellmech_model.rGetStateVariables()[4] );
        }

        std::vector<std::vector<double> > data;
        data.push_back(times);
        data.push_back(active_tensions);
        data.push_back(ca_trop);
        data.push_back(z);
        data.push_back(Ca_Is);
        //data.push_back(Q1); //Qi all zero as dlamdt=0
        //data.push_back(Q2);
        //data.push_back(Q3);

        SimpleDataWriter writer("TestNhsContractionModel", "no_lam.dat", data);

        // having looked at the results to see if they look sensible...
        // The active tension in this case increases rapidly for 100ms, then plateaus for 300s,
        // then decreases slowly.
        //
        // 2 hardcoded tests as the answers look correct, on plateau and decreasing stage
        TS_ASSERT_DELTA(times[20000], 200, 1e-2);
        TS_ASSERT_DELTA(active_tensions[20000], 54.99, 1e-1);

        TS_ASSERT_DELTA(times[50000], 500, 1e-2);
        TS_ASSERT_DELTA(active_tensions[50000], 28.25, 1e-1);
    }


    /*
     *  A test which tests the NHS model with varying lambda
     *
     *  There are two ways to force an active tension in the NHS model - changing
     *  the CaI concentration, or changing lam and dlam_dt. Here keep CaI fixed (by
     *  solving a Lr91 model with no stimulus, and specify non-zero lambda. We specify
     *  a functional form of lambda rather than get it from a mechanics model
     */
    void TestNhsConstantCaVaryingStretch() throw(Exception)
    {
        boost::shared_ptr<ZeroStimulus> p_zero_stimulus(new ZeroStimulus);
        double end_time = 100.0;

        double min_lam = 0.85;

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        LuoRudyIModel1991OdeSystem electrophys_model(p_solver, p_zero_stimulus);
        NhsContractionModel cellmech_model;

        // find out if electrophys model has CaTrop
        unsigned Ca_i_index = electrophys_model.GetStateVariableIndex("CaI");
        bool has_Ca_trop = false;
        unsigned Ca_trop_index=0;

        try
        {
            Ca_trop_index = electrophys_model.GetStateVariableIndex("CaTrop");
            has_Ca_trop = true;
        }
        catch(Exception& e)
        {
            has_Ca_trop = false;
        }

        std::vector<double> times;
        std::vector<double> active_tensions;
        std::vector<double> ca_trop;
        std::vector<double> z;
        std::vector<double> Q1;
        std::vector<double> Q2;
        std::vector<double> Q3;

        // time loop
        for(double current_time = 0; current_time<end_time; current_time+=HeartConfig::Instance()->GetOdeTimeStep())
        {
            // solve electrophys model
            electrophys_model.Compute(current_time, current_time+HeartConfig::Instance()->GetOdeTimeStep());

            // get CaI
            double Ca_I = electrophys_model.rGetStateVariables()[Ca_i_index];

            // set lam, dlam_dt and Ca_i in the cellular mechanics model
            double lam = MyLam(current_time, end_time, min_lam);
            double dlam_dt = MyLamDeriv(current_time, end_time, min_lam);

            cellmech_model.SetStretchAndStretchRate(lam, dlam_dt);
            ContractionModelInputParameters input_parameters;
            input_parameters.intracellularCalciumConcentration = Ca_I;
            cellmech_model.SetInputParameters(input_parameters);

            // solve the cellular mechanics model
            cellmech_model.RunDoNotUpdate(current_time, current_time+HeartConfig::Instance()->GetOdeTimeStep(), HeartConfig::Instance()->GetOdeTimeStep());
            cellmech_model.UpdateStateVariables();


            // IF electrophy model has CaTrop, get CaTrop from
            // cellular mechanics models and send to electrophy model

            TS_ASSERT_EQUALS(has_Ca_trop, false);//NOTE - The following code is never activated....
            if(has_Ca_trop)
            {
                std::vector<double>& electrophys_model_state_vars = electrophys_model.rGetStateVariables();
                electrophys_model_state_vars[Ca_trop_index] = cellmech_model.GetCalciumTroponinValue();
            }

            times.push_back(current_time);
            active_tensions.push_back( cellmech_model.GetActiveTension() );
            ca_trop.push_back( cellmech_model.rGetStateVariables()[0] );
            z.push_back( cellmech_model.rGetStateVariables()[1] );
            Q1.push_back( cellmech_model.rGetStateVariables()[2] );
            Q2.push_back( cellmech_model.rGetStateVariables()[3] );
            Q3.push_back( cellmech_model.rGetStateVariables()[4] );
        }

        std::vector<std::vector<double> > data;
        data.push_back(times);
        data.push_back(active_tensions);
        data.push_back(ca_trop);
        data.push_back(z);
        data.push_back(Q1);
        data.push_back(Q2);
        data.push_back(Q3);

        SimpleDataWriter writer("TestNhsContractionModel", "specified_lambda.dat", data, false);

        // having looked at the results to see if they look sensible...
        // The active tension in this case increases rapidly for 100ms, then plateaus for 300s,
        // then decreases slowly.
        //
        // 2 hardcoded tests as the answers look correct, on plateau and decreasing stage
        TS_ASSERT_DELTA(times[3000], 30, 1e-2);
        TS_ASSERT_DELTA(active_tensions[3000],  0.0662, 1e-3);

        TS_ASSERT_DELTA(times[6000], 60, 1e-2);
        TS_ASSERT_DELTA(active_tensions[6000], -0.0031, 1e-3);
    }

    void TestRunDoNotUpdateEtcUsingKerchoff() throw(Exception)
    {
        Kerchoffs2003ContractionModel kerchoffs_model;
        kerchoffs_model.SetStretchAndStretchRate(0.85,0.0);

        ContractionModelInputParameters input_params;
        input_params.voltage = 50;
        kerchoffs_model.SetInputParameters(input_params);

        double initial_state_var = kerchoffs_model.rGetStateVariables()[0];
        double initial_active_tension = kerchoffs_model.GetActiveTension();

        TS_ASSERT_DELTA(initial_active_tension, 0.0, 1e-12);

        kerchoffs_model.RunDoNotUpdate(0,1,0.01);

        // check the state variable has not changed.
        TS_ASSERT_DELTA(initial_state_var, kerchoffs_model.rGetStateVariables()[0], 1e-10);

        // get the active tension for next time
        double new_active_tension = kerchoffs_model.GetNextActiveTension();

        // update
        kerchoffs_model.UpdateStateVariables();
        TS_ASSERT_DELTA(new_active_tension, kerchoffs_model.GetActiveTension(), 1e-10);

        double new_state_var = kerchoffs_model.rGetStateVariables()[0];

        // run using second model and euler solver
        Kerchoffs2003ContractionModel kerchoffs_model2;
        kerchoffs_model2.SetStretchAndStretchRate(0.85,0.0);
        kerchoffs_model2.SetInputParameters(input_params);
        EulerIvpOdeSolver euler_solver;

        euler_solver.SolveAndUpdateStateVariable(&kerchoffs_model2, 0, 1, 0.01);

        kerchoffs_model2.mTime = 1.0;

        TS_ASSERT_DELTA(new_state_var, kerchoffs_model2.rGetStateVariables()[0], 1e-10);
        TS_ASSERT_DELTA(new_active_tension, kerchoffs_model2.GetActiveTension(), 1e-10);

        // test RunAndUpdate
        Kerchoffs2003ContractionModel kerchoffs_model3;
        kerchoffs_model3.SetStretchAndStretchRate(0.85,0.0);
        kerchoffs_model3.SetInputParameters(input_params);
        kerchoffs_model3.RunAndUpdate(0, 1, 0.01);
        kerchoffs_model3.mTime = 1.0;
        TS_ASSERT_DELTA(kerchoffs_model2.rGetStateVariables()[0], kerchoffs_model3.rGetStateVariables()[0], 1e-10);
        TS_ASSERT_DELTA(kerchoffs_model2.GetActiveTension(),      kerchoffs_model3.GetActiveTension(), 1e-10);
    }

    void TestKerchoffs2003ContractionModelConstantStretch() throw(Exception)
    {
        Kerchoffs2003ContractionModel kerchoffs_model;
        TS_ASSERT_EQUALS(kerchoffs_model.IsStretchDependent(), true);
        TS_ASSERT_EQUALS(kerchoffs_model.IsStretchRateDependent(), false);

        EulerIvpOdeSolver euler_solver;

        kerchoffs_model.SetStretchAndStretchRate(0.85,0.0);

        ContractionModelInputParameters input_params;
        input_params.intracellularCalciumConcentration = DOUBLE_UNSET;

        std::vector<double> times;
        std::vector<double> active_tensions;

        TimeStepper stepper(0, 1000, 1.0);  //ms
        times.push_back(stepper.GetTime());
        active_tensions.push_back(kerchoffs_model.GetActiveTension());

        while(!stepper.IsTimeAtEnd())
        {
            // specify a step-change voltage since this model gets activated at V=0 and deactivated at V=-70
            if( (stepper.GetTime()>100) && (stepper.GetTime()<600) )
            {
                input_params.voltage = 50;
            }
            else
            {
                input_params.voltage = -90;
            }
            kerchoffs_model.SetInputParameters(input_params);

            kerchoffs_model.RunDoNotUpdate(stepper.GetTime(), stepper.GetNextTime(), 0.01);
            kerchoffs_model.UpdateStateVariables();

            times.push_back(stepper.GetTime());
            active_tensions.push_back(kerchoffs_model.GetActiveTension());

            stepper.AdvanceOneTimeStep();
        }

        std::vector<std::vector<double> > data;
        data.push_back(times);
        data.push_back(active_tensions);
        SimpleDataWriter writer("TestKerchoffContractionModel", "constant_lam.dat", data);

        // visualise to verify validity..

        // hardcoded test, somewhere near the peak
        TS_ASSERT_DELTA(times[274], 273, 1e-2);
        TS_ASSERT_DELTA(active_tensions[274], 2.1633, 1e-2);

        TS_ASSERT_DELTA(active_tensions.back(), 0.0,  1e-2);
    }

    void TestKerchoffs2003ContractionModelVaryingStetch() throw(Exception)
    {
        Kerchoffs2003ContractionModel kerchoffs_model;
        TS_ASSERT_EQUALS(kerchoffs_model.IsStretchDependent(), true);
        TS_ASSERT_EQUALS(kerchoffs_model.IsStretchRateDependent(), false);

        EulerIvpOdeSolver euler_solver;

        kerchoffs_model.SetStretchAndStretchRate(0.85,0.0);

        ContractionModelInputParameters input_params;
        input_params.intracellularCalciumConcentration = DOUBLE_UNSET;

        std::vector<double> times;
        std::vector<double> stretches;
        std::vector<double> active_tensions;

        TimeStepper stepper(0, 600, 1.0);  //ms

        while(!stepper.IsTimeAtEnd())
        {
            // specify a step-change voltage since this model gets activated at V=0 and deactivated at V=-70
            if( (stepper.GetTime()>0) && (stepper.GetTime()<500) )
            {
                input_params.voltage = 50;
            }
            else
            {
                input_params.voltage = -90;
            }
            kerchoffs_model.SetInputParameters(input_params);

            // linearly decrease from 1 to 0.85 in first 100ms then increase up to 1 at 500ms
            double stretch = stepper.GetTime()<200 ? 1.0 - (stepper.GetTime()*0.15)/200 : 0.85 + ((stepper.GetTime()-200)*0.15)/300;

            kerchoffs_model.SetStretch(stretch);

            kerchoffs_model.RunDoNotUpdate(stepper.GetTime(), stepper.GetNextTime(), 0.01);
            kerchoffs_model.UpdateStateVariables();

            times.push_back(stepper.GetTime());
            stretches.push_back(stretch);
            active_tensions.push_back(kerchoffs_model.GetActiveTension());

            stepper.AdvanceOneTimeStep();
        }

        std::vector<std::vector<double> > data;
        data.push_back(times);
        data.push_back(stretches);
        data.push_back(active_tensions);
        SimpleDataWriter writer("TestKerchoffContractionModel", "linear_lam.dat", data, false);

        // visualise to verify validity.. // EMTODO2

        // hardcoded test, somewhere near the two peaks
        TS_ASSERT_DELTA(times[100], 100, 1e-2);
        TS_ASSERT_DELTA(active_tensions[100], 17.0696, 1e-2);
        TS_ASSERT_DELTA(times[250], 250, 1e-2);
        TS_ASSERT_DELTA(active_tensions[250], 3.3597,  1e-2);

        TS_ASSERT_DELTA(active_tensions.back(), 0.0,  1e-2);
    }


    void TestNash2004ContractionLaw() throw(Exception)
    {
        Nash2004ContractionModel nash_model;
        TS_ASSERT_EQUALS(nash_model.IsStretchDependent(), false);
        TS_ASSERT_EQUALS(nash_model.IsStretchRateDependent(), false);

        // coverage, doesn't use these..
        nash_model.SetStretchAndStretchRate(0.99,0.0);

        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(-25.5, 2.0, 0.0));

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        LuoRudyIModel1991OdeSystem electrophys_model(p_solver, p_stimulus);

        std::vector<double> times;
        std::vector<double> active_tensions;
        std::vector<double> voltages;

        double end_time = 1000.0;
        double dt = 0.01;

        // time loop
        for(double time = 0; time<end_time; time+=dt)
        {
            // solve electrophys model
            electrophys_model.Compute(time, time+dt);

            ContractionModelInputParameters input_parameters;
            input_parameters.voltage = electrophys_model.GetVoltage();
            nash_model.SetInputParameters(input_parameters);

            nash_model.RunAndUpdate(time, time+dt, dt);

            times.push_back(time);
            voltages.push_back(electrophys_model.GetVoltage());
            active_tensions.push_back(nash_model.GetActiveTension());
        }

        std::vector<std::vector<double> > data;
        data.push_back(times);
        data.push_back(active_tensions);
        data.push_back(voltages);
        SimpleDataWriter writer("TestNash2004ContractionModel", "ta.dat", data);

        //visualise results

        // somewhere around the peak
        TS_ASSERT_DELTA(times[5000],50,1e-3);
        TS_ASSERT_DELTA(active_tensions[5000],35.150,1e-2);

        // towards the end
        TS_ASSERT_DELTA(times[100000],1000,1e-3);
        TS_ASSERT_DELTA(active_tensions[100000],0.2496,1e-2);
    }
};
#endif /*TESTCONTRACTIONMODELS_HPP_*/
