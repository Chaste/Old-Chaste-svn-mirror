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


#ifndef TESTELECTROMECHANICCALELLULARMODELS_HPP_
#define TESTELECTROMECHANICCALELLULARMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include "NhsCellularMechanicsOdeSystem.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "SimpleDataWriter.hpp"
#include "ZeroStimulus.hpp"


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


class TestElectroMechanicCellularModels : public CxxTest::TestSuite
{
public:
    /*
     *  A test which couples the NHS model with the Lr91 model, via CaI
     *
     *  There are two ways to force an active tension in the NHS model - changing
     *  the CaI concentration, or changing lam and dlam_dt. Here we keep lam=1, dlam_dt=0
     *  and stimulate an action potential in a Lr91 cell, in order to get
     *  active tension in the NHS model
     */
    void TestSingleCellWithForcedByCaI() throw(Exception)
    {
        // setup
        double magnitude =  -25.5;
        double duration  =   2.0;  // ms
        double when      =   0.0;  // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        double end_time = 1000.0;

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        LuoRudyIModel1991OdeSystem electrophys_model(p_solver, p_stimulus);
        NhsCellularMechanicsOdeSystem cellmech_model;

        // find out if electrophys model has CaTrop
        unsigned Ca_i_index = electrophys_model.GetStateVariableNumberByName("CaI");
        bool has_Ca_trop = false;
        unsigned Ca_trop_index=0;

        try
        {
            Ca_trop_index = electrophys_model.GetStateVariableNumberByName("CaTrop");
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
            input_parameters.IntracellularCalciumConcentrations = Ca_I;    
            cellmech_model.SetInputParameters(input_parameters);

            // solve the cellular mechanics model
            p_solver->SolveAndUpdateStateVariable(&cellmech_model, current_time, current_time+HeartConfig::Instance()->GetOdeTimeStep(), HeartConfig::Instance()->GetOdeTimeStep());

            // IF electrophy model has CaTrop, get CaTrop from
            // cellular mechanics models and send to electrophy model
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

        SimpleDataWriter writer("TestElectroMechanicCellularModels", "no_lam.dat", data);

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
    void TestSingleCellWithSpecifiedLambda() throw(Exception)
    {
        boost::shared_ptr<ZeroStimulus> p_zero_stimulus(new ZeroStimulus);
        double end_time = 100.0;

        double min_lam = 0.85;

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        LuoRudyIModel1991OdeSystem electrophys_model(p_solver, p_zero_stimulus);
        NhsCellularMechanicsOdeSystem cellmech_model;

        // find out if electrophys model has CaTrop
        unsigned Ca_i_index = electrophys_model.GetStateVariableNumberByName("CaI");
        bool has_Ca_trop = false;
        unsigned Ca_trop_index=0;

        try
        {
            Ca_trop_index = electrophys_model.GetStateVariableNumberByName("CaTrop");
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
            input_parameters.IntracellularCalciumConcentrations = Ca_I;    
            cellmech_model.SetInputParameters(input_parameters);

            // solve the cellular mechanics model
            p_solver->SolveAndUpdateStateVariable(&cellmech_model, current_time, current_time+HeartConfig::Instance()->GetOdeTimeStep(), HeartConfig::Instance()->GetOdeTimeStep());

            // IF electrophy model has CaTrop, get CaTrop from
            // cellular mechanics models and send to electrophy model
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

        SimpleDataWriter writer("TestElectroMechanicCellularModels", "specified_lambda.dat", data, false);

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
};
#endif /*TESTELECTROMECHANICCALELLULARMODELS_HPP_*/
