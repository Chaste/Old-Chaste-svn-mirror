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

#ifndef TESTPEREGOCELLMODELS_HPP_
#define TESTPEREGOCELLMODELS_HPP_

#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "HeartConfig.hpp"
#include "PeregoLuoRudyIModel1991OdeSystem.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

#include "SimpleDataWriter.hpp"


class TestPeregoCellModels : public CxxTest::TestSuite
{
public:

    // This test checks that the EvaluatePredictedValues method computes the right thing
    // We compare against values from Chris' Matlab code
    //
    // This method does not call Compute - it has the time loop hardcoded
    void TestPeregoCellModelPredictor(void) throw(Exception)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 1.99; // ms
        double when = 0.0;

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver()); // Not actually used; just for compliance with the standard Chaste AbstractCardiacCell
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        // Create a luo Rudy cell set up for Perego-like solving
        PeregoLuoRudyIModel1991OdeSystem lr91_perego(p_solver, p_stimulus);


        std::vector<double> r_initial_values;

        // State variables at time 0
        r_initial_values.push_back(0.9804713); // h
        r_initial_values.push_back(0.98767124); // j
        r_initial_values.push_back(0.00187018); // m
        r_initial_values.push_back(0.0002); // Cai (mMol)
        r_initial_values.push_back(-83.853); // V (mV)
        r_initial_values.push_back(0.00316354); // d
        r_initial_values.push_back(0.99427859); // f
        r_initial_values.push_back(0.16647703); // X

        std::vector<double> r_previous_yvalues = r_initial_values;
        std::vector<double> r_predicted_values(r_previous_yvalues.size());
        r_predicted_values = r_initial_values;

        //as we don't have implemented all the functions yet, we make a time loop here
        for (unsigned i=0; i<1000; i++)
        {
            double time = i/100.0;
            // Initialise the previous state variables to be a copy of the current state variables
            std::vector<double> r_previous_yvalues = r_predicted_values;
            //predict the next value
            lr91_perego.EvaluatePredictedValues(r_previous_yvalues, r_predicted_values, time);
        }

        //values from Chris' code
        std::vector<double> matlab_answers;
        matlab_answers.push_back(1.6229e-27);
        matlab_answers.push_back(0.0776);
        matlab_answers.push_back(0.9987);
        matlab_answers.push_back(0.0013);
        matlab_answers.push_back(12.5641);
        matlab_answers.push_back(0.4017);
        matlab_answers.push_back(0.9764);
        matlab_answers.push_back(0.1919);

        //check that after 1000 time steps the results are the same.
        for(unsigned i=0; i<matlab_answers.size(); i++)
        {
            TS_ASSERT_DELTA(matlab_answers[i],r_predicted_values[i],1e-3);
        }

        //for coverage test I ionic against hardcoded value, possibly re-evaluate when tissue simulations are in place
        TS_ASSERT_DELTA(lr91_perego.GetIIonic(), 0.003,1e-4);
    }

    // This method does not call Compute - it has the time loop hardcoded
    void TestPeregoCellModelPredictorAndCorrector(void) throw(Exception)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 1.99; // ms
        double when = 0.0;

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver()); // Not actually used; just for compliance with the standard Chaste AbstractCardiacCell
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        // Create a luo Rudy cell set up for Perego-like solving
        PeregoLuoRudyIModel1991OdeSystem lr91_perego(p_solver, p_stimulus);


        std::vector<double> r_initial_values;

        // State variables at time 0
        r_initial_values.push_back(0.9804713); // h
        r_initial_values.push_back(0.98767124); // j
        r_initial_values.push_back(0.00187018); // m
        r_initial_values.push_back(0.0002); // Cai (mMol)
        r_initial_values.push_back(-83.853); // V (mV)
        r_initial_values.push_back(0.00316354); // d
        r_initial_values.push_back(0.99427859); // f
        r_initial_values.push_back(0.16647703); // X

        std::vector<double> r_previous_yvalues = r_initial_values;
        std::vector<double> r_predicted_values(r_previous_yvalues.size());
        std::vector<double> r_corrected_values(r_previous_yvalues.size());
        r_predicted_values = r_initial_values;
        r_corrected_values = r_initial_values;

        std::vector<double> times;
        std::vector<double> output[8];

        //as we don't have implemented all the functions yet, we make a time loop here
        for (unsigned i=0; i<1000; i++)
        {
            double time = i/100.0;
            // Initialise the previous state variables to be a copy of the current state variables
            std::vector<double> r_previous_yvalues = r_corrected_values;
            //predict the next value
            lr91_perego.EvaluatePredictedValues(r_previous_yvalues, r_predicted_values, time);
            lr91_perego.EvaluateCorrectedValues(r_predicted_values, r_corrected_values, time);
            times.push_back(time);
            for(unsigned i=0; i<8; i++)
            {
                output[i].push_back(r_corrected_values[i]);
            }
        }


        std::vector<std::vector<double> > data;
        data.push_back(times);
        for(unsigned i=0; i<8; i++)
        {
            data.push_back(output[i]);
        }
        SimpleDataWriter writer("TestPeregoLr91", "output.dat", data, false);


        //values from Chris' code
        std::vector<double> matlab_answers;
        matlab_answers.push_back(0.0);
        matlab_answers.push_back(0.0776);
        matlab_answers.push_back(0.9987);
        matlab_answers.push_back(0.0013);
        matlab_answers.push_back(12.5825);
        matlab_answers.push_back(0.4017);
        matlab_answers.push_back(0.9764);
        matlab_answers.push_back(0.1920);

        //check that after 1000 time steps the results are the same.
        for(unsigned i=0; i<matlab_answers.size(); i++)
        {
            TS_ASSERT_DELTA(matlab_answers[i],r_corrected_values[i],1e-3);
        }
    }


    void TestCompute(void) throw(Exception)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.00; //ms
        double when = 0.0;
        
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver()); // Not actually used; just for compliance with the standard Chaste AbstractCardiacCell
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        // Create a luo Rudy cell set up for Perego-like solving
        PeregoLuoRudyIModel1991OdeSystem lr91_perego(p_solver, p_stimulus, false);

        OdeSolution solutions = lr91_perego.Compute(0.0, 10.0);

        TS_ASSERT_EQUALS(solutions.rGetTimes()[0], 0.0);
        TS_ASSERT_DELTA(solutions.rGetTimes().back(), 10.0, 1e-9);
        TS_ASSERT_EQUALS(solutions.rGetTimes().size(), 1001u);

        //values from Chris' code
        std::vector<double> matlab_answers;
        matlab_answers.push_back(0.0);
        matlab_answers.push_back(0.0776);
        matlab_answers.push_back(0.9987);
        matlab_answers.push_back(0.0013);
        matlab_answers.push_back(12.5621);
        matlab_answers.push_back(0.4017);
        matlab_answers.push_back(0.9764);
        matlab_answers.push_back(0.1920);

        //check that after 1000 time steps the results are the same.
        for(unsigned i=0; i<matlab_answers.size(); i++)
        {
            TS_ASSERT_DELTA(matlab_answers[i],solutions.rGetSolutions().back()[i],1e-3);
        }
    }

    void TestCompareToStandardLuoRudy(void) throw(Exception)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 1.99; //ms
        double when = 0.0;

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver()); // Not actually used; just for compliance with the standard Chaste AbstractCardiacCell
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.001, 0.01, 0.01);

        // Compute the solution with the Perego nonadaptive predictor-corrector scheme
        // Create a luo Rudy cell set up for Perego-like solving
        PeregoLuoRudyIModel1991OdeSystem lr91_perego(p_solver, p_stimulus, false);
        OdeSolution solutions_perego = lr91_perego.Compute(0.0, 10.0);

        LuoRudyIModel1991OdeSystem lr91(p_solver, p_stimulus);
        OdeSolution solutions = lr91.Compute(0.0, 10.0);

        std::vector<double> tol(8);

        // Tolerance to error in each variable is based on a visual inspection of the difference
        // between the two AP traces and the maximum values of the differences between Perego and
        // forward Euler over time in this case.
        tol[0]=6e-2;
        tol[1]=4e-3;
        tol[2]=4e-2;
        tol[3]=1e-4;
        tol[4]=1;
        tol[5]=2e-3;
        tol[6]=3e-4;
        tol[7]=1e-3;

        TS_ASSERT_EQUALS(solutions_perego.rGetTimes().size(), solutions_perego.rGetSolutions().size());

        for(unsigned i=0; i<solutions_perego.rGetSolutions().size(); i++)
        {
            TS_ASSERT_DELTA(solutions.rGetTimes()[i],solutions_perego.rGetTimes()[i],1e-4);
            for(unsigned j=0; j<8; j++)
            {
                TS_ASSERT_DELTA(solutions.rGetSolutions()[i][j],solutions_perego.rGetSolutions()[i][j],tol[j]);

            }
        }
    }

    void TestEvaluateErrors(void) throw(Exception)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 1.99; //ms
        double when = 0.0;

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver()); // Not actually used; just for compliance with the standard Chaste AbstractCardiacCell
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        PeregoLuoRudyIModel1991OdeSystem lr91_perego(p_solver, p_stimulus, true);

        std::vector<double> errors(8); // Gets filled by method below
        std::vector<double> predicted_solution(8);
        std::vector<double> corrected_solution(8);

        std::vector<double> previous_corrected_solution(8);

        predicted_solution[0] = 0.980452408227077;
        predicted_solution[1] = 0.987668868231472;
        predicted_solution[2] = 0.001957111536055;
        predicted_solution[3] = 1.999029977193713e-4;
        predicted_solution[4] = -81.916308145968557;
        predicted_solution[5] = 0.003164107621352;
        predicted_solution[6] = 0.994286712717492;
        predicted_solution[7] = 0.166424185294100;

        corrected_solution[0] = 0.980408224875755;
        corrected_solution[1] = 0.987661246226870;
        corrected_solution[2] = 0.002189916542301;
        corrected_solution[3] = 1.99901905229788e-4;
        corrected_solution[4] = -81.929678049406846;
        corrected_solution[5] = 0.003165648122979;
        corrected_solution[6] = 0.994286670014978;
        corrected_solution[7] = 0.166424998680492;

        previous_corrected_solution[0] = 0.980470316683752;
        previous_corrected_solution[1] = 0.987671157351160;
        previous_corrected_solution[2] = 0.001896834011472;
        previous_corrected_solution[3] = 0.199901905229788e-3;
        previous_corrected_solution[4] = -83.598370438354308;
        previous_corrected_solution[5] = 0.003163567176509;
        previous_corrected_solution[6] = 0.994279658014422;
        previous_corrected_solution[7] = 0.166470067496488;
        

        double current_time = 0.02;

        // These need to be set correctly (usually by previous run)

        // Set the member variables to fake the predictor-corrector loop that would occur on a normal solve
        lr91_perego.mIsTheCorrectorStep=false;
        lr91_perego.mIsTheErrorEvaluationStep=false;
        //Compute parameters (done in the child class, will modify the member variables here)
        lr91_perego.ComputeSystemParameters(previous_corrected_solution, current_time);

        lr91_perego.mIsTheErrorEvaluationStep = true;
        //Compute parameters (done in the child class, will modify the member variables here)
        lr91_perego.ComputeSystemParameters(corrected_solution, current_time);
        lr91_perego.mIsTheErrorEvaluationStep = false;

        // Check that ComputeParameters works properly
        TS_ASSERT_DELTA(lr91_perego.ma_error[0], -0.18480, 1e-5);
        TS_ASSERT_DELTA(lr91_perego.ma_current[0], -0.23399, 1e-5);

        TS_ASSERT_EQUALS(lr91_perego.ma_error.size(), lr91_perego.ma_current.size());
        
        std::vector<double> matlab_ma_error(8);
        std::vector<double> matlab_ma_current(8);
        matlab_ma_error[0] = -0.18480;
        matlab_ma_error[1] = -0.0478996;
        matlab_ma_error[2] = -1.376890822e2;
        matlab_ma_error[3] = 0; //not a gating variable
        matlab_ma_error[4] = 0; //not a gating variable
        matlab_ma_error[5] = -0.11642;
        matlab_ma_error[6] = -0.01847;
        matlab_ma_error[7] = -0.00417;
        
        matlab_ma_current[0] = -0.23399;
        matlab_ma_current[1] = -0.05668;
        matlab_ma_current[2] = -1.60144412e2;
        matlab_ma_current[3] = 0; //not a gating variable
        matlab_ma_current[4] = 0; //not a gating variable
        matlab_ma_current[5] = -0.12097;
        matlab_ma_current[6] = -0.01872;
        matlab_ma_current[7] = -0.00433;       
        
        
        for(unsigned i=0; i<lr91_perego.ma_error.size();i++)
        {
            if(i!=3 && i!=4) // test is only run the gating variables; 3=calcium, 4=voltage
            {
                TS_ASSERT_DELTA(lr91_perego.ma_error[i], matlab_ma_error[i], 1e-5);
                TS_ASSERT_DELTA(lr91_perego.ma_current[i], matlab_ma_current[i], 1e-5);
            }
        }
        
        lr91_perego.mThetaP = -3.309930740729620;
        lr91_perego.mThetaC = -0.134988456788270;
        lr91_perego.mLocalTimeStep = 0.066198614814592;
        lr91_perego.EvaluateErrors(errors,predicted_solution,corrected_solution,current_time);

        std::vector<double> matlab_errors(8);
        matlab_errors[0] = 6.66198e-6;
        matlab_errors[1] = 1.13020e-6;
        matlab_errors[2] = 0.0050271;
        matlab_errors[3] = 1.61148e-10;
        matlab_errors[4] = 0.00197213;
        matlab_errors[5] = 2.29757e-7;
        matlab_errors[6] = 6.29993e-9;
        matlab_errors[7] = 1.19986e-7;

        TS_ASSERT_EQUALS(errors.size(), predicted_solution.size());


        for (unsigned i=0; i<errors.size(); i++)
        {
            TS_ASSERT_DELTA(errors[i],matlab_errors[i],(1e-4)*matlab_errors[i]);
        }

    }

    void TestAdaptiveTimesteps(void) throw(Exception)
    {
        std::cout << std::setprecision(8);
        
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 1.99; //ms
        double when = 0.0;

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver()); // Not actually used; just for compliance with the standard Chaste AbstractCardiacCell
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        // Compute the solution with the Perego nonadaptive predictor-corrector scheme
        // Create a luo Rudy cell set up for Perego-like solving
        PeregoLuoRudyIModel1991OdeSystem lr91_perego(p_solver, p_stimulus, true);
        OdeSolution solutions_perego = lr91_perego.Compute(0.0, 10.0);

        TS_ASSERT_EQUALS(solutions_perego.rGetTimes().size(), solutions_perego.rGetSolutions().size());

        std::vector<double> matlab_timesteps(11);

        matlab_timesteps[0] = 0.00;
        matlab_timesteps[1] = 0.01;
        matlab_timesteps[2] = 0.066198614814592;
        matlab_timesteps[3] = 0.066603712934993;
        matlab_timesteps[4] = 0.068085035666819;
        matlab_timesteps[5] = 0.069279542552099;
        matlab_timesteps[6] = 0.070507181748773;
        matlab_timesteps[7] = 0.071772602703771;
        matlab_timesteps[8] = 0.073092397564555;
        matlab_timesteps[9] = 0.074486709902735;
        matlab_timesteps[10] = 0.075980375254532;

        for(unsigned i=0; i<matlab_timesteps.size(); i++)
        {
            if (i>0)
            {   // Sum up the times
                matlab_timesteps[i] += matlab_timesteps[i-1];
            }
            TS_ASSERT_DELTA(matlab_timesteps[i],solutions_perego.rGetTimes()[i],1e-5);
        }
        TS_ASSERT_DELTA(solutions_perego.rGetTimes().back(),10.0,1e-9);
        
    }

    void TestFullActionPotential(void) throw(Exception)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 1.99; //ms
        double when = 0.0;

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver()); // Not actually used; just for compliance with the standard Chaste AbstractCardiacCell
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        // Compute the solution with the Perego nonadaptive predictor-corrector scheme
        // Create a luo Rudy cell set up for Perego-like solving
        PeregoLuoRudyIModel1991OdeSystem lr91_perego(p_solver, p_stimulus, true);
        lr91_perego.SetToleranceWeight(1e-2);
        OdeSolution solutions_perego = lr91_perego.Compute(0.0, 350.0);
        
        solutions_perego.WriteToFile("TestPeregoFullAp", "results","ms");
        
        unsigned final_time = solutions_perego.rGetTimes().size()-2u;
        //we check that the time step before the last one is the same as in the matlab code.
        TS_ASSERT_DELTA(solutions_perego.rGetTimes()[final_time],3.497529355*1e2,1e-6);
        
        PeregoLuoRudyIModel1991OdeSystem lr91_perego_2(p_solver, p_stimulus, true);

        lr91_perego_2.SetToleranceWeight(1e-3);
        OdeSolution solutions_perego_2 = lr91_perego_2.Compute(0.0, 350.0);
        final_time = solutions_perego_2.rGetTimes().size()-2u;
        
        TS_ASSERT_DELTA(solutions_perego_2.rGetTimes()[final_time],3.499228143*1e2,1e-6);
    }


};


#endif /*TESTPEREGOCELLMODELS_HPP_*/
