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
#include "Debug.hpp"


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
        
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        // Create a luo Rudy cell set up for Perego-like solving
        PeregoLuoRudyIModel1991OdeSystem lr91_perego( p_stimulus);

        
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
        
        //for coverage test I ionic against hardcoded value
        TS_ASSERT_DELTA(lr91_perego.GetIIonic(), 0.003,1e-4);
    }
    
    // This method does not call Compute - it has the time loop hardcoded
    void TestPeregoCellModelPredictorAndCorrector(void) throw(Exception)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 1.99; // ms
        double when = 0.0; 
        
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        // Create a luo Rudy cell set up for Perego-like solving
        PeregoLuoRudyIModel1991OdeSystem lr91_perego( p_stimulus);

        
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
        double duration  = 1.99;//for this test it will be number of time steps
        double when = 0.0; 
        
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        // Create a luo Rudy cell set up for Perego-like solving
        PeregoLuoRudyIModel1991OdeSystem lr91_perego( p_stimulus);
        
        OdeSolution solutions = lr91_perego.Compute(0.0, 10.0);
        
        TS_ASSERT_EQUALS(solutions.rGetTimes()[0], 0.0);
        TS_ASSERT_DELTA(solutions.rGetTimes().back(), 10.0, 1e-12);
        TS_ASSERT_EQUALS(solutions.rGetTimes().size(), 1001u);

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
            TS_ASSERT_DELTA(matlab_answers[i],solutions.rGetSolutions().back()[i],1e-3);
        }
    }
    
    void TestCompareToStandardLuoRudy(void) throw(Exception)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 1.99;//for this test it will be number of time steps
        double when = 0.0; 
       
        
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        // Compute the solution with the Perego nonadaptive predictor-corrector scheme
        // Create a luo Rudy cell set up for Perego-like solving
        PeregoLuoRudyIModel1991OdeSystem lr91_perego( p_stimulus);
        OdeSolution solutions_perego = lr91_perego.Compute(0.0, 10.0);
        
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        LuoRudyIModel1991OdeSystem lr91(p_solver, p_stimulus);
        OdeSolution solutions = lr91.Compute(0.0, 10.0);
        
        std::vector<double> tol(8);
        
        
        // Tolerance to error in each variable is based on a visual inspection of the difference
        // between the two AP traces and the maximum values of the differences between Perego and 
        // forward Euler over time in this case.
        tol[0]=6e-2;
        tol[1]=4e-3;
        tol[2]=4e-2;
        
        tol[3]=1e-5;
        tol[4]=6;
        
        tol[5]=2e-3;
        tol[6]=3e-4;
        tol[7]=4e-4;
        
        for(unsigned i=0; i<solutions.rGetTimes().size(); i++)
        {
            TS_ASSERT_DELTA(solutions.rGetTimes()[i],solutions_perego.rGetTimes()[i],1e-12);
            for(unsigned j=0; j<8; j++)
            {
                TS_ASSERT_DELTA(solutions.rGetSolutions()[i][j],solutions_perego.rGetSolutions()[i][j],tol[j]);
            }
        }
        
        solutions.WriteToFile("TestPeregoLr91Compare","lr91_standard",&lr91,"ms");
        solutions_perego.WriteToFile("TestPeregoLr91Compare","lr91_perego",&lr91_perego,"ms", 1, false);
    }
};


#endif /*TESTPEREGOCELLMODELS_HPP_*/
