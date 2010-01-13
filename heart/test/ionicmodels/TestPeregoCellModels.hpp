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

class TestPeregoCellModels : public CxxTest::TestSuite
{
public:

    //This tst checks that the EvaluatePredictedGates method computes the right thing
    // We compare against values from Chris' Matlab code
    void TestPeregoCellModelPredictor(void) throw(Exception)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 200 ;//for this test it will be number of time steps
        double when = 0.0; 
        
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        // Create a luo Rudy cell set up for Perego-like solving
        PeregoLuoRudyIModel1991OdeSystem lr91_perego( p_stimulus);

        
        std::vector<double> rInitialValues;
             
        // State variables at time 0
        rInitialValues.push_back(0.9804713); // h
        rInitialValues.push_back(0.98767124); // j
        rInitialValues.push_back(0.00187018); // m
        rInitialValues.push_back(0.0002); // Cai (mMol)
        rInitialValues.push_back(-83.853); // V (mV)
        rInitialValues.push_back(0.00316354); // d
        rInitialValues.push_back(0.99427859); // f
        rInitialValues.push_back(0.16647703); // X

        std::vector<double> rPreviousYValues = rInitialValues;
        std::vector<double> rPredictedValues(rPreviousYValues.size());
        rPredictedValues = rInitialValues;
        
        //as we don't have implemented all the functions yet, we make a time loop here
        for (unsigned time_step=0;time_step<1000;time_step++)
        {        
            // Initialise the previous state variables to be a copy of the current state variables
            std::vector<double> rPreviousYValues = rPredictedValues;
            //predict the next value
            lr91_perego.EvaluatePredictedGates(rPreviousYValues, rPredictedValues, time_step);
        }   
        
        //values from Chris' code
        std::vector<double> MatlabAnswers;      
        MatlabAnswers.push_back(1.7e-27);
        MatlabAnswers.push_back(0.0778);
        MatlabAnswers.push_back(0.9987);
        MatlabAnswers.push_back(0.0013);
        MatlabAnswers.push_back(12.71);
        MatlabAnswers.push_back(0.4015);
        MatlabAnswers.push_back(0.9765);
        MatlabAnswers.push_back(0.1921);
        
        //check that after 1000 time steps the results are the same.
        for(unsigned i=0; i<MatlabAnswers.size(); i++)
        {
            TS_ASSERT_DELTA(MatlabAnswers[i],rPredictedValues[i],1e-3);
        }
        
        //for coverage test I ionic against hardcoded value
        ///\TODO make a more meaningful test for this 
        TS_ASSERT_DELTA(lr91_perego.GetIIonic(), 0.003,1e-4);
    }
    
};


#endif /*TESTPEREGOCELLMODELS_HPP_*/
