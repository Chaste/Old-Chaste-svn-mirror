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
#ifndef TESTODEBASEDCELLCYCLEMODELS_HPP_
#define TESTODEBASEDCELLCYCLEMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "TysonNovakCellCycleModel.hpp"
#include "IngeWntSwatCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "AbstractCancerTestSuite.hpp"


/**
 * This class contains tests for methods on classes
 * inheriting from AbstractOdeBasedCellCycleModel.
 */
class TestOdeBasedCellCycleModels : public AbstractCancerTestSuite
{
public:

    void TestTysonNovakCellCycleModel() throw(Exception)
    {
        // Set up
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_timesteps = 50;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, num_timesteps);

        double standard_divide_time = 75.19/60.0;

        // Test TysonNovakCellCycleModel methods for a healthy cell
        TysonNovakCellCycleModel* p_cell_model = new TysonNovakCellCycleModel;
        p_cell_model->SetBirthTime(p_simulation_time->GetTime());
        TissueCell cell(STEM, HEALTHY, p_cell_model);

        // Test the cell is ready to divide at the right time
        for (unsigned i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            bool result = p_cell_model->ReadyToDivide();

            if (time > standard_divide_time)
            {
                TS_ASSERT_EQUALS(result, true);
            }
            else
            {
                TS_ASSERT_EQUALS(result, false);
            }
        }

        // Test ODE solution
        std::vector<double> proteins = p_cell_model->GetProteinConcentrations();
        TS_ASSERT(proteins.size()==6);
        TS_ASSERT_DELTA(proteins[0], 0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(proteins[1], 0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(proteins[2], 1.54216806705641, 1e-1);
        TS_ASSERT_DELTA(proteins[3], 1.40562614481544, 2e-2);
        TS_ASSERT_DELTA(proteins[4], 0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(proteins[5], 0.95328206604519, 2e-2);

        // For coverage, we also test TysonNovakCellCycleModel methods for a mutant cell
        p_cell_model->ResetForDivision();
        TysonNovakCellCycleModel* p_cell_model2 = static_cast<TysonNovakCellCycleModel*> (p_cell_model->CreateCellCycleModel());
        TissueCell stem_cell_2(STEM, APC_ONE_HIT, p_cell_model2);

        // Test the cell is ready to divide at the right time
        for (unsigned i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            bool result = p_cell_model->ReadyToDivide();
            bool result2 = p_cell_model2->ReadyToDivide();

            if (time> 2.0* standard_divide_time)
            {
                TS_ASSERT_EQUALS(result, true);
                TS_ASSERT_EQUALS(result2, true);
            }
            else
            {
                TS_ASSERT_EQUALS(result, false);
                TS_ASSERT_EQUALS(result2, false);
            }
        }

        // Test ODE solution
        proteins = p_cell_model->GetProteinConcentrations();
        TS_ASSERT_EQUALS(proteins.size(), 6u);
        TS_ASSERT_DELTA(proteins[0],0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(proteins[1],0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(proteins[2],1.54216806705641, 1e-1);
        TS_ASSERT_DELTA(proteins[3],1.40562614481544, 1e-1);
        TS_ASSERT_DELTA(proteins[4],0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(proteins[5],0.9662, 1e-2);
    }

    /**
     * Test for Tyson & Novak self-cycling cells without having their
     * initial conditions reset. When using CVODE, the cell cycle model
     * resets itself by halving the mass of the cell.
     * When not using CVODE, the cell cycle model resets its initial
     * conditions, since the oscillatory solution computed using the Chaste
     * ODE solver is not stable.
     */
    void TestTysonNovakCellCycleModelSolver() throw(Exception)
    {
        // Set up simulation time
        unsigned num_timesteps = 100000;
        double standard_divide_time = 75.19/60.0;
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(100.1*standard_divide_time, num_timesteps);

        // Create cell cycle model and associated cell
        TysonNovakCellCycleModel* p_repeating_cell_model = new TysonNovakCellCycleModel;
        TissueCell tyson_novak_cell(STEM, APC_ONE_HIT, p_repeating_cell_model);

        // Run through the cell cycle model for a certain duration
        // and test how many times it has stopped for division
        unsigned num_divisions = 0u;
        for (unsigned i=0; i<num_timesteps; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            bool result = p_repeating_cell_model->ReadyToDivide();
//                std::vector<double> proteins = p_repeating_cell_model->GetProteinConcentrations();
//                out << SimulationTime::Instance()->GetTime() << "\t";
//                for (unsigned j=0 ; j<proteins.size(); j++)
//                {
//                    out << proteins[j] << "\t";
//                }
//                out << "\n" << std::flush;

            if (result)
            {
                p_repeating_cell_model->ResetForDivision();
                p_repeating_cell_model->SetBirthTime(SimulationTime::Instance()->GetTime());
                num_divisions++;
            }
        }
        std::cout << num_divisions << "\n" << std::flush;
        TS_ASSERT_LESS_THAN(num_divisions,102u);
        TS_ASSERT_LESS_THAN(99u, num_divisions);
//            out.close();
        /*
         * Matlab code for plotting the output commented above:
         * cdchaste
         * data = load('TN_output.txt');
         * figure
         * for i=1:6
         *   subplot(3,2,i)
         *   plot(data(:,1),data(:,1+i))
         * end
         */
    }

    /**
     * In this test we use a WntCellCycleModel and begin with a steady-state
     * Wnt concentration of 1.0. Under such circumstances, the cell cycle model
     * would normally go into S phase at time t=5.971. Instead, we reduce the
     * Wnt concentration linearly to zero over the time interval 1<t<2, and the
     * cell doesn't divide.
     */
    void TestWntCellCycleModelForVaryingWntStimulus() throw(Exception)
    {
        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double end_time = 10.0 + CancerParameters::Instance()->GetMDuration(); // hours
        int num_timesteps = 1000*(int)end_time;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps); // 15.971 hours to go into S phase

        // Set up Wnt concentration
        double wnt_level = 1.0;
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Create cell cycle model and associated cell
        WntCellCycleModel* p_cell_model = new WntCellCycleModel();
        TissueCell stem_cell(STEM, HEALTHY, p_cell_model);
        stem_cell.InitialiseCellCycleModel();

        // When using a WntCellCycleModel, there is no such thing as
        // a 'stem cell'. Cell type is changed to transit or
        // differentiated, depending on the Wnt concentration, when
        // InitialiseCellCycleModel() is called.
        TS_ASSERT_EQUALS(stem_cell.GetCellType(), TRANSIT);

        // Progress through the cell cycle
        for (int i=0; i<num_timesteps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();
            bool result = p_cell_model->ReadyToDivide();

            // The Wnt concentration reduces from 1 to 0 over
            // the interval 1<t<2 (at beginning of G1 phase)
            if (time <= 2.0)
            {
                wnt_level = 2.0-time;
            }
            else
            {
                wnt_level = 0.0;
            }
            WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

            // Test that the cell cycle model does not stop for division
            TS_ASSERT_EQUALS(result, false);
        }

        // Test ODE solution
        double tol = 1e-5;

        std::vector<double> test_results = p_cell_model->GetProteinConcentrations();
#ifdef CHASTE_CVODE
        TS_ASSERT_DELTA(test_results[0], 0.7329922345, tol);
#else
        TS_ASSERT_DELTA(test_results[0], 7.330036281693106e-01, tol);
#endif //CHASTE_CVODE
        TS_ASSERT_DELTA(test_results[1], 1.715690244022676e-01, tol);
        TS_ASSERT_DELTA(test_results[2], 6.127460817296076e-02, tol);
        TS_ASSERT_DELTA(test_results[3], 1.549402358669023e-07, tol);
        TS_ASSERT_DELTA(test_results[4], 4.579067802591843e-08, tol);
        TS_ASSERT_DELTA(test_results[5], 9.999999999999998e-01, tol);
        TS_ASSERT_DELTA(test_results[6], 0.5*7.415537855270896e-03, tol);
        TS_ASSERT_DELTA(test_results[7], 0.5*7.415537855270896e-03, tol);
        TS_ASSERT_DELTA(test_results[8], 0.0, tol);

        // Test that, since the cell now experiences a low Wnt concentration,
        // it has indeed changed cell type to differentiated
        TS_ASSERT_EQUALS(stem_cell.GetCellType(), DIFFERENTIATED);

        double diff = 1.0;
        test_results[6] = test_results[6] + diff;

        p_cell_model->SetProteinConcentrationsForTestsOnly(1.0, test_results);

        test_results = p_cell_model->GetProteinConcentrations();

        TS_ASSERT_DELTA(test_results[6], diff + 0.5*7.415537855270896e-03, 1e-5);
        TS_ASSERT_DELTA(test_results[5], 0.9999999999999998, 1e-5);

        // Tidy up
        WntConcentration::Destroy();
    }


    void TestIngeWntSwatCellCycleModel() throw(Exception)
    {
        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double end_time = 30; // hours
        int num_timesteps = 100*(int)end_time;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps); // 15.971 hours to go into S phase

        // Set up Wnt concentration
        double wnt_level = 1.0;
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Cover exception - when constructing an instance of IngeWntSwatCellCycleModel,
        // we must pass in an hypothesis number (1 or 2)
        TS_ASSERT_THROWS_ANYTHING(IngeWntSwatCellCycleModel model(0));

        // Create cell cycle model and associated cell
        IngeWntSwatCellCycleModel* p_cell_model = new IngeWntSwatCellCycleModel(1);

        // Test that member variables are set correctly
        TS_ASSERT_EQUALS(p_cell_model->GetHypothesis(), 1u);

        TissueCell stem_cell(STEM, HEALTHY, p_cell_model);

        // Coverage of cell cycle model copying without an ODE system set up
        TissueCell stem_cell2 = stem_cell;
        TS_ASSERT_EQUALS(stem_cell2.GetMutationState(), HEALTHY);

        stem_cell.InitialiseCellCycleModel();

        // When using a WntCellCycleModel, there is no such thing as
        // a 'stem cell'. Cell type is changed to transit or
        // differentiated, depending on the Wnt concentration, when
        // InitialiseCellCycleModel() is called.
        TS_ASSERT_EQUALS(stem_cell.GetCellType(), TRANSIT);

        WntConcentration::Instance()->SetConstantWntValueForTesting(1.0);

        double tol = 1e-4;
#ifdef CHASTE_CVODE
        const double expected_g1_duration = 6.18252;
        tol  = 1e-5;
#else
        const double expected_g1_duration = 6.1959;
#endif //CHASTE_CVODE

        // Progress through the cell cycle under a constant Wnt concentration
        for (int i=0; i<21*num_timesteps/30.0; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // Call ReadyToDivide on the cell, then test the results
            // of calling ReadyToDivide on the model and test (in
            // CheckReadyToDivideAndPhaseIsUpdated).
            stem_cell.ReadyToDivide();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model, expected_g1_duration);
        }

        // Test that the cell cycle model stopped for division correctly
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 21.0, 1e-4);
        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), true);

        std::vector<double> test_results = p_cell_model->GetProteinConcentrations();

        // Test ODE solution (the correct values were found to an accuracy
        // of around 1e-6 using the Matlab solver ode15s)
        TS_ASSERT_DELTA(test_results[0] , 2.93699961539512e-01 , 2*10*tol);
        TS_ASSERT_DELTA(test_results[1] , 1.000000000000000, 2*100*tol);
        TS_ASSERT_DELTA(test_results[2] , 2.40050625298734 , 2*10*tol);
        TS_ASSERT_DELTA(test_results[3] , 1.39281551739568 , 2*10*tol);
        TS_ASSERT_DELTA(test_results[4] , 1.35836451056026e-01 , 2*10*tol);
        TS_ASSERT_DELTA(test_results[5] , 1.428571428571429e-01, tol);
        TS_ASSERT_DELTA(test_results[6] , 2.857142857142857e-02, tol);
        TS_ASSERT_DELTA(test_results[7] , 2.120643654085205e-01, tol);
        TS_ASSERT_DELTA(test_results[8] , 1.439678172957377e+01, 10*tol);
        TS_ASSERT_DELTA(test_results[9] ,                     0, tol);
        TS_ASSERT_DELTA(test_results[10],                     0, tol);
        TS_ASSERT_DELTA(test_results[11],                     0, tol);
        TS_ASSERT_DELTA(test_results[12], 1.000000000000002e+01, tol);
        TS_ASSERT_DELTA(test_results[13], 1.028341552112414e+02, 100*tol);
        TS_ASSERT_DELTA(test_results[14],                     0, tol);
        TS_ASSERT_DELTA(test_results[15], 2.499999999999999e+01, tol);
        TS_ASSERT_DELTA(test_results[16], 1.439678172957377e+01, 10*tol);
        TS_ASSERT_DELTA(test_results[17],                     0, tol);
        TS_ASSERT_DELTA(test_results[18],                     0, tol);
        TS_ASSERT_DELTA(test_results[19],                     0, tol);
        TS_ASSERT_DELTA(test_results[20], 2.235636835087684, tol);
        TS_ASSERT_DELTA(test_results[21], 1.000000000000000, tol);

        // The cell cycle model acts as if it was divided at time = 16.1877. This
        // is fine as the cell cycle model dictates the division time, not when
        // the cell is actually divided.
        TissueCell daughter_cell = stem_cell.Divide();
        AbstractCellCycleModel* p_cell_model2 = daughter_cell.GetCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model->GetCurrentCellCyclePhase(), M_PHASE);
        TS_ASSERT_EQUALS(p_cell_model2->GetCurrentCellCyclePhase(), M_PHASE);

        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_cell_model2->ReadyToDivide(), false);

        // Check that the first five protein levels have been reset and the rest are the same as before.
        tol = 1e-4;
        test_results = p_cell_model->GetProteinConcentrations();
        TS_ASSERT_DELTA(test_results[0] , 2.631865125420296e-01, 10*tol);
        TS_ASSERT_DELTA(test_results[1] , 2.678271949808561e-01, 10*tol);
        TS_ASSERT_DELTA(test_results[2] , 2.389956081120099, 10*tol);
        TS_ASSERT_DELTA(test_results[3] , 1.390258620103223, 10*tol);
        TS_ASSERT_DELTA(test_results[4] , 1.218603203963113e-01, 10*tol);
        TS_ASSERT_DELTA(test_results[5] , 1.428571428571429e-01, tol);
        TS_ASSERT_DELTA(test_results[6] , 2.857142857142857e-02, tol);
        TS_ASSERT_DELTA(test_results[7] , 2.120643654085205e-01, tol);
        TS_ASSERT_DELTA(test_results[8] , 1.439678172957377e+01, 10*tol);
        TS_ASSERT_DELTA(test_results[9] ,                     0, tol);
        TS_ASSERT_DELTA(test_results[10],                     0, tol);
        TS_ASSERT_DELTA(test_results[11],                     0, tol);
        TS_ASSERT_DELTA(test_results[12], 1.000000000000002e+01, tol);
        TS_ASSERT_DELTA(test_results[13], 1.028341552112414e+02, 100*tol);
        TS_ASSERT_DELTA(test_results[14],                     0, tol);
        TS_ASSERT_DELTA(test_results[15], 2.499999999999999e+01, tol);
        TS_ASSERT_DELTA(test_results[16], 1.439678172957377e+01, 10*tol);
        TS_ASSERT_DELTA(test_results[17],                     0, tol);
        TS_ASSERT_DELTA(test_results[18],                     0, tol);
        TS_ASSERT_DELTA(test_results[19],                     0, tol);
        TS_ASSERT_DELTA(test_results[20], 2.235636835087684, tol);
        TS_ASSERT_DELTA(test_results[21], 1.000000000000000, tol);

        // Now progress through the cell cycle under a decreasing Wnt concentration
        for (int i=0; i<9*num_timesteps/30.0; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();
            bool result = p_cell_model->ReadyToDivide();
            bool result2 = p_cell_model2->ReadyToDivide();

            if (time <= 22.0)
            {
                wnt_level = 22.0 - time;
            }
            else
            {
                wnt_level = 0.0;
            }
            WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

            TS_ASSERT_EQUALS(result, false);
            TS_ASSERT_EQUALS(result2, false);
        }
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 30.0, 1e-4);

        test_results = p_cell_model->GetProteinConcentrations();

        // Test ODE solution (the tolerances for some values are deliberately
        // loose, so that we don't need different answers for CVODE and RK4)
        TS_ASSERT_DELTA(test_results[0] , 0.3648, 2e-3);
        TS_ASSERT_DELTA(test_results[1] , 1.000, 1e-2);
        TS_ASSERT_DELTA(test_results[2] , 1.4955, 1e-2);
        TS_ASSERT_DELTA(test_results[3] , 0.8125, 2e-2);
        TS_ASSERT_DELTA(test_results[4] , 0.0996, 1e-2);
        TS_ASSERT_DELTA(test_results[5] , 0.6666, 1e-4);
        TS_ASSERT_DELTA(test_results[6] , 0.0666, 1e-4);
        TS_ASSERT_DELTA(test_results[7] , 0.7311, 2e-3);
        TS_ASSERT_DELTA(test_results[8] , 6.0299, 4e-2);
        TS_ASSERT_DELTA(test_results[9] , 0, 1e-4);
        TS_ASSERT_DELTA(test_results[10], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[11], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[12], 17.5407, 1e-2);
        TS_ASSERT_DELTA(test_results[13], 75.5926, .5);
        TS_ASSERT_DELTA(test_results[14], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[15], 29.5666, 1e-2);
        TS_ASSERT_DELTA(test_results[16], 7.1333, 1e-1);
        TS_ASSERT_DELTA(test_results[17], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[18], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[19], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[20], 1.5991, 2e-2);
        TS_ASSERT_DELTA(test_results[21], 0.0000, 1e-4);

        // Test that, since the cell now experiences a low Wnt concentration,
        // it has indeed changed cell type to differentiated
        TS_ASSERT_EQUALS(stem_cell.GetCellType(), DIFFERENTIATED);

        // Test beta catenin levels

        // membrane_beta_cat = Ca + Ma
        double membrane_beta_cat = test_results[13]+test_results[14];

        // cytoplasmic_beta_cat = Cu + Co + Cc + Mo + Mc
        double cytoplasm_beta_cat = test_results[7] + test_results[8]
                          + test_results[9] + test_results[10]+test_results[11];

        // nuclear_beta_cat = Cot + Cct + Mot + Mct
        double nuclear_beta_cat = test_results[16] + test_results[17]
                                    + test_results[18] + test_results[19];

        TS_ASSERT_DELTA(p_cell_model->GetMembraneBoundBetaCateninLevel(), membrane_beta_cat, 1e-4);
        TS_ASSERT_DELTA(p_cell_model->GetCytoplasmicBetaCateninLevel(), cytoplasm_beta_cat, 1e-4);
        TS_ASSERT_DELTA(p_cell_model->GetNuclearBetaCateninLevel(), nuclear_beta_cat, 1e-4);

        // Tidy up
        WntConcentration::Destroy();
    }


    void TestWntCellCycleModelForAPCSingleHit() throw(Exception)
    {
        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        int num_timesteps = 500;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps); // 15.971 hours to go into S phase

        // Set up Wnt concentration
        double wnt_level = 1.0;
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Create cell cycle model and associated cell
        WntCellCycleModel* p_cell_model = new WntCellCycleModel();
        TissueCell stem_cell(STEM, HEALTHY, p_cell_model);
        stem_cell.InitialiseCellCycleModel();

        double SG2M_duration = CancerParameters::Instance()->GetSG2MDuration();
        TS_ASSERT_THROWS_NOTHING(WntCellCycleModel cell_model_3());

        // Create another cell cycle model and associated cell
        WntCellCycleModel* p_cell_model_1 = new WntCellCycleModel();
        TissueCell stem_cell_1(STEM, APC_ONE_HIT, p_cell_model_1);
        stem_cell_1.InitialiseCellCycleModel();

        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 4.804 hrs and then finish dividing
        // 10 hours later at 14.804 hours.
#ifdef CHASTE_CVODE
        double expected_g1_duration = 4.79772;
#else
        double expected_g1_duration = 4.8084;
#endif //CHASTE_CVODE
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_1, expected_g1_duration);
        }

        p_cell_model_1->ResetForDivision();
#ifdef CHASTE_CVODE
        expected_g1_duration = 4.80678;
#else
        expected_g1_duration = 4.8084;
#endif //CHASTE_CVODE

        TS_ASSERT_DELTA(SG2M_duration, 10.0, 1e-5);
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_1, expected_g1_duration);
        }

        // Tidy up
        WntConcentration::Destroy();
    }


    void TestWntCellCycleModelForBetaCatSingleHit() throw(Exception)
    {
        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        int num_timesteps = 500;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps); // 15.971 hours to go into S phase

        // Set up Wnt concentration
        double wnt_level = 0.0;
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Create cell cycle model and associated cell
        WntCellCycleModel* p_cell_model = new WntCellCycleModel();
        TissueCell stem_cell(STEM, BETA_CATENIN_ONE_HIT, p_cell_model);
        stem_cell.InitialiseCellCycleModel();

        TS_ASSERT_THROWS_NOTHING(WntCellCycleModel cell_model_3());

        // Create another cell cycle model and associated cell
        WntCellCycleModel* p_cell_model_1 = new WntCellCycleModel();
        TissueCell stem_cell_1(STEM, BETA_CATENIN_ONE_HIT, p_cell_model_1);
        stem_cell_1.InitialiseCellCycleModel();

        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 7.82 hrs and then finish dividing
        // 10 hours later at 17.82 hours.
#ifdef CHASTE_CVODE
        double expected_g1_duration = 7.81718;
#else
        double expected_g1_duration = 7.8342;
#endif //CHASTE_CVODE

        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_THROWS_ANYTHING(p_cell_model_1->UpdateCellType());
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_1, expected_g1_duration);
        }

        p_cell_model_1->ResetForDivision();
#ifdef CHASTE_CVODE
        expected_g1_duration = 7.81873;
#else
        expected_g1_duration = 7.8342;
#endif //CHASTE_CVODE

        // Test progress through the cell cycle
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_1, expected_g1_duration);
        }

        // Tidy up
        WntConcentration::Destroy();
    }


    void TestWntCellCycleModelForAPCDoubleHit() throw(Exception)
    {
        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        int num_timesteps = 500;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps); // 15.971 hours to go into S phase

        // Set up Wnt concentration
        double wnt_level = 0.738; // the Wnt concentrationshouldn't matter for a cell with APC double hit
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Create cell cycle model and associated cell
        WntCellCycleModel* p_cell_model_1 = new WntCellCycleModel();
        TissueCell stem_cell_1(STEM, APC_TWO_HIT, p_cell_model_1);
        stem_cell_1.InitialiseCellCycleModel();

        // Create another cell cycle model and associated cell
        WntCellCycleModel* p_cell_model_2 = new WntCellCycleModel();
        TissueCell stem_cell_2(STEM, APC_TWO_HIT, p_cell_model_2);
        stem_cell_2.InitialiseCellCycleModel();

        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 3.943 hrs and then finish dividing
        // 10 hours later at 13.9435 hours.
#ifdef CHASTE_CVODE
        double expected_g1_duration = 3.93959;
#else
        double expected_g1_duration = 3.9455;
#endif //CHASTE_CVODE

        // Test progress through the cell cycle
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_2, expected_g1_duration);
        }

        p_cell_model_2->ResetForDivision();

#ifdef CHASTE_CVODE
        expected_g1_duration = 3.94529;
#else
        expected_g1_duration = 3.9455;
#endif //CHASTE_CVODE

        // Test progress through the cell cycle
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_2, expected_g1_duration);
        }

        // Tidy up
        WntConcentration::Destroy();
    }


    void TestWntCellCycleModelForConstantWntStimulusHealthyCell() throw(Exception)
    {
        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        int num_timesteps = 500;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase

        // Set up Wnt concentration
        double wnt_level = 1.0;
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Create cell cycle model and associated cell
        WntCellCycleModel* p_cell_model_1 = new WntCellCycleModel();
        TissueCell stem_cell_1(STEM, HEALTHY, p_cell_model_1);
        stem_cell_1.InitialiseCellCycleModel();

        // Create another cell cycle model and associated cell
        WntCellCycleModel* p_cell_model_2 = new WntCellCycleModel();
        TissueCell stem_cell_2(STEM, HEALTHY, p_cell_model_2);
        stem_cell_2.InitialiseCellCycleModel();

        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 5.971 hrs and then finish dividing
        // 10 hours later at 15.971 hours.
#ifdef CHASTE_CVODE
        double expected_g1_duration = 5.96441;
#else
        double expected_g1_duration = 5.9782;
#endif //CHASTE_CVODE

        // Test progress through the cell cycle
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_2, expected_g1_duration);
        }

        p_cell_model_2->ResetForDivision();
#ifdef CHASTE_CVODE
        expected_g1_duration = 5.96016;
#else
        expected_g1_duration = 5.9782;
#endif //CHASTE_CVODE

        // Test progress through the cell cycle
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_2, expected_g1_duration);
        }

        // Tidy up
        WntConcentration::Destroy();
    }


    void TestStochasticWntCellCycleModel() throw (Exception)
    {
        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        int num_timesteps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(20, num_timesteps);// 15.971 hours to go into S phase

        // Set up Wnt concentration
        double wnt_level = 1.0;
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Create cell cycle model and associated cell
        StochasticWntCellCycleModel* p_cell_model = new StochasticWntCellCycleModel();
        TissueCell stem_cell(STEM, HEALTHY, p_cell_model);
        stem_cell.InitialiseCellCycleModel();

        // A WntCellCycleModel does this:
        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 5.971 hrs and then finish dividing
        // 10 hours later at 15.971 hours.
        //
        // A StochasticWntCellCycleModel does this:
        // divides at the same time with a random normal distribution
        // for the SG2M time (default 10) in this case 9.0676

        // Test progress through the cell cycle
        for (int i=0; i<num_timesteps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();
            bool result = p_cell_model->ReadyToDivide();

            if (time < 5.971 + 9.0676)
            {
                TS_ASSERT_EQUALS(result, false);
            }
            else
            {
                TS_ASSERT_EQUALS(result, true);
            }
        }

        // Tidy up
        WntConcentration::Destroy();
    }


    void TestArchiveTysonNovakCellCycleModels()
    {
        // Set up
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "tyson_novak_cell_cycle.arch";

        {
            // Set up simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(100.0, 1);

            // Create cell cycle model and associated cell
            TysonNovakCellCycleModel* p_model = new TysonNovakCellCycleModel;

            p_simulation_time->IncrementTimeOneStep();

            TissueCell cell(TRANSIT, HEALTHY, p_model);
            cell.InitialiseCellCycleModel();

            TS_ASSERT_EQUALS(p_model->ReadyToDivide(),true);

            p_model->SetBirthTime(-1.0);

            // Create an ouput archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive cell
            TissueCell* const p_cell = &cell;
            output_arch << p_cell;

            SimulationTime::Destroy();
        }

        {
            // Set up simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            TissueCell* p_cell;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_cell;

            AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();

            // Test archiving
            TS_ASSERT_EQUALS(p_model->ReadyToDivide(),true);
            TS_ASSERT_DELTA(p_model->GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(),101.0,1e-12);
            delete p_cell;
        }
    }


    void TestArchiveWntCellCycleModel()
    {
        // Set up
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "wnt_cell_cycle.arch";
        WntConcentration::Instance()->SetConstantWntValueForTesting(1.0);

        {
            // Set up simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(16, 2);

            // Create cell cycle model and associated cell
            WntCellCycleModel* p_cell_model = new WntCellCycleModel();
            TissueCell stem_cell(STEM, HEALTHY, p_cell_model);
            stem_cell.InitialiseCellCycleModel();

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(),false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(),true);

            // Should be in G2 after a couple of timesteps
            TS_ASSERT_EQUALS(p_cell_model->GetCurrentCellCyclePhase(), G_TWO_PHASE);

            stem_cell.GetCellCycleModel()->SetBirthTime(-1.0);

            // Create an ouput archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive cell
            TissueCell* const p_cell = &stem_cell;
            output_arch << p_cell;

            SimulationTime::Destroy();
        }

        {
            // Set up simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            CancerParameters* p_inst1 = CancerParameters::Instance();

            p_inst1->SetSDuration(101.0);

            TissueCell* p_cell;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_cell;

            // Test archiving
            AbstractCellCycleModel* p_cell_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());

            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), true);
            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_cell_model->GetAge(), 17.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetSG2MDuration(), 10.0, 1e-12);
            TS_ASSERT_EQUALS(p_cell_model->GetCurrentCellCyclePhase(), G_TWO_PHASE);

            delete p_cell;
        }

        // Tidy up
        WntConcentration::Destroy();
    }


    void TestArchiveIngeWntSwatCellCycleModel()
    {
        // Set up
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "inge_wnt_swat_cell_cycle.arch";
        WntConcentration::Instance()->SetConstantWntValueForTesting(1.0);

        {
            // Set up simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(17, 2);

            // Create cell cycle model and associated cell
            IngeWntSwatCellCycleModel* p_cell_model = new IngeWntSwatCellCycleModel(1);
            TissueCell stem_cell(STEM, HEALTHY, p_cell_model);
            stem_cell.InitialiseCellCycleModel();

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(),false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(),true);

            stem_cell.GetCellCycleModel()->SetBirthTime(-1.0);

            // Create an ouput archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive cell
            TissueCell* const p_cell = &stem_cell;
            output_arch << p_cell;

            SimulationTime::Destroy();
        }

        {
            // Set up simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            CancerParameters* p_inst1 = CancerParameters::Instance();

            p_inst1->SetSDuration(101.0);

            TissueCell* p_cell;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_cell;

            // Test archiving
            AbstractCellCycleModel* p_cell_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());

            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), true);
            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_cell_model->GetAge(), 18.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetSG2MDuration(), 10.0, 1e-12);
            TS_ASSERT_EQUALS(p_cell_model->GetCurrentCellCyclePhase(), G_TWO_PHASE);

            delete p_cell;
        }

        // Tidy up
        WntConcentration::Destroy();
    }


    void TestArchiveStochasticWntCellCycleModels()
    {
        // In this case the first cycle time will be 5.971+9.0676 = 15.0386
        // note that the S-G2-M time is assigned when the cell finishes G1
        //(i.e. at time 5.971 here so the model has to be archived BEFORE that.

        // Set up
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "stochastic_wnt_cell_cycle.arch";
        WntConcentration::Instance()->SetConstantWntValueForTesting(1.0);

        {
            // In this test the RandomNumberGenerator in existence

            // Set up simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(16.0, 1000);

            // Create cell cycle model and associated cell
            StochasticWntCellCycleModel* p_stoc_model = new StochasticWntCellCycleModel();
            TissueCell stoc_cell(STEM, HEALTHY, p_stoc_model);
            stoc_cell.InitialiseCellCycleModel();

            // Create another cell cycle model and associated cell
            WntCellCycleModel* p_wnt_model = new WntCellCycleModel();
            TissueCell wnt_cell(STEM, HEALTHY, p_wnt_model);
            wnt_cell.InitialiseCellCycleModel();

            p_simulation_time->IncrementTimeOneStep(); // 5.5

            while (p_simulation_time->GetTime() < 4.0)
            {
                p_simulation_time->IncrementTimeOneStep();
            }

            TS_ASSERT_EQUALS(stoc_cell.GetCellCycleModel()->ReadyToDivide(),false);
            TS_ASSERT_EQUALS(wnt_cell.GetCellCycleModel()->ReadyToDivide(),false);
            TS_ASSERT_EQUALS(stoc_cell.GetCellCycleModel()->GetCurrentCellCyclePhase(), G_ONE_PHASE);

            // When the above tests are included here they pass, so we
            // also put them after the load to see if they still pass.

            // Create an ouput archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive cells
            TissueCell* const p_wnt_cell = &wnt_cell;
            TissueCell* const p_stoc_cell = &stoc_cell;

            output_arch << p_stoc_cell;
            output_arch << p_wnt_cell;

            SimulationTime::Destroy();
        }

        {
            // Set up simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(16.0, 2);

            CancerParameters* p_inst1 = CancerParameters::Instance();

            p_inst1->SetSDuration(101.0);

            TissueCell* p_stoc_cell;
            TissueCell* p_wnt_cell;

            std::vector<double> cell_cycle_influence1;
            cell_cycle_influence1.push_back(1.0);

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_stoc_cell;
            input_arch >> p_wnt_cell;

            // Test archiving
            TS_ASSERT_EQUALS(p_stoc_cell->GetCellCycleModel()->GetCurrentCellCyclePhase(), G_ONE_PHASE);

            // Check - stochastic should divide at 15.03
            // Wnt should divide at 15.971
            while (p_simulation_time->GetTime() < 15.0)
            {
                p_simulation_time->IncrementTimeOneStep();
            }
            TS_ASSERT_EQUALS(p_stoc_cell->GetCellCycleModel()->ReadyToDivide(), false);
            TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->ReadyToDivide(), false);

            while (p_simulation_time->GetTime() < 15.5)
            {
                p_simulation_time->IncrementTimeOneStep();
            }
            TS_ASSERT_EQUALS(p_stoc_cell->GetCellCycleModel()->ReadyToDivide(), true); // only for stochastic
            TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->ReadyToDivide(), false);

            while (p_simulation_time->GetTime() < 16.0)
            {
                p_simulation_time->IncrementTimeOneStep();
            }
            TS_ASSERT_EQUALS(p_stoc_cell->GetCellCycleModel()->ReadyToDivide(), true);
            TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->ReadyToDivide(), true);

            TS_ASSERT_DELTA(p_stoc_cell->GetCellCycleModel()->GetBirthTime(), 0.0, 1e-12);
            TS_ASSERT_DELTA(p_stoc_cell->GetCellCycleModel()->GetAge(), 16.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetSG2MDuration(), 10.0, 1e-12);

            delete p_stoc_cell;
            delete p_wnt_cell;
        }

        // Tidy up
        WntConcentration::Destroy();
    }
    
    void TestCopyingCells() throw(Exception)
    {
        TysonNovakCellCycleModel *p_original_cell_cycle= new TysonNovakCellCycleModel;
        TissueCell cell(STEM, HEALTHY, p_original_cell_cycle);
        
        // These changed from default
        p_original_cell_cycle->mLastTime = 42.01;
        p_original_cell_cycle->mDivideTime = -42.01;
        p_original_cell_cycle->mFinishedRunningOdes = true;
        p_original_cell_cycle->mG2PhaseStartTime = 124;
        p_original_cell_cycle->mBirthTime = 765;
        p_original_cell_cycle->mCurrentCellCyclePhase = S_PHASE;
        p_original_cell_cycle->mG1Duration = 987;
        p_original_cell_cycle->mReadyToDivide = true;
        
        TissueCell cell2 = cell;
        TysonNovakCellCycleModel *p_new_cell_cycle = static_cast<TysonNovakCellCycleModel* > (cell2.GetCellCycleModel());
        TS_ASSERT_EQUALS(&cell2, p_new_cell_cycle->GetCell() );
        TS_ASSERT_DIFFERS(&cell, &cell2);
        TS_ASSERT_DIFFERS(p_original_cell_cycle, p_new_cell_cycle);
        TS_ASSERT_EQUALS(p_original_cell_cycle->mLastTime, p_new_cell_cycle->mLastTime);
        TS_ASSERT_EQUALS(p_original_cell_cycle->mDivideTime, p_new_cell_cycle->mDivideTime);
        TS_ASSERT_EQUALS(p_original_cell_cycle->mFinishedRunningOdes, p_new_cell_cycle->mFinishedRunningOdes);
        TS_ASSERT_EQUALS(p_original_cell_cycle->mG2PhaseStartTime, p_new_cell_cycle->mG2PhaseStartTime);
        TS_ASSERT_EQUALS(p_original_cell_cycle->mBirthTime, p_new_cell_cycle->mBirthTime);
        TS_ASSERT_EQUALS(p_original_cell_cycle->mCurrentCellCyclePhase, p_new_cell_cycle->mCurrentCellCyclePhase);
        TS_ASSERT_EQUALS(p_original_cell_cycle->mG1Duration, p_new_cell_cycle->mG1Duration);
        TS_ASSERT_EQUALS(p_original_cell_cycle->mReadyToDivide, p_new_cell_cycle->mReadyToDivide);
        
    }
    
};


#endif /*TESTODEBASEDCELLCYCLEMODELS_HPP_*/
