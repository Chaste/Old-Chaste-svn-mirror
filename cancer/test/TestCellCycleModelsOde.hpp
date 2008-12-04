/*

Copyright (C) University of Oxford, 2008

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
#ifndef TESTCELLCYCLEMODELSODE_HPP_
#define TESTCELLCYCLEMODELSODE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "Alarcon2004OxygenBasedCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "IngeWntSwatCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestOdeCellCycleModels : public AbstractCancerTestSuite
{
public:

    void TestTysonNovakCellCycleModel(void) throw(Exception)
    {
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        int num_timesteps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, num_timesteps);

        double standard_divide_time = 75.19/60.0;

        TysonNovakCellCycleModel* p_cell_model = new TysonNovakCellCycleModel;
        // Coverage
        p_cell_model->SetBirthTime(p_simulation_time->GetDimensionalisedTime());
        TissueCell cell(STEM, HEALTHY, p_cell_model);

        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();

            bool result = p_cell_model->ReadyToDivide();

            if (time>standard_divide_time)
            {
                TS_ASSERT(result==true);
            }
            else
            {
                TS_ASSERT(result==false);
            }
        }

        std::vector<double> proteins = p_cell_model->GetProteinConcentrations();

        TS_ASSERT(proteins.size()==6);

        TS_ASSERT_DELTA(proteins[0],0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(proteins[1],0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(proteins[2],1.54216806705641, 1e-2);
        TS_ASSERT_DELTA(proteins[3],1.40562614481544, 1e-2);
        TS_ASSERT_DELTA(proteins[4],0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(proteins[5],0.95328206604519, 1e-2);

        p_cell_model->ResetForDivision();
        TysonNovakCellCycleModel *p_cell_model2 = static_cast <TysonNovakCellCycleModel*> (p_cell_model->CreateCellCycleModel());

        TissueCell stem_cell_2(STEM, APC_ONE_HIT, p_cell_model2);

        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();

            bool result = p_cell_model->ReadyToDivide();
            bool result2 = p_cell_model2->ReadyToDivide();

            if (time> 2.0* standard_divide_time)
            {
                TS_ASSERT_EQUALS(result,true);
                TS_ASSERT_EQUALS(result2,true);
            }
            else
            {
                TS_ASSERT_EQUALS(result,false);
                TS_ASSERT_EQUALS(result2,false);
            }
        }

        proteins = p_cell_model->GetProteinConcentrations();

        TS_ASSERT_EQUALS(proteins.size(),6u);

        TS_ASSERT_DELTA(proteins[0],0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(proteins[1],0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(proteins[2],1.54216806705641, 1e-2);
        TS_ASSERT_DELTA(proteins[3],1.40562614481544, 1e-1);
        TS_ASSERT_DELTA(proteins[4],0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(proteins[5],0.95328206604519, 1e-2);
    }


    void TestWntCellCycleModelForVaryingWntStimulus() throw(Exception)
    {
        // Here we have a system at rest at Wnt = 1.0 - it would normally go into S phase at 5.971.
        // Instead we reduce Wnt linearly over 1<t<2 to zero and the cell doesn't divide.
        SimulationTime *p_simulation_time = SimulationTime::Instance();

        double end_time = 10.0 + CancerParameters::Instance()->GetMDuration(); //hours
        int num_timesteps = 1000*(int)end_time;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);// 15.971 hours to go into S phase

        double wnt_level = 1.0;
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

        WntCellCycleModel* p_cell_model = new WntCellCycleModel();

        TissueCell stem_cell(STEM, HEALTHY, p_cell_model);

        stem_cell.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(stem_cell.GetCellType(),TRANSIT);

        for (int i=0; i<num_timesteps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            bool result = p_cell_model->ReadyToDivide();

            // Reduces from 1 to 0 over the interval 1<t<2
            // (at beginning of G1 phase)
            if (time <= 2.0)
            {
                wnt_level = 2.0-time;
            }
            else
            {
                wnt_level = 0.0;
            }
            WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

            TS_ASSERT(result==false);
        }

        std::vector<double> test_results = p_cell_model->GetProteinConcentrations();
        TS_ASSERT_DELTA(test_results[0] , 7.330036281693106e-01 , 1e-5);
        TS_ASSERT_DELTA(test_results[1] , 1.715690244022676e-01 , 1e-5);
        TS_ASSERT_DELTA(test_results[2] , 6.127460817296076e-02 , 1e-5);
        TS_ASSERT_DELTA(test_results[3] , 1.549402358669023e-07 , 1e-5);
        TS_ASSERT_DELTA(test_results[4] , 4.579067802591843e-08 , 1e-5);
        TS_ASSERT_DELTA(test_results[5] , 9.999999999999998e-01 , 1e-5);
        TS_ASSERT_DELTA(test_results[6] , 0.5*7.415537855270896e-03 , 1e-5);
        TS_ASSERT_DELTA(test_results[7] , 0.5*7.415537855270896e-03 , 1e-5);
        TS_ASSERT_DELTA(test_results[8] , 0.0 , 1e-6);

        TS_ASSERT_EQUALS(stem_cell.GetCellType(), DIFFERENTIATED);

        double diff = 1.0;
        test_results[6] = test_results[6] + diff;

        p_cell_model->SetProteinConcentrationsForTestsOnly(1.0, test_results);

        test_results = p_cell_model->GetProteinConcentrations();

        TS_ASSERT_DELTA(test_results[6] , diff + 0.5*7.415537855270896e-03 , 1e-5);
        TS_ASSERT_DELTA(test_results[5] , 9.999999999999998e-01 , 1e-5);

        WntConcentration::Destroy();
    }


    void TestIngeWntSwatCellCycleModel() throw(Exception)
    {
        // Here we have a system at rest at Wnt = 1.0 - it would normally go into S phase at 5.971.
        // Instead we reduce Wnt linearly over 0<t<1 to zero and the cell doesn't divide.
        SimulationTime *p_simulation_time = SimulationTime::Instance();

        double end_time = 30; //hours
        int num_timesteps = 1000*(int)end_time;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);// 15.971 hours to go into S phase

        double wnt_level = 1.0;
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Cover exception
        TS_ASSERT_THROWS_ANYTHING(IngeWntSwatCellCycleModel model(0));

        IngeWntSwatCellCycleModel* p_cell_model = new IngeWntSwatCellCycleModel(1);
        TS_ASSERT(p_cell_model->UsesBetaCat());
        TS_ASSERT_EQUALS(p_cell_model->GetHypothesis(), 1u);

        TissueCell stem_cell(STEM, HEALTHY, p_cell_model);

        // Coverage of cell cycle model copying without an ODE system set up
        TissueCell stem_cell2 = stem_cell;
        TS_ASSERT_EQUALS(stem_cell2.GetMutationState(), HEALTHY);

        stem_cell.InitialiseCellCycleModel();

        // Check the Inge model has changed the cell type correctly.
        TS_ASSERT_EQUALS(stem_cell.GetCellType(),TRANSIT);
        WntConcentration::Instance()->SetConstantWntValueForTesting(1.0);
        for (int i=0; i<21*num_timesteps/30.0; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // Call ReadyToDivide on the cell, then test the results
            // of calling ReadyToDivide on the model and test (in
            // CheckReadyToDivideAndPhaseIsUpdated).
            stem_cell.ReadyToDivide();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model,6.1877);
        }

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetDimensionalisedTime(), 21.0, 1e-4);
        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), true);

        std::vector<double> test_results = p_cell_model->GetProteinConcentrations();

        TS_ASSERT_DELTA(test_results[0] , 2.937528298476724e-01, 1e-3);
        TS_ASSERT_DELTA(test_results[1] , 1.000000000000000e+00, 1e-2);
        TS_ASSERT_DELTA(test_results[2] , 2.400571020059542e+00, 1e-3);
        TS_ASSERT_DELTA(test_results[3] , 1.392891502782046e+00, 1e-3);
        TS_ASSERT_DELTA(test_results[4] , 1.358708742498684e-01, 1e-3);
        TS_ASSERT_DELTA(test_results[5] , 1.428571428571429e-01, 1e-4);
        TS_ASSERT_DELTA(test_results[6] , 2.857142857142857e-02, 1e-4);
        TS_ASSERT_DELTA(test_results[7] , 2.120643654085205e-01, 1e-4);
        TS_ASSERT_DELTA(test_results[8] , 1.439678172957377e+01, 1e-3);
        TS_ASSERT_DELTA(test_results[9] ,                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[10],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[11],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[12], 1.000000000000002e+01, 1e-4);
        TS_ASSERT_DELTA(test_results[13], 1.028341552112414e+02, 1e-2);
        TS_ASSERT_DELTA(test_results[14],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[15], 2.499999999999999e+01, 1e-4);
        TS_ASSERT_DELTA(test_results[16], 1.439678172957377e+01, 1e-3);
        TS_ASSERT_DELTA(test_results[17],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[18],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[19],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[20], 2.235636835087684e+00, 1e-4);
        TS_ASSERT_DELTA(test_results[21], 1.000000000000000e+00, 1e-4);

        // Acts as if it was divided at time = 16.1877... which is OK
        // (cell cycle model dictates division time, not when the cell is manually
        // divided)
        TissueCell daughter_cell = stem_cell.Divide();
        AbstractCellCycleModel* p_cell_model2 = daughter_cell.GetCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model->GetCurrentCellCyclePhase(), M_PHASE);
        TS_ASSERT_EQUALS(p_cell_model2->GetCurrentCellCyclePhase(), M_PHASE);

        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_cell_model2->ReadyToDivide(), false);

        // Check first 5 protein levels have been reset and the rest are the same.
        test_results = p_cell_model->GetProteinConcentrations();
        TS_ASSERT_DELTA(test_results[0] , 2.631865125420296e-01, 1e-3);
        TS_ASSERT_DELTA(test_results[1] , 2.678271949808561e-01, 1e-3);
        TS_ASSERT_DELTA(test_results[2] , 2.389956081120099e+00, 1e-3);
        TS_ASSERT_DELTA(test_results[3] , 1.390258620103223e+00, 1e-3);
        TS_ASSERT_DELTA(test_results[4] , 1.218603203963113e-01, 1e-3);
        TS_ASSERT_DELTA(test_results[5] , 1.428571428571429e-01, 1e-4);
        TS_ASSERT_DELTA(test_results[6] , 2.857142857142857e-02, 1e-4);
        TS_ASSERT_DELTA(test_results[7] , 2.120643654085205e-01, 1e-4);
        TS_ASSERT_DELTA(test_results[8] , 1.439678172957377e+01, 1e-3);
        TS_ASSERT_DELTA(test_results[9] ,                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[10],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[11],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[12], 1.000000000000002e+01, 1e-4);
        TS_ASSERT_DELTA(test_results[13], 1.028341552112414e+02, 1e-2);
        TS_ASSERT_DELTA(test_results[14],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[15], 2.499999999999999e+01, 1e-4);
        TS_ASSERT_DELTA(test_results[16], 1.439678172957377e+01, 1e-3);
        TS_ASSERT_DELTA(test_results[17],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[18],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[19],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[20], 2.235636835087684e+00, 1e-4);
        TS_ASSERT_DELTA(test_results[21], 1.000000000000000e+00, 1e-4);

        for (int i=0; i<9*num_timesteps/30.0; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
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
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetDimensionalisedTime(), 30.0, 1e-4);

        test_results = p_cell_model->GetProteinConcentrations();

        TS_ASSERT_DELTA(test_results[0] , 0.3636, 1e-3);
        TS_ASSERT_DELTA(test_results[1] , 0.9900, 1e-2);
        TS_ASSERT_DELTA(test_results[2] , 1.4996, 1e-3);
        TS_ASSERT_DELTA(test_results[3] , 0.8181, 1e-4);
        TS_ASSERT_DELTA(test_results[4] , 0.1000, 1e-4);
        TS_ASSERT_DELTA(test_results[5] , 0.6666, 1e-4);
        TS_ASSERT_DELTA(test_results[6] , 0.0666, 1e-4);
        TS_ASSERT_DELTA(test_results[7] , 0.7311, 1e-3);
        TS_ASSERT_DELTA(test_results[8] , 6.0481, 1e-3);
        TS_ASSERT_DELTA(test_results[9] , 0, 1e-4);
        TS_ASSERT_DELTA(test_results[10], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[11], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[12], 17.5432, 1e-4);
        TS_ASSERT_DELTA(test_results[13], 75.8316, 1e-2);
        TS_ASSERT_DELTA(test_results[14], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[15], 29.5728, 1e-4);
        TS_ASSERT_DELTA(test_results[16], 7.1564, 1e-3);
        TS_ASSERT_DELTA(test_results[17], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[18], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[19], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[20], 1.6048, 1e-4);
        TS_ASSERT_DELTA(test_results[21], 0.0000, 1e-4);

        TS_ASSERT_EQUALS(stem_cell.GetCellType(), DIFFERENTIATED);

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

        WntConcentration::Destroy();
    }


    void TestWntCellCycleModelForAPCSingleHit(void) throw(Exception)
    {
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        int num_timesteps = 500;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps); // 15.971 hours to go into S phase

        double wnt_level = 1.0;
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

        WntCellCycleModel* p_cell_model = new WntCellCycleModel();

        TissueCell stem_cell(STEM, HEALTHY, p_cell_model);

        stem_cell.InitialiseCellCycleModel();

        double SG2MDuration = CancerParameters::Instance()->GetSG2MDuration();
        TS_ASSERT_THROWS_NOTHING(WntCellCycleModel cell_model_3());

        WntCellCycleModel* p_cell_model_1 = new WntCellCycleModel();

        TissueCell stem_cell_1(STEM, APC_ONE_HIT, p_cell_model_1);
        stem_cell_1.InitialiseCellCycleModel();

        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 4.804 hrs and then finish dividing
        // 10 hours later at 14.804 hours.
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_1, 4.804);
        }

        p_cell_model_1->ResetForDivision();
        double second_cycle_start = p_cell_model_1->GetBirthTime();

        TS_ASSERT_DELTA(SG2MDuration, 10.0, 1e-5);
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            bool result = p_cell_model_1->ReadyToDivide();
            if (time< second_cycle_start + 4.804 + SG2MDuration)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }

        WntConcentration::Destroy();
    }


    void TestWntCellCycleModelForBetaCatSingleHit(void) throw(Exception)
    {
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        int num_timesteps = 500;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps); // 15.971 hours to go into S phase

        double wnt_level = 0.0;
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

        WntCellCycleModel* p_cell_model = new WntCellCycleModel();

        TissueCell stem_cell(STEM, BETA_CATENIN_ONE_HIT, p_cell_model);
        stem_cell.InitialiseCellCycleModel();

        TS_ASSERT_THROWS_NOTHING(WntCellCycleModel cell_model_3());

        WntCellCycleModel* p_cell_model_1 = new WntCellCycleModel();

        TissueCell stem_cell_1(STEM, BETA_CATENIN_ONE_HIT, p_cell_model_1);
        stem_cell_1.InitialiseCellCycleModel();

        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 7.82 hrs and then finish dividing
        // 10 hours later at 17.82 hours.
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_THROWS_ANYTHING(p_cell_model_1->UpdateCellType());
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_1, 7.82);
        }

        p_cell_model_1->ResetForDivision();

        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_1, 7.82);
        }

        WntConcentration::Destroy();
    }


    void TestWntCellCycleModelForAPCDoubleHit(void) throw(Exception)
    {
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        int num_timesteps = 500;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps); // 15.971 hours to go into S phase

        double wnt_level = 0.738; // This shouldn't matter for this kind of cell!
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

        WntCellCycleModel* p_cell_model_1 = new WntCellCycleModel();

        TissueCell stem_cell_1(STEM, APC_TWO_HIT, p_cell_model_1);
        stem_cell_1.InitialiseCellCycleModel();

        WntCellCycleModel* p_cell_model_2 = new WntCellCycleModel();

        TissueCell stem_cell_2(STEM, APC_TWO_HIT, p_cell_model_2);
        stem_cell_2.InitialiseCellCycleModel();

        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 3.943 hrs and then finish dividing
        // 10 hours later at 13.9435 hours.
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_2, 3.9435);
        }

        p_cell_model_2->ResetForDivision();

        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_2, 3.9435);
        }

        WntConcentration::Destroy();
    }


    void TestWntCellCycleModelForConstantWntStimulusHealthyCell(void) throw(Exception)
    {
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        int num_timesteps = 500;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase

        double wnt_level = 1.0;
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

        WntCellCycleModel* p_cell_model_1 = new WntCellCycleModel();

        TissueCell stem_cell_1(STEM, HEALTHY, p_cell_model_1);
        stem_cell_1.InitialiseCellCycleModel();

        WntCellCycleModel* p_cell_model_2 = new WntCellCycleModel();

        TissueCell stem_cell_2(STEM, HEALTHY, p_cell_model_2);
        stem_cell_2.InitialiseCellCycleModel();

        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 5.971 hrs and then finish dividing
        // 10 hours later at 15.971 hours.
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_2, 5.971);
        }

        p_cell_model_2->ResetForDivision();

        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cell_model_2, 5.971);
        }

        WntConcentration::Destroy();
    }


    void TestStochasticWntCellCycleModel() throw (Exception)
    {
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        int num_timesteps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(20, num_timesteps);// 15.971 hours to go into S phase

        double wnt_level = 1.0;
        WntConcentration::Instance()->SetConstantWntValueForTesting(wnt_level);

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

        for (int i=0; i<num_timesteps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            bool result = p_cell_model->ReadyToDivide();

            if (time < 5.971 + 9.0676)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }

        WntConcentration::Destroy();
    }


    void TestAlarcon2004OxygenBasedCellCycleModel() throw(Exception)
    {
        CancerParameters::Instance()->SetHepaOneParameters();

        // Set up SimulationTime
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(20.0, 2);

        // Set up oxygen_concentration
        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        // Create model
        Alarcon2004OxygenBasedCellCycleModel* p_cell_model = new Alarcon2004OxygenBasedCellCycleModel();

        // Create cell
        TissueCell cell(STEM, HEALTHY, p_cell_model);

        // Coverage of cell cycle model copying without an ODE system set up
        TissueCell stem_cell2 = cell;
        TS_ASSERT_EQUALS(stem_cell2.GetMutationState(), HEALTHY);

        cell.InitialiseCellCycleModel();

        // Check oxygen concentration is correct in cell cycle model
        TS_ASSERT_DELTA(p_cell_model->GetProteinConcentrations()[5], 1.0, 1e-5);
        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), false);

        // Divide a cell
        Alarcon2004OxygenBasedCellCycleModel *p_cell_model2 = static_cast <Alarcon2004OxygenBasedCellCycleModel*> (p_cell_model->CreateCellCycleModel());

        TissueCell cell2(STEM, HEALTHY, p_cell_model2);

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),false)
        TS_ASSERT_EQUALS(p_cell_model2->ReadyToDivide(),false);

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),true)
        TS_ASSERT_EQUALS(p_cell_model2->ReadyToDivide(),true);

        TS_ASSERT_THROWS_NOTHING(p_cell_model->ResetForDivision());

        CellwiseData<2>::Destroy();
    }


    void TestArchiveTysonNovakCellCycleModels()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "tyson_novak_cell_cycle.arch";

        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(100.0, 1);

            TysonNovakCellCycleModel* p_model = new TysonNovakCellCycleModel;

            p_simulation_time->IncrementTimeOneStep();

            TissueCell cell(TRANSIT, HEALTHY, p_model);            cell.InitialiseCellCycleModel();

            TS_ASSERT_EQUALS(p_model->ReadyToDivide(),true);

            p_model->SetBirthTime(-1.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            TissueCell* const p_cell = &cell;

            output_arch << p_cell;

            SimulationTime::Destroy();
        }

        {
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

            // Check
            TS_ASSERT_EQUALS(p_model->ReadyToDivide(),true);
            TS_ASSERT_DELTA(p_model->GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(),101.0,1e-12);
            delete p_cell;
        }
    }


    void TestArchiveWntCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "wnt_cell_cycle.arch";
        WntConcentration::Instance()->SetConstantWntValueForTesting(1.0);

        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(16, 2);

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

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            TissueCell* const p_cell = &stem_cell;

            output_arch << p_cell;
            SimulationTime::Destroy();
        }

        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            CancerParameters *inst1 = CancerParameters::Instance();

            inst1->SetSDuration(101.0);

            TissueCell* p_cell;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_cell;

            // Check
            AbstractCellCycleModel* p_cell_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());

            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),true);
            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_DELTA(p_cell_model->GetAge(),17.0,1e-12);
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
            TS_ASSERT_EQUALS(p_cell_model->GetCurrentCellCyclePhase(), G_TWO_PHASE);
            delete p_cell;
        }

        WntConcentration::Destroy();
    }


    void TestArchiveIngeWntSwatCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "inge_wnt_swat_cell_cycle.arch";
        WntConcentration::Instance()->SetConstantWntValueForTesting(1.0);

        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(17, 2);

            IngeWntSwatCellCycleModel* p_cell_model = new IngeWntSwatCellCycleModel(1);

            TissueCell stem_cell(STEM, HEALTHY, p_cell_model);
            stem_cell.InitialiseCellCycleModel();

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(),false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(),true);

            stem_cell.GetCellCycleModel()->SetBirthTime(-1.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            TissueCell* const p_cell = &stem_cell;

            output_arch << p_cell;
            SimulationTime::Destroy();
        }

        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            CancerParameters *inst1 = CancerParameters::Instance();

            inst1->SetSDuration(101.0);

            TissueCell* p_cell;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_cell;

            // Check
            AbstractCellCycleModel* p_cell_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());

            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),true);
            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_DELTA(p_cell_model->GetAge(),18.0,1e-12);
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
            TS_ASSERT_EQUALS(p_cell_model->GetCurrentCellCyclePhase(), G_TWO_PHASE);

            delete p_cell;
        }

        WntConcentration::Destroy();
    }


    void TestArchiveStochasticWntCellCycleModels()
    {
        // In this case the first cycle time will be 5.971+9.0676 = 15.0386
        // note that the S-G2-M time is assigned when the cell finishes G1
        //(i.e. at time 5.971 here so the model has to be archived BEFORE that.
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "stochastic_wnt_cell_cycle.arch";
        WntConcentration::Instance()->SetConstantWntValueForTesting(1.0);

        // Create an ouput archive
        {   // In this test the RandomNumberGenerator in existence
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(16.0, 1000);

            StochasticWntCellCycleModel* p_stoc_model = new StochasticWntCellCycleModel();

            TissueCell stoc_cell(STEM, HEALTHY, p_stoc_model);
            stoc_cell.InitialiseCellCycleModel();

            WntCellCycleModel* p_wnt_model = new WntCellCycleModel();
            TissueCell wnt_cell(STEM, HEALTHY, p_wnt_model);
            wnt_cell.InitialiseCellCycleModel();

            p_simulation_time->IncrementTimeOneStep(); // 5.5

            while (p_simulation_time->GetDimensionalisedTime() < 4.0)
            {
                p_simulation_time->IncrementTimeOneStep();
            }

            TS_ASSERT_EQUALS(stoc_cell.GetCellCycleModel()->ReadyToDivide(),false);
            TS_ASSERT_EQUALS(wnt_cell.GetCellCycleModel()->ReadyToDivide(),false);
            TS_ASSERT_EQUALS(stoc_cell.GetCellCycleModel()->GetCurrentCellCyclePhase(), G_ONE_PHASE);

            // When these are included here they pass - so are moved down into
            // after load to see if they still pass.

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            TissueCell* const p_wnt_cell = &wnt_cell;
            TissueCell* const p_stoc_cell = &stoc_cell;

            output_arch << p_stoc_cell;
            output_arch << p_wnt_cell;
            SimulationTime::Destroy();
        }

        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(16.0, 2);

            CancerParameters *inst1 = CancerParameters::Instance();

            inst1->SetSDuration(101.0);

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

            TS_ASSERT_EQUALS(p_stoc_cell->GetCellCycleModel()->GetCurrentCellCyclePhase(), G_ONE_PHASE);

            // Check - stochastic should divide at 15.03
            // Wnt should divide at 15.971
            while (p_simulation_time->GetDimensionalisedTime() < 15.0)
            {
                p_simulation_time->IncrementTimeOneStep();
            }
            TS_ASSERT_EQUALS(p_stoc_cell->GetCellCycleModel()->ReadyToDivide(),false);
            TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->ReadyToDivide(),false);

            while (p_simulation_time->GetDimensionalisedTime() < 15.5)
            {
                p_simulation_time->IncrementTimeOneStep();
            }
            TS_ASSERT_EQUALS(p_stoc_cell->GetCellCycleModel()->ReadyToDivide(),true);// only for stochastic
            TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->ReadyToDivide(),false);

            while (p_simulation_time->GetDimensionalisedTime() < 16.0)
            {
                p_simulation_time->IncrementTimeOneStep();
            }
            TS_ASSERT_EQUALS(p_stoc_cell->GetCellCycleModel()->ReadyToDivide(),true);
            TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->ReadyToDivide(),true);

            TS_ASSERT_DELTA(p_stoc_cell->GetCellCycleModel()->GetBirthTime(),0.0,1e-12);
            TS_ASSERT_DELTA(p_stoc_cell->GetCellCycleModel()->GetAge(),16.0,1e-12);
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);

            delete p_stoc_cell;
            delete p_wnt_cell;
        }
        WntConcentration::Destroy();
    }


    void TestArchiveAlarcon2004OxygenBasedCellCycleModels()
    {
        CancerParameters::Instance()->SetHepaOneParameters();

        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "alarcon_cell_cycle.arch";

        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 2);

            Alarcon2004OxygenBasedCellCycleModel* p_cell_model = new Alarcon2004OxygenBasedCellCycleModel();

            TissueCell cell(STEM, HEALTHY, p_cell_model);
            cell.InitialiseCellCycleModel();
            cell.GetCellCycleModel()->SetBirthTime(-10.0);

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(cell.GetCellCycleModel()->ReadyToDivide(),false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(cell.GetCellCycleModel()->ReadyToDivide(),true);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            TissueCell* const p_cell = &cell;

            output_arch << p_cell;
            SimulationTime::Destroy();
        }

        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            CancerParameters *inst1 = CancerParameters::Instance();

            inst1->SetSDuration(101.0);

            TissueCell* p_cell;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_cell;

            // Check
            AbstractCellCycleModel* p_cell_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());

            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),true);
            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(),-10.0,1e-12);
            TS_ASSERT_DELTA(p_cell_model->GetAge(),20.0,1e-12);
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
            delete p_cell;
        }

        CellwiseData<2>::Destroy();
    }

};

#endif /*TESTCELLCYCLEMODELSODE_HPP_*/
