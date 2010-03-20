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
#ifndef TESTSINGLEODEWNTCELLCYCLEMODEL_HPP_
#define TESTSINGLEODEWNTCELLCYCLEMODEL_HPP_

#include <cxxtest/TestSuite.h>
#include <fstream>

// Must be included before any other cell_based headers
#include "TissueSimulationArchiver.hpp"
#include "SingleOdeWntCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "LabelledCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"

class TestSingleOdeWntCellCycleModel : public AbstractCellBasedTestSuite
{
public:

    void TestCorrectBehaviour() throw(Exception)
    {
        // Set up the simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double end_time = 1200.0;
        unsigned num_timesteps = 10*(unsigned)end_time;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

        // Set up Wnt concentration
        double wnt_level = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Create cell mutation state
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

        // Create cell cycle model
        SingleOdeWntCellCycleModel* p_cycle_model = new SingleOdeWntCellCycleModel();
        p_cycle_model->SetDimension(2);

        // Construct a cell with this cell cycle model and cell mutation state
        TissueCell cell(STEM, p_state, p_cycle_model);
        cell.InitialiseCellCycleModel();

        // Test the cell cycle model is behaving correctly
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // Stem cell should have been changed into a transit cell by wnt cell cycle model
            TS_ASSERT_EQUALS(cell.GetCellProliferativeType(), TRANSIT);

            // The number for the G1 duration is taken from the first random number generated
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, 1.0676);
        }

        double steady_beta_cat_at_wnt_equals_1 = p_cycle_model->GetBetaCateninConcentration();

#ifdef CHASTE_CVODE
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_1, 143.4962, 1e-4);
#else
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_1, 143.5119, 1e-4);
#endif

        // Divide the cell
        TS_ASSERT_EQUALS(cell.ReadyToDivide(), true);
        TissueCell cell2 = cell.Divide();
        boost::shared_ptr<AbstractCellMutationState> p_labelled(new LabelledCellMutationState);
        cell.SetMutationState(p_labelled);

        SingleOdeWntCellCycleModel* p_cycle_model2 = static_cast<SingleOdeWntCellCycleModel*> (cell2.GetCellCycleModel());

        // Test that both cells have inherited the same Beta-catenin concentration.
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_1, p_cycle_model->GetBetaCateninConcentration(), 1e-12);
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_1, p_cycle_model2->GetBetaCateninConcentration(), 1e-12);

        // Now reduce the Wnt concentration
        wnt_level = 0.2;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        double division_time = SimulationTime::Instance()->GetTime();

        // The numbers for the G1 durations are taken from
        // the first two random numbers generated
        double new_g1_duration = 3.16316;
        double new_g1_duration2 = 1.2712;
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, new_g1_duration);
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model2, new_g1_duration2);
        }

        TS_ASSERT_EQUALS(cell.GetCellProliferativeType(), TRANSIT);
        TS_ASSERT_EQUALS(cell2.GetCellProliferativeType(), TRANSIT);

        // Test that both cells have inherited the same beta-catenin concentration
        double steady_beta_cat_at_wnt_equals_0_2 = 113.4683;
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_0_2, p_cycle_model->GetBetaCateninConcentration(), 1e-3);
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_0_2, p_cycle_model2->GetBetaCateninConcentration(), 1e-3);

        p_cycle_model->ResetForDivision();
        p_cycle_model2->ResetForDivision();

        // Test that both cells have still inherited the same ceta-catenin concentration
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_0_2, p_cycle_model->GetBetaCateninConcentration(), 1e-3);
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_0_2, p_cycle_model2->GetBetaCateninConcentration(), 1e-3);

        division_time = SimulationTime::Instance()->GetTime();

        // Now reduce the Wnt concentration so only mutant cells divide
        wnt_level = 0.1;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Introduce a mutation (no immediate effect)
        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);
        cell2.SetMutationState(p_apc1);
        TS_ASSERT(!p_cycle_model->ReadyToDivide());
        TS_ASSERT(!p_cycle_model2->ReadyToDivide());

        // Coverage
        cell2.SetMutationState(p_apc2);
        TS_ASSERT(!p_cycle_model->ReadyToDivide());

        cell2.SetMutationState(p_bcat1);
        TS_ASSERT(!p_cycle_model->ReadyToDivide());

        // The numbers for the G1 durations are taken from the next two random numbers generated
        new_g1_duration = 1.22037;
        new_g1_duration2 = 0.74699;

        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, new_g1_duration);
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model2, new_g1_duration2);
        }

        TS_ASSERT_DELTA(91.6693, p_cycle_model->GetBetaCateninConcentration(), 1e-3);
#ifdef CHASTE_CVODE
        TS_ASSERT_DELTA(355.9138, p_cycle_model2->GetBetaCateninConcentration(), 1e-3);
#else
        TS_ASSERT_DELTA(355.7790, p_cycle_model2->GetBetaCateninConcentration(), 1e-3);
#endif

        TS_ASSERT_DELTA( p_cycle_model->GetBetaCateninDivisionThreshold(), 100, 1e-9);
        TS_ASSERT_DELTA(p_cycle_model2->GetBetaCateninDivisionThreshold(), 100, 1e-9);

        TS_ASSERT_EQUALS(cell.GetCellProliferativeType(), DIFFERENTIATED);
        TS_ASSERT_EQUALS(cell2.GetCellProliferativeType(), TRANSIT);

        // Coverage of 1D

        // Set up SimulationTime
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(30.0, 2);

        // Instantiate 1D Wnt concentration
        wnt_level = 1.0;
        WntConcentration<1>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Create cell cycle model and cell and test behaviour
        SingleOdeWntCellCycleModel* p_cell_model_1d = new SingleOdeWntCellCycleModel;
        p_cell_model_1d->SetDimension(1);
        p_cell_model_1d->SetUseCellProliferativeTypeDependentG1Duration();

        TS_ASSERT_EQUALS(p_cell_model_1d->GetDimension(), 1u);

        TissueCell stem_cell_1d(STEM, p_state, p_cell_model_1d);
        stem_cell_1d.InitialiseCellCycleModel();

        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(stem_cell_1d.ReadyToDivide(), false);

        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(stem_cell_1d.ReadyToDivide(), true);

        TissueCell daughter_1d = stem_cell_1d.Divide();

        // Coverage of 3D

        // Set up SimulationTime
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(20.0, 2);

        // Instantiate 3D Wnt concentration
        WntConcentration<3>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Create cell cycle model and cell and test behaviour
        SingleOdeWntCellCycleModel* p_cell_model_3d = new SingleOdeWntCellCycleModel;
        p_cell_model_3d->SetDimension(3);

        TS_ASSERT_EQUALS(p_cell_model_3d->GetDimension(), 3u);

        TissueCell stem_cell_3d(STEM, p_state, p_cell_model_3d);
        stem_cell_3d.InitialiseCellCycleModel();

        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(stem_cell_3d.ReadyToDivide(), false);

        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(stem_cell_3d.ReadyToDivide(), true);

        TissueCell daughter_3d = stem_cell_3d.Divide();

        // Tidy up
        WntConcentration<1>::Destroy();
        WntConcentration<2>::Destroy();
        WntConcentration<3>::Destroy();
    }

    // This test currently fails as archiving is not yet correctly implemented for SingleOdeWntCellCycleModel
    // (see #1239)
    void DONTTestArchiving()
    {
        TissueConfig* p_params = TissueConfig::Instance();

        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "single_ode_wnt.arch";

        double random_number_test = 0;

        // Create an output archive
        {
            // The number for the G1 duration is taken from the first random number generated
            double g1_duration = 1.0676;

            // Set up the simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            double end_time = g1_duration + p_params->GetSG2MDuration() + 5.0;
            unsigned num_timesteps = 50;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

            // Set up the Wnt concentration for testing
            WntConcentration<1>::Instance()->SetConstantWntValueForTesting(0.7);

            // Create cell cycle model and associated cell
            SimpleWntCellCycleModel* p_cell_model = new SingleOdeWntCellCycleModel;
            p_cell_model->SetDimension(1);
            p_cell_model->SetBirthTime(-1.0);

            boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

            TissueCell stem_cell(STEM, p_healthy_state, p_cell_model);
            stem_cell.InitialiseCellCycleModel();

            while (p_cell_model->GetAge() < g1_duration + p_params->GetSG2MDuration()
                    - p_simulation_time->GetTimeStep()) // minus one to match birth time.
            {
                p_simulation_time->IncrementTimeOneStep();
                CheckReadyToDivideAndPhaseIsUpdated(p_cell_model, g1_duration);
            }

            // Wnt should change this to a transit cell
            TS_ASSERT_EQUALS(stem_cell.GetCellProliferativeType(), TRANSIT);
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(), false);
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->GetCurrentCellCyclePhase(), G_TWO_PHASE);

            // Archive the cell
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            TissueCell* const p_cell = &stem_cell;
            output_arch << p_cell;

            TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT(stem_cell.GetCellCycleModel()->ReadyToDivide());

            // Tidy up
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            random_number_test = p_gen->ranf();
            SimulationTime::Destroy();
        }

        {
            // Set up SimulationTime
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            TissueConfig* p_inst1 = TissueConfig::Instance();
            p_inst1->SetSDuration(101.0);

            TissueCell* p_cell;

            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            p_gen->Reseed(36);

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_cell;

            // Check
            AbstractCellCycleModel* p_cell_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());

            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), true);

            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetSG2MDuration(), 10.0, 1e-12);

            TS_ASSERT_DELTA(p_gen->ranf(), random_number_test, 1e-7);
            TS_ASSERT_EQUALS((static_cast<SingleOdeWntCellCycleModel*>(p_cell_model))->GetDimension(), 1u);

            // Tidy up
            delete p_cell;
            SimulationTime::Destroy();
        }
    }
};


#endif /* TESTSINGLEODEWNTCELLCYCLEMODEL_HPP_ */
