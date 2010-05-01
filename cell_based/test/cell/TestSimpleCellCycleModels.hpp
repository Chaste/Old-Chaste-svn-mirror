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
#ifndef TESTSIMPLECELLCYCLEMODELS_HPP_
#define TESTSIMPLECELLCYCLEMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "StochasticOxygenBasedCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "LabelledCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"

class TestSimpleCellCycleModels : public AbstractCellBasedTestSuite
{
public:

    void TestFixedDurationGenerationBasedCellCycleModel() throw(Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            2.0*(p_params->GetStemCellG1Duration()
                  +p_params->GetSG2MDuration()     ), num_steps);

        TS_ASSERT_THROWS_NOTHING(FixedDurationGenerationBasedCellCycleModel model3);

        FixedDurationGenerationBasedCellCycleModel* p_stem_model = new FixedDurationGenerationBasedCellCycleModel;
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        TissueCell stem_cell(STEM, p_healthy_state, p_stem_model);
        stem_cell.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(),M_PHASE);
        TS_ASSERT_EQUALS(p_stem_model->GetGeneration(), 0u);
        TS_ASSERT_EQUALS(p_stem_model->CanCellTerminallyDifferentiate(), true);

        TS_ASSERT_EQUALS(stem_cell.GetCellProliferativeType(),STEM);

        FixedDurationGenerationBasedCellCycleModel* p_transit_model = new FixedDurationGenerationBasedCellCycleModel;
        TissueCell transit_cell(TRANSIT, p_healthy_state, p_transit_model);
        transit_cell.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(transit_cell.GetCellProliferativeType(),TRANSIT);
        TS_ASSERT_EQUALS(p_transit_model->GetGeneration(), 0u);

        FixedDurationGenerationBasedCellCycleModel* p_diff_model = new FixedDurationGenerationBasedCellCycleModel;
        TissueCell diff_cell(DIFFERENTIATED, p_healthy_state, p_diff_model);
        diff_cell.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(diff_cell.GetCellProliferativeType(),DIFFERENTIATED);
        TS_ASSERT_EQUALS(p_diff_model->GetGeneration(), 0u);

        for (unsigned i=0; i<num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, p_params->GetStemCellG1Duration());
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, p_params->GetTransitCellG1Duration());
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 100);
        }

        TS_ASSERT_DELTA(p_stem_model->GetAge(), p_simulation_time->GetTime(), 1e-9);
        TS_ASSERT_DELTA(p_transit_model->GetAge(), p_simulation_time->GetTime(), 1e-9);
        TS_ASSERT_DELTA(p_diff_model->GetAge(), p_simulation_time->GetTime(), 1e-9);

        double hepa_one_cell_birth_time = p_simulation_time->GetTime();

        p_params->SetHepaOneParameters();
        FixedDurationGenerationBasedCellCycleModel* p_hepa_one_model = new FixedDurationGenerationBasedCellCycleModel;
        TissueCell hepa_one_cell(STEM, p_healthy_state, p_hepa_one_model);
        hepa_one_cell.InitialiseCellCycleModel();

        for (unsigned i=0; i<num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, p_params->GetHepaOneCellG1Duration());
        }

        TS_ASSERT_DELTA(p_hepa_one_model->GetAge() + hepa_one_cell_birth_time, p_simulation_time->GetTime(), 1e-9);
    }


    void TestStochasticDurationGenerationBasedCellCycleModel() throw(Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            2.0*(p_params->GetStemCellG1Duration() + p_params->GetSG2MDuration()), num_steps);

        TS_ASSERT_THROWS_NOTHING(StochasticDurationGenerationBasedCellCycleModel cell_model3);

        StochasticDurationGenerationBasedCellCycleModel* p_stem_model = new StochasticDurationGenerationBasedCellCycleModel;
        StochasticDurationGenerationBasedCellCycleModel* p_transit_model = new StochasticDurationGenerationBasedCellCycleModel;
        StochasticDurationGenerationBasedCellCycleModel* p_diff_model = new StochasticDurationGenerationBasedCellCycleModel;

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        TissueCell stem_cell(STEM, p_healthy_state,  p_stem_model);
        stem_cell.InitialiseCellCycleModel();

        TissueCell transit_cell(TRANSIT, p_healthy_state, p_transit_model);
        transit_cell.InitialiseCellCycleModel();

        TissueCell diff_cell(DIFFERENTIATED, p_healthy_state, p_diff_model);
        diff_cell.InitialiseCellCycleModel();

        for (unsigned i=0; i<num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The numbers for the G1 durations below are taken from the first three
            // random numbers generated
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, 4.36075);
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, 1.78877);
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132);  // any old number
        }

        p_params->SetHepaOneParameters();
        StochasticDurationGenerationBasedCellCycleModel* p_hepa_one_model = new StochasticDurationGenerationBasedCellCycleModel;
        TissueCell hepa_one_cell(STEM, p_healthy_state, p_hepa_one_model);
        hepa_one_cell.InitialiseCellCycleModel();

        for (unsigned i=0; i< num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, 4.1324);
        }
    }


    void TestSimpleWntCellCycleModel() throw(Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();

        // Set up the simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        double end_time = 60.0;
        unsigned num_timesteps = 1000*(unsigned)end_time;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

        // Set up the Wnt concentration
        double wnt_level = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        SimpleWntCellCycleModel* p_cycle_model = new SimpleWntCellCycleModel;
        TS_ASSERT_EQUALS(p_cycle_model->GetDimension(), 0u);
        TS_ASSERT_EQUALS(p_cycle_model->CanCellTerminallyDifferentiate(), false);

        // Test the dimension must be 1, 2 or 3
        TS_ASSERT_THROWS_THIS(p_cycle_model->SetDimension(4), "Dimension must be 1, 2 or 3");

        // Test the set/get dimension methods
        p_cycle_model->SetDimension(2);
        TS_ASSERT_EQUALS(p_cycle_model->GetDimension(), 2u);

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        TissueCell cell(STEM, p_healthy_state, p_cycle_model);
        cell.InitialiseCellCycleModel();

        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The number for the G1 duration is taken from
            // the first random number generated
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, 1.0676);
        }
        // Stem cell should have been changed into a transit cell by wnt cell cycle model
        TS_ASSERT_EQUALS(cell.GetCellProliferativeType(), TRANSIT);

        // Divide the cell
        TS_ASSERT_EQUALS(cell.ReadyToDivide(), true);
        TissueCell cell2 = cell.Divide();
        boost::shared_ptr<AbstractCellMutationState> p_lab_mutation(new LabelledCellMutationState);
        cell.SetMutationState(p_lab_mutation);

        SimpleWntCellCycleModel* p_cycle_model2 = static_cast<SimpleWntCellCycleModel*> (cell2.GetCellCycleModel());

        // Now reduce the Wnt concentration
        wnt_level = 0.7;
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

        p_cycle_model->ResetForDivision();
        p_cycle_model2->ResetForDivision();

        division_time = SimulationTime::Instance()->GetTime();

        // Now reduce the Wnt concentration so only beta-cat or APC2 hit cells divide.
        wnt_level = 0.15;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        boost::shared_ptr<AbstractCellMutationState> p_apc1_mutation(new ApcOneHitCellMutationState);
        cell.SetMutationState(p_apc1_mutation);
        boost::shared_ptr<AbstractCellMutationState> p_bcat_mutation(new BetaCateninOneHitCellMutationState);
        cell2.SetMutationState(p_bcat_mutation);

        TS_ASSERT(!p_cycle_model->ReadyToDivide());
        TS_ASSERT(!p_cycle_model2->ReadyToDivide());

        // Coverage...
        boost::shared_ptr<AbstractCellMutationState> p_apc2_mutation(new ApcTwoHitCellMutationState);
        cell.SetMutationState(p_apc2_mutation);
        TS_ASSERT(!p_cycle_model->ReadyToDivide());
        cell.SetMutationState(p_apc1_mutation);
        TS_ASSERT(!p_cycle_model->ReadyToDivide());

        // The numbers for the G1 durations are taken from
        // the first two random numbers generated
        new_g1_duration = 1.22037;
        new_g1_duration2 = 0.74699;

        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, new_g1_duration);
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model2, new_g1_duration2);
        }

        TS_ASSERT_EQUALS(cell.GetCellProliferativeType(), DIFFERENTIATED);
        TS_ASSERT_EQUALS(cell2.GetCellProliferativeType(), TRANSIT);

        // For coverage...
        SimpleWntCellCycleModel* p_cycle_model1 = new SimpleWntCellCycleModel;
        p_cycle_model1->SetDimension(2);

        TissueCell cell1(DIFFERENTIATED, p_healthy_state, p_cycle_model1);
        cell1.InitialiseCellCycleModel();

        SimpleWntCellCycleModel* p_another_cycle_model = new SimpleWntCellCycleModel;
        p_another_cycle_model->SetDimension(2);

        TissueCell another_cell(STEM, p_healthy_state, p_another_cycle_model);
        another_cell.InitialiseCellCycleModel();
        // ...end of coverage

        /*
         * Test the case of a radial Wnt concentration
         */

        p_params->Reset();
        RandomNumberGenerator::Instance()->Reseed(0);

        // Set up the Wnt concentration
        wnt_level = p_params->GetWntStemThreshold() + 0.01;
        WntConcentration<2>::Destroy();
        WntConcentration<2>::Instance()->SetType(RADIAL);
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Set up a cell cycle model and cell
        SimpleWntCellCycleModel* p_cycle_model4 = new SimpleWntCellCycleModel;
        p_cycle_model4->SetDimension(2);

        TissueCell cell4(STEM, p_healthy_state,  p_cycle_model4);
        cell4.InitialiseCellCycleModel();

        // Test the GetCurrentCellCyclePhase() and ReadyToDivide() methods
        double first_g1_duration = 1.0676;
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The number for the G1 duration is taken from
            // the first random number generated
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model4, first_g1_duration);
        }

        // We should still have a stem cell since the WntConcentration exceeds mRadialWntThreshold
        TS_ASSERT_EQUALS(cell4.GetCellProliferativeType(), STEM);

        // Divide the cell
        TS_ASSERT_EQUALS(cell4.ReadyToDivide(), true);
        TS_ASSERT_EQUALS(cell4.GetCellProliferativeType(), STEM);
        TissueCell cell5 = cell4.Divide();
        TS_ASSERT_EQUALS(cell4.GetCellProliferativeType(), STEM);
        TS_ASSERT_EQUALS(cell5.GetCellProliferativeType(), TRANSIT);
        cell2.SetMutationState(p_lab_mutation);

        // Now reduce the Wnt concentration
        wnt_level = p_params->GetWntStemThreshold() - 0.01;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // The numbers for the G1 durations are taken from
        // the first two random numbers generated
        new_g1_duration = 3.16316;
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model4, new_g1_duration);
        }

        TS_ASSERT_DELTA(WntConcentration<2>::Instance()->GetWntLevel(cell4), wnt_level, 1e-12);
        TS_ASSERT_EQUALS(cell4.GetCellProliferativeType(), TRANSIT);
        TS_ASSERT_EQUALS(cell5.GetCellProliferativeType(), TRANSIT);


        // Coverage of 1D

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(30.0, 2);

        WntConcentration<1>::Instance()->SetConstantWntValueForTesting(wnt_level);
        SimpleWntCellCycleModel* p_cell_model_1d = new SimpleWntCellCycleModel;
        p_cell_model_1d->SetDimension(1);
        p_cell_model_1d->SetUseCellProliferativeTypeDependentG1Duration();

        TS_ASSERT_EQUALS(p_cell_model_1d->GetDimension(), 1u);

        TissueCell stem_cell_1d(STEM, p_healthy_state, p_cell_model_1d);
        stem_cell_1d.InitialiseCellCycleModel();

        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(stem_cell_1d.ReadyToDivide(), false);

        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(stem_cell_1d.ReadyToDivide(), true);

        TissueCell daughter_1d = stem_cell_1d.Divide();

        // Coverage of 3D

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(20.0, 2);

        WntConcentration<3>::Instance()->SetConstantWntValueForTesting(wnt_level);
        SimpleWntCellCycleModel* p_cell_model_3d = new SimpleWntCellCycleModel;
        p_cell_model_3d->SetDimension(3);

        TS_ASSERT_EQUALS(p_cell_model_3d->GetDimension(), 3u);

        TissueCell stem_cell_3d(STEM, p_healthy_state, p_cell_model_3d);
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


    void TestSimpleOxygenBasedCellCycleModel() throw(Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();
        p_params->SetHepaOneParameters();

        // Check that mCurrentHypoxiaOnsetTime and mCurrentHypoxicDuration are
        // updated correctly
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

        SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
        p_model->SetDimension(2);
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        TissueCell cell(STEM, p_state, p_model);

        cell.InitialiseCellCycleModel();

        // Set up constant oxygen_concentration
        std::vector<double> low_oxygen_concentration;
        std::vector<double> high_oxygen_concentration;
        low_oxygen_concentration.push_back(0.0);
        high_oxygen_concentration.push_back(1.0);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);

        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);

        p_simulation_time->IncrementTimeOneStep(); // t=1.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(high_oxygen_concentration);

        p_simulation_time->IncrementTimeOneStep(); // t=2.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 2.0, 1e-12);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);
        p_simulation_time->IncrementTimeOneStep(); // t=3.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 2.0, 1e-12);

        // Set up SimulationTime
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            4.0*(p_params->GetHepaOneCellG1Duration()
                  +p_params->GetSG2MDuration()     ), num_steps);

        // Set up constant oxygen concentration
        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        TS_ASSERT_THROWS_NOTHING(SimpleOxygenBasedCellCycleModel model);

        // Create cell cycle model
        SimpleOxygenBasedCellCycleModel* p_hepa_one_model = new SimpleOxygenBasedCellCycleModel;
        SimpleOxygenBasedCellCycleModel* p_diff_model = new SimpleOxygenBasedCellCycleModel;
        p_hepa_one_model->SetDimension(2);
        p_diff_model->SetDimension(2);

        // Create cell
        TissueCell hepa_one_cell(STEM, p_state, p_hepa_one_model);
        hepa_one_cell.InitialiseCellCycleModel();

        TissueCell diff_cell(DIFFERENTIATED, p_state, p_diff_model);
        diff_cell.InitialiseCellCycleModel();

        // Check that the cell cycle phase and ready to divide
        // are updated correctly
        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model->GetCurrentCellCyclePhase(),M_PHASE);

        TS_ASSERT_EQUALS(p_diff_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_diff_model->GetCurrentCellCyclePhase(),G_ZERO_PHASE);

        for (unsigned i=0; i<num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // Note that we need to pass in the updated G1 duration
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, p_hepa_one_model->GetG1Duration());
        }

        TS_ASSERT_DELTA(p_hepa_one_model->GetAge(), p_simulation_time->GetTime(), 1e-9);
        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(), true);

        // Check that cell division correctly resets the cell cycle phase
        TS_ASSERT_EQUALS(hepa_one_cell.ReadyToDivide(), true);
        TissueCell hepa_one_cell2 = hepa_one_cell.Divide();
        SimpleOxygenBasedCellCycleModel* p_hepa_one_model2 = static_cast <SimpleOxygenBasedCellCycleModel*>(hepa_one_cell2.GetCellCycleModel());

        TS_ASSERT_EQUALS(p_hepa_one_model2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model2->GetCurrentCellCyclePhase(), M_PHASE);

        // Set up SimulationTime
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0*TissueConfig::Instance()->GetCriticalHypoxicDuration(), num_steps);

        // Create a cell with a simple oxygen-based cell cycle model
        SimpleOxygenBasedCellCycleModel* p_cell_model = new SimpleOxygenBasedCellCycleModel;
        p_cell_model->SetDimension(2);
        TissueCell apoptotic_cell(STEM, p_state, p_cell_model);

        // Set up constant oxygen_concentration
        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);

        // Force the cell to be apoptotic
        for (unsigned i=0; i<num_steps; i++)
        {
            TS_ASSERT(apoptotic_cell.GetCellProliferativeType()!=APOPTOTIC ||
                      p_simulation_time->GetTime() >= TissueConfig::Instance()->GetCriticalHypoxicDuration());
            p_simulation_time->IncrementTimeOneStep();

            // Note that we need to pass in the updated G1 duration
            apoptotic_cell.ReadyToDivide();
        }

        // Test that the cell type is updated to be APOPTOTIC
        TS_ASSERT(apoptotic_cell.GetCellProliferativeType()==APOPTOTIC);
        TS_ASSERT_EQUALS(p_cell_model->GetCurrentHypoxicDuration(), 2.04);

        // Tidy up
        CellwiseData<2>::Destroy();

        // For coverage, create a 1D model
        CellwiseData<1>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        SimpleOxygenBasedCellCycleModel* p_cell_model1d = new SimpleOxygenBasedCellCycleModel;
        p_cell_model1d->SetDimension(1);
        TissueCell cell1d(STEM, p_state, p_cell_model1d);

        cell1d.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model1d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<1>::Destroy();

        // For coverage, create a 3D model
        CellwiseData<3>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        SimpleOxygenBasedCellCycleModel* p_cell_model3d = new SimpleOxygenBasedCellCycleModel;
        p_cell_model3d->SetDimension(3);
        TissueCell cell3d(STEM, p_state, p_cell_model3d);

        cell3d.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model3d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<3>::Destroy();
    }


    void TestStochasticOxygenBasedCellCycleModel() throw(Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();
        p_params->SetHepaOneParameters();

        // Check that mCurrentHypoxiaOnsetTime and mCurrentHypoxicDuration
        // are updated correctly
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);


        StochasticOxygenBasedCellCycleModel* p_model = new StochasticOxygenBasedCellCycleModel;
        p_model->SetDimension(2);
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        TissueCell cell(STEM, p_state, p_model);

        cell.InitialiseCellCycleModel();

        // Set up constant oxygen_concentration
        std::vector<double> low_oxygen_concentration;
        std::vector<double> high_oxygen_concentration;
        low_oxygen_concentration.push_back(0.0);
        high_oxygen_concentration.push_back(1.0);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);

        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);

        p_simulation_time->IncrementTimeOneStep(); // t=1.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(high_oxygen_concentration);

        p_simulation_time->IncrementTimeOneStep(); // t=2.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 2.0, 1e-12);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);
        p_simulation_time->IncrementTimeOneStep(); // t=3.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 2.0, 1e-12);

        // Set up simulation time
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            4.0*(p_params->GetHepaOneCellG1Duration()
                  +p_params->GetSG2MDuration()     ), num_steps);

        // Set up constant oxygen_concentration
        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        TS_ASSERT_THROWS_NOTHING(StochasticOxygenBasedCellCycleModel model);

        // Create cell cycle model
        StochasticOxygenBasedCellCycleModel* p_hepa_one_model = new StochasticOxygenBasedCellCycleModel;
        p_hepa_one_model->SetDimension(2);
        StochasticOxygenBasedCellCycleModel* p_diff_model = new StochasticOxygenBasedCellCycleModel;
        p_diff_model->SetDimension(2);

        // Create cell
        TissueCell hepa_one_cell(STEM, p_state, p_hepa_one_model);
        hepa_one_cell.InitialiseCellCycleModel();

        TissueCell diff_cell(DIFFERENTIATED, p_state, p_diff_model);
        diff_cell.InitialiseCellCycleModel();

        // Check that the cell cycle phase and ready to divide
        // are updated correctly
        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model->GetCurrentCellCyclePhase(), M_PHASE);

        TS_ASSERT_EQUALS(p_diff_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_diff_model->GetCurrentCellCyclePhase(), G_ZERO_PHASE);

        for (unsigned i=0; i<num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // Note that we need to pass in the updated G1 duration
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, p_hepa_one_model->GetG1Duration(), p_hepa_one_model->GetG2Duration());
        }

        TS_ASSERT_DELTA(p_hepa_one_model->GetAge(), p_simulation_time->GetTime(), 1e-9);
        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(), true);

        // Coverage
        TS_ASSERT_EQUALS(hepa_one_cell.ReadyToDivide(), true);
        TissueCell hepa_one_cell_divide = hepa_one_cell.Divide();

        // Check that cell division correctly resets the cell cycle phase
        StochasticOxygenBasedCellCycleModel* p_hepa_one_model2 = static_cast <StochasticOxygenBasedCellCycleModel*> (p_hepa_one_model->CreateCellCycleModel());

        TissueCell hepa_one_cell2(STEM, p_state, p_hepa_one_model2);
        TS_ASSERT_EQUALS(p_hepa_one_model2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model2->GetCurrentCellCyclePhase(), M_PHASE);

        // Set up SimulationTime
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0*TissueConfig::Instance()->GetCriticalHypoxicDuration(), num_steps);

        // Create a cell with a simple oxygen-based cell cycle model
        StochasticOxygenBasedCellCycleModel* p_cell_model = new StochasticOxygenBasedCellCycleModel;
        p_cell_model->SetDimension(2);
        TissueCell apoptotic_cell(STEM, p_state, p_cell_model);

        apoptotic_cell.InitialiseCellCycleModel();

        // Set up constant oxygen_concentration
        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);

        // Force the cell to be apoptotic
        for (unsigned i=0; i<num_steps; i++)
        {
            TS_ASSERT(apoptotic_cell.GetCellProliferativeType()!=APOPTOTIC ||
                      p_simulation_time->GetTime() >= TissueConfig::Instance()->GetCriticalHypoxicDuration());
            p_simulation_time->IncrementTimeOneStep();

            // Note that we need to pass in the updated G1 duration
            apoptotic_cell.ReadyToDivide();
        }

        // Test that the cell type is updated to be APOPTOTIC
        TS_ASSERT_EQUALS(apoptotic_cell.GetCellProliferativeType(), APOPTOTIC);
        TS_ASSERT_EQUALS(p_cell_model->GetCurrentHypoxicDuration(), 2.04);


        // Coverage
        TissueConfig::Instance()->SetMinimumGapDuration(1e20);

        StochasticOxygenBasedCellCycleModel* p_cell_model2 = new StochasticOxygenBasedCellCycleModel;
        p_cell_model2->SetDimension(2);
        TissueCell cell2(STEM, p_state, p_cell_model2);
        cell2.InitialiseCellCycleModel();

        TS_ASSERT_DELTA(p_cell_model2->GetG2Duration(), 1e20, 1e-4);

        // Tidy up
        CellwiseData<2>::Destroy();

        // For coverage, create a 1D model
        CellwiseData<1>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        StochasticOxygenBasedCellCycleModel* p_cell_model1d = new StochasticOxygenBasedCellCycleModel;
        p_cell_model1d->SetDimension(1);
        TissueCell cell1d(STEM, p_state, p_cell_model1d);

        cell1d.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model1d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<1>::Destroy();

        // For coverage, create a 3D model
        CellwiseData<3>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        StochasticOxygenBasedCellCycleModel* p_cell_model3d = new StochasticOxygenBasedCellCycleModel;
        p_cell_model3d->SetDimension(3);
        TissueCell cell3d(STEM, p_state, p_cell_model3d);

        cell3d.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model3d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<3>::Destroy();
    }


    void TestArchiveFixedDurationGenerationBasedCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "fixed_cell_cycle.arch";

        // Create an output archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 4);
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel;
            boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

            TissueCell cell(TRANSIT, p_healthy_state, p_model);
            cell.InitialiseCellCycleModel();

            p_simulation_time->IncrementTimeOneStep();
            p_simulation_time->IncrementTimeOneStep();

            p_model->SetBirthTime(-1.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Update cell phase
            p_model->ReadyToDivide();

            TissueCell* const p_cell = &cell;

            output_arch << p_cell;
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), G_ONE_PHASE);

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

            // Check private data has been restored correctly.
            TS_ASSERT_DELTA(p_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(), 2.5, 1e-12);
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), G_ONE_PHASE);

            delete p_cell;
        }
    }


    void TestArchiveStochasticDurationGenerationBasedCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "stoch_cell_cycle.arch";

        double random_number_test = 0;

        // Create an output archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            StochasticDurationGenerationBasedCellCycleModel* p_model = new StochasticDurationGenerationBasedCellCycleModel;

            boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);
            TissueCell cell(TRANSIT,  p_healthy_state, p_model);
            cell.InitialiseCellCycleModel();
            cell.SetBirthTime(-1.1);
            p_simulation_time->IncrementTimeOneStep();
            p_simulation_time->IncrementTimeOneStep();

            cell.ReadyToDivide(); // updates phases

            TissueCell* const p_cell = &cell;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            TS_ASSERT_DELTA(TissueConfig::Instance()->GetSDuration(), 5.0, 1e-12);

            output_arch << p_cell;

            TS_ASSERT_DELTA(p_model->GetBirthTime(), -1.1, 1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(), 2.1, 1e-12);
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), G_ONE_PHASE);

            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            random_number_test = p_gen->ranf();
            SimulationTime::Destroy();
        }

        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            p_gen->Reseed(128);

            TissueCell* p_cell;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            TissueConfig* p_inst1 = TissueConfig::Instance();

            p_inst1->SetSDuration(101.0);

            // Restore from the archive
            input_arch >> p_cell;

            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(), random_number_test, 1e-7);

            AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();

            // Check
            TS_ASSERT_DELTA(p_model->GetBirthTime(), -1.1, 1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(), 2.1, 1e-12);
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), G_ONE_PHASE);

            TS_ASSERT_DELTA(p_inst1->GetSDuration(), 5.0, 1e-12);

            // Tidy up
            delete p_cell;
        }
    }


    void TestArchiveSimpleWntCellCycleModel()
    {
        TissueConfig* p_params = TissueConfig::Instance();

        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "simple_wnt_cell_cycle.arch";

        // Set up the Wnt concentration
        double wnt_level = 1.0;
        WntConcentration<1>::Instance()->SetConstantWntValueForTesting(wnt_level);

        double random_number_test = 0;

        // Create an output archive
        {
            // Set up the simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();

            // The number for the G1 duration is taken from
            // the first random number generated
            double g1_duration = 1.0676;

            double end_time = g1_duration + p_params->GetSG2MDuration() + 5.0;

            unsigned num_timesteps = 50;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

            // Set up the Wnt concentration for testing
            WntConcentration<1>::Instance()->SetConstantWntValueForTesting(0.7);

            // Create cell cycle model and associated cell
            SimpleWntCellCycleModel* p_cell_model = new SimpleWntCellCycleModel;
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

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            TissueCell* const p_cell = &stem_cell;

            output_arch << p_cell;

            TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT(stem_cell.GetCellCycleModel()->ReadyToDivide());

            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            random_number_test = p_gen->ranf();
            SimulationTime::Destroy();
        }
        {
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
            TS_ASSERT_EQUALS((static_cast<SimpleWntCellCycleModel*>(p_cell_model))->GetDimension(), 1u);

            // Tidy up
            SimulationTime::Destroy();
            delete p_cell;
        }

        /*
         * Test the case of a radial Wnt concentration
         */

        p_params->Reset();
        RandomNumberGenerator::Instance()->Reseed(0);

        OutputFileHandler handler2("archive", false);
        archive_filename = handler2.GetOutputDirectoryFullPath() + "crypt_projection_cell_cycle.arch";

        // Set up the Wnt concentration
        wnt_level = p_params->GetWntStemThreshold() - 0.01;
        WntConcentration<2>::Destroy();
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        random_number_test = 0;

        // Create an output archive
        {
            // Set up the simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();

            // The number for the G1 duration is taken from
            // the first random number generated
            double g1_duration = 1.0676;

            double end_time = g1_duration + p_params->GetSG2MDuration() + 5.0;

            unsigned num_timesteps = 50;
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

            SimpleWntCellCycleModel* p_cell_model = new SimpleWntCellCycleModel;
            p_cell_model->SetDimension(2);
            p_cell_model->SetBirthTime(-1.0);

            boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

            TissueCell cell(STEM, p_healthy_state, p_cell_model);
            cell.InitialiseCellCycleModel();

            // Run to division age minus one time step to match birth time
            while (p_cell_model->GetAge() < g1_duration + p_params->GetSG2MDuration()
                                            - p_simulation_time->GetTimeStep())
            {
                p_simulation_time->IncrementTimeOneStep();
                CheckReadyToDivideAndPhaseIsUpdated(p_cell_model, g1_duration);
            }

            // Wnt should change this to a transit cell
            TS_ASSERT_EQUALS(cell.GetCellProliferativeType(), TRANSIT);
            TS_ASSERT_EQUALS(cell.GetCellCycleModel()->ReadyToDivide(), false);
            TS_ASSERT_EQUALS(cell.GetCellCycleModel()->GetCurrentCellCyclePhase(), G_TWO_PHASE);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            TissueCell* const p_cell = &cell;

            output_arch << p_cell;

            TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT(cell.GetCellCycleModel()->ReadyToDivide());

            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            random_number_test = p_gen->ranf();

            RandomNumberGenerator::Destroy();
            SimulationTime::Destroy();
        }
        {
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
            TS_ASSERT_EQUALS(p_cell->GetCellProliferativeType(), TRANSIT);

            p_simulation_time->IncrementTimeOneStep();

            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), true);
            TS_ASSERT_EQUALS(p_cell->GetCellProliferativeType(), TRANSIT);

            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_inst1->GetSG2MDuration(), 10.0, 1e-12);

            TS_ASSERT_DELTA(p_gen->ranf(), random_number_test, 1e-7);

            delete p_cell;
        }

        // Tidy up
        WntConcentration<1>::Destroy();
        WntConcentration<2>::Destroy();
    }


    void TestArchiveSimpleOxygenBasedCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "oxygen_based_cell_cycle.arch";

        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<1>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        // Create an output archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
            p_model->SetDimension(1);

            p_simulation_time->IncrementTimeOneStep();

            p_model->SetBirthTime(-1.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << static_cast<const SimpleOxygenBasedCellCycleModel&>(*p_model);

            // Tidy up
            delete p_model;
            SimulationTime::Destroy();
        }

        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
            p_model->SetDimension(2);
            p_model->SetBirthTime(-2.0);

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> *p_model;

            // Check that archiving worked correctly
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), M_PHASE);
            TS_ASSERT_EQUALS(p_model->GetDimension(), 1u);

            TS_ASSERT_DELTA(p_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(), 1.5, 1e-12);

            // Tidy up
            delete p_model;
            CellwiseData<1>::Destroy();
        }
    }


    void TestArchiveStochasticOxygenBasedCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "stochastic_oxygen_based_cell_cycle.arch";

        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<3>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        // Create an output archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // Create cell cycle model and associated cell
            StochasticOxygenBasedCellCycleModel* p_cell_model = new StochasticOxygenBasedCellCycleModel;
            p_cell_model->SetDimension(3);
            boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
            TissueCell cell(STEM, p_state, p_cell_model);

            cell.InitialiseCellCycleModel();
            cell.GetCellCycleModel()->SetBirthTime(-1.0);

            p_simulation_time->IncrementTimeOneStep();

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive cell
            TissueCell* const p_cell = &cell;
            output_arch << p_cell;

            // Tidy up
            SimulationTime::Destroy();
        }

        {
            // Set up simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            TissueConfig* inst1 = TissueConfig::Instance();

            inst1->SetSDuration(101.0);

            TissueCell* p_cell;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_cell;

            // Check that archiving worked correctly
            StochasticOxygenBasedCellCycleModel* p_model = static_cast<StochasticOxygenBasedCellCycleModel*> (p_cell->GetCellCycleModel());

            TS_ASSERT_EQUALS(p_cell, p_model->GetCell());
            TS_ASSERT_EQUALS(p_model->GetDimension(), 3u);
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), M_PHASE);

            TS_ASSERT_DELTA(p_model->GetBirthTime(), -1.0, 1e-4);
            TS_ASSERT_DELTA(p_model->GetAge(), 1.5, 1e-4);
            TS_ASSERT_DELTA(p_model->GetG2Duration(), 3.0676, 1e-4); // first random number generated

            // Tidy up
            delete p_cell;
        }
    }
};

#endif /*TESTSIMPLECELLCYCLEMODELS_HPP_*/
