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

#include "SingleOdeWntCellCycleModelCellsGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"

class TestCellBasedSingleOdeWnt : public AbstractCellBasedTestSuite
{
public:

    void TestSingleOdeWntCellCycleModel() throw(Exception)
    {
        // Set up the simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        double end_time = 1200.0;
        unsigned num_timesteps = 10*(unsigned)end_time;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);
        std::vector<double> initial_protein_concs;
        initial_protein_concs.push_back(0.0);
        initial_protein_concs.push_back(0.0);
        initial_protein_concs.push_back(1.0);

        // Set up the Wnt concentration
        double wnt_level = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);
        unsigned dimension = 2;
        TissueCell cell(STEM, HEALTHY, new SingleOdeWntCellCycleModel(dimension)); // Deal with memory properly
        SingleOdeWntCellCycleModel* p_cycle_model = static_cast<SingleOdeWntCellCycleModel*> (cell.GetCellCycleModel());

        cell.InitialiseCellCycleModel(); // Associates a cell with a cell cycle model.

        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            // Stem cell should have been changed into a transit cell by wnt cell cycle model
            TS_ASSERT_EQUALS(cell.GetCellProliferativeType(), TRANSIT);
            // The number for the G1 duration is taken from
            // the first random number generated
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, 1.0676);
        }

        double steady_beta_cat_at_wnt_equals_1 = p_cycle_model->GetBetaCateninConcentration();
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_1, 143.5119, 1e-4);

        // Divide the cell
        TS_ASSERT_EQUALS(cell.ReadyToDivide(), true);
        TissueCell cell2 = cell.Divide();
        cell.SetMutationState(LABELLED);

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

        // Test that both cells have inherited the same Beta-catenin concentration.
        double steady_beta_cat_at_wnt_equals_0_2 = 113.4683;
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_0_2, p_cycle_model->GetBetaCateninConcentration(), 1e-3);
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_0_2, p_cycle_model2->GetBetaCateninConcentration(), 1e-3);

        p_cycle_model->ResetForDivision();
        p_cycle_model2->ResetForDivision();

        // Test that both cells have inherited the same Beta-catenin concentration.
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_0_2, p_cycle_model->GetBetaCateninConcentration(), 1e-3);
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_0_2, p_cycle_model2->GetBetaCateninConcentration(), 1e-3);

        division_time = SimulationTime::Instance()->GetTime();

        // Now reduce the Wnt concentration so only mutant cells divide.
        wnt_level = 0.1;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Introduce a mutation (no immediate effect)
        cell2.SetMutationState(APC_ONE_HIT);
        TS_ASSERT(!p_cycle_model->ReadyToDivide());
        TS_ASSERT(!p_cycle_model2->ReadyToDivide());

        // Coverage...
        cell2.SetMutationState(APC_TWO_HIT);
        TS_ASSERT(!p_cycle_model->ReadyToDivide());
        cell2.SetMutationState(BETA_CATENIN_ONE_HIT);
        TS_ASSERT(!p_cycle_model->ReadyToDivide());

        // The numbers for the G1 durations are taken from
        // the next two random numbers generated
        new_g1_duration = 1.22037;
        new_g1_duration2 = 0.74699;

        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, new_g1_duration);
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model2, new_g1_duration2);
        }

        TS_ASSERT_DELTA(91.6693, p_cycle_model->GetBetaCateninConcentration(), 1e-3);
        TS_ASSERT_DELTA(355.7790, p_cycle_model2->GetBetaCateninConcentration(), 1e-3);

        TS_ASSERT_DELTA( p_cycle_model->GetBetaCateninDivisionThreshold(), 100, 1e-9);
        TS_ASSERT_DELTA(p_cycle_model2->GetBetaCateninDivisionThreshold(), 100, 1e-9);

        TS_ASSERT_EQUALS(cell.GetCellProliferativeType(), DIFFERENTIATED);
        TS_ASSERT_EQUALS(cell2.GetCellProliferativeType(), TRANSIT);

        // Coverage of 1D

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(30.0, 2);

        wnt_level = 1.0;
        WntConcentration<1>::Instance()->SetConstantWntValueForTesting(wnt_level);
        dimension = 1u;
        SingleOdeWntCellCycleModel* p_cell_model_1d = new SingleOdeWntCellCycleModel(dimension, true);

        TS_ASSERT_EQUALS(p_cell_model_1d->GetDimension(), 1u);

        TissueCell stem_cell_1d(STEM, HEALTHY, p_cell_model_1d);
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
        dimension = 3;
        WntConcentration<3>::Instance()->SetConstantWntValueForTesting(wnt_level);
        SingleOdeWntCellCycleModel* p_cell_model_3d = new SingleOdeWntCellCycleModel(dimension);

        TS_ASSERT_EQUALS(p_cell_model_3d->GetDimension(), 3u);

        TissueCell stem_cell_3d(STEM, HEALTHY, p_cell_model_3d);
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
};


#endif /* TESTSINGLEODEWNTCELLCYCLEMODEL_HPP_ */
