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
#ifndef TESTTISSUECELL_HPP_
#define TESTTISSUECELL_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>
#include <iostream>

#include "AbstractCellBasedTestSuite.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "TissueCell.hpp"
#include "CellMutationStateRegistry.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"


class TestTissueCell: public AbstractCellBasedTestSuite
{
public:

    void TestTissueCellConstructor() throw(Exception)
    {
    	// Set up SimulationTime
    	SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(200, 20);

        // Create a cell mutation state
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

        // Create a cell cycle model
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(TRANSIT);

        // Create a cell
        TissueCellPtr p_cell(new TissueCell(p_state, p_model));
        p_cell->SetBirthTime(-0.5);

        // Test members of cell directly
        TS_ASSERT_EQUALS(p_cell->GetCellId(), 0u);
        TS_ASSERT_DELTA(p_cell->GetBirthTime(), -0.5, 1e-6);

        // Test members of cell through cell cycle model
        TS_ASSERT_EQUALS(p_cell->GetCellCycleModel()->GetCell()->GetCellId(), 0u);
        TS_ASSERT_DELTA(p_cell->GetCellCycleModel()->GetCell()->GetBirthTime(), -0.5, 1e-6);
    }


    void TestCellsAgeingCorrectly() throw(Exception)
    {
        // These lines are added to cover the exception case that a cell is
        // created without simulation time being set up...
        FixedDurationGenerationBasedCellCycleModel fixed_model;
        SimulationTime::Destroy();

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        TS_ASSERT_THROWS_THIS(TissueCellPtr p_bad_cell(new TissueCell(p_healthy_state, &fixed_model)),
                              "TissueCell is setting up a cell cycle model but SimulationTime has not been set up");

        // Cell wasn't created - count should be zero
        TS_ASSERT_EQUALS(p_healthy_state->GetCellCount(), 0u);

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

        TS_ASSERT_THROWS_THIS(TissueCellPtr p_stem_cell(new TissueCell(p_healthy_state, NULL)),
                              "Cell cycle model is null");

        // Cell wasn't created - count should be zero
        TS_ASSERT_EQUALS(p_healthy_state->GetCellCount(), 0u);

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_stem_cell(new TissueCell(p_healthy_state, p_model));
        p_stem_cell->InitialiseCellCycleModel();

        // Cell was created - count should be one
        TS_ASSERT_EQUALS(p_healthy_state->GetCellCount(), 1u);

        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_stem_cell->GetAge(), 0.5);

        p_simulation_time->IncrementTimeOneStep();
        p_stem_cell->SetBirthTime(p_simulation_time->GetTime());
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_stem_cell->GetAge(), 1.0);

        TS_ASSERT_EQUALS(p_stem_cell->IsDead(), false);
        p_stem_cell->Kill();
        TS_ASSERT_EQUALS(p_stem_cell->IsDead(), true);

        // Coverage of operator equals.
        FixedDurationGenerationBasedCellCycleModel* p_model2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model2->SetCellProliferativeType(STEM);
        TissueCellPtr p_live_cell(new TissueCell(p_healthy_state, p_model2));

        TS_ASSERT_EQUALS(p_live_cell->IsDead(), false);
        p_live_cell = p_stem_cell;
        TS_ASSERT_EQUALS(p_live_cell->IsDead(), true);
    }


    void TestCellDivision()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(54.0, 9);

        // We are going to start at t=0 and jump up in steps of 6.0
        TissueConfig* p_params = TissueConfig::Instance();

        // This test needs particular cell cycle times
        TS_ASSERT_DELTA(p_params->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetSG2MDuration(), 10.0, 1e-12);

        p_simulation_time->IncrementTimeOneStep();//t=6

        // Cover bad cell cycle model

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        TS_ASSERT_THROWS_THIS(TissueCellPtr p_bad_cell2(new TissueCell(p_healthy_state, NULL)), "Cell cycle model is null");

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_stem_cell(new TissueCell(p_healthy_state, p_model));
        p_stem_cell->InitialiseCellCycleModel();

        // Test coverage of operator=
        FixedDurationGenerationBasedCellCycleModel* p_model2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model2->SetCellProliferativeType(TRANSIT);
        TissueCellPtr p_other_cell(new TissueCell(p_healthy_state, p_model2));
        p_other_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_other_cell->GetCellCycleModel()->GetCellProliferativeType(), TRANSIT);
        p_other_cell = p_stem_cell;
        TS_ASSERT_EQUALS(p_other_cell->GetCellCycleModel()->GetCellProliferativeType(), STEM);

        // Back to the test
        p_simulation_time->IncrementTimeOneStep(); //t=12
        p_simulation_time->IncrementTimeOneStep(); //t=18
        p_simulation_time->IncrementTimeOneStep(); //t=24

        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);

        p_simulation_time->IncrementTimeOneStep();//t=30

        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);

        // Create transit progeny of stem
        TissueCellPtr p_daughter_cell = p_stem_cell->Divide();

        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(static_cast<FixedDurationGenerationBasedCellCycleModel*>(p_daughter_cell->GetCellCycleModel())->GetGeneration(), 1u);
        TS_ASSERT_EQUALS(p_daughter_cell->GetCellCycleModel()->GetCellProliferativeType(), TRANSIT);
        TS_ASSERT_DELTA(p_daughter_cell->GetAge(), 0, 1e-9);

        p_simulation_time->IncrementTimeOneStep(); //t=36

        TS_ASSERT_EQUALS(p_daughter_cell->ReadyToDivide(), false);

        p_simulation_time->IncrementTimeOneStep(); //t=42

        TS_ASSERT_EQUALS(p_daughter_cell->ReadyToDivide(), true);

        // Create transit progeny of transit
        TissueCellPtr p_grandaughter_cell = p_daughter_cell->Divide();

        p_simulation_time->IncrementTimeOneStep(); //t=48

        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_grandaughter_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_daughter_cell->ReadyToDivide(), false);

        // Stem cell ready to divide again
        p_simulation_time->IncrementTimeOneStep(); //t=54
        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);

        // Both grandaughter and daughter cells should be ready to divide
        TS_ASSERT_EQUALS(p_grandaughter_cell->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_daughter_cell->ReadyToDivide(), true);
    }


    /*
     * ReadyToDivide() now calls UpdateCellProliferativeType() where appropriate.
     * (at the moment in Wnt-dependent cells).
     */
    void TestUpdateCellProliferativeTypes() throw (Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(200, 20);
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_stem_cell(new TissueCell(p_healthy_state, p_model));
        p_stem_cell->InitialiseCellCycleModel();
        p_stem_cell->ReadyToDivide();

        TS_ASSERT_EQUALS(p_stem_cell->GetCellCycleModel()->GetCellProliferativeType(),STEM);

        p_stem_cell->GetCellCycleModel()->SetCellProliferativeType(TRANSIT);

        p_stem_cell->ReadyToDivide();

        TS_ASSERT_EQUALS(p_stem_cell->GetCellCycleModel()->GetCellProliferativeType(),TRANSIT);

        // Test a Wnt dependent cell
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(0.0);

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        p_cell_cycle_model1->SetCellProliferativeType(TRANSIT);
        TissueCellPtr p_wnt_cell(new TissueCell(p_healthy_state, p_cell_cycle_model1));

        TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->GetCellProliferativeType(),TRANSIT);

        p_wnt_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->GetCellProliferativeType(),DIFFERENTIATED);

        p_wnt_cell->ReadyToDivide();

        TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->GetCellProliferativeType(),DIFFERENTIATED);

        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(1.0);

        // Go forward through time
        for (unsigned i=0; i<20; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }

        p_wnt_cell->ReadyToDivide();

        TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->GetCellProliferativeType(),TRANSIT);

        // Tidy up
        WntConcentration<2>::Destroy();
    }


    void Test0DBucket()
    {
        double end_time = 61.0;
        int time_steps = 61;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, time_steps);
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_stem_cell(new TissueCell(p_healthy_state, p_model));
        p_stem_cell->InitialiseCellCycleModel();

        std::vector<TissueCellPtr> cells;
        std::vector<TissueCellPtr> newly_born;
        std::vector<unsigned> stem_cells(time_steps);
        std::vector<unsigned> transit_cells(time_steps);
        std::vector<unsigned> differentiated_cells(time_steps);
        std::vector<double> times(time_steps);

        cells.push_back(p_stem_cell);
        std::vector<TissueCellPtr>::iterator cell_iterator;

        unsigned i=0;
        while (p_simulation_time->GetTime()< end_time)
        {
            // Produce the offspring of the cells
            p_simulation_time->IncrementTimeOneStep();
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                if ((*cell_iterator)->ReadyToDivide())
                {
                    newly_born.push_back((*cell_iterator)->Divide());
                }
                ++cell_iterator;
            }

            // Copy offspring in newly_born vector to cells vector
            cell_iterator = newly_born.begin();
            while (cell_iterator < newly_born.end())
            {
                cells.push_back(*cell_iterator);
                ++cell_iterator;
            }
            newly_born.clear();

            // Update cell counts
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                switch ((*cell_iterator)->GetCellCycleModel()->GetCellProliferativeType())
                {
                    case STEM:
                        stem_cells[i]++;
                        break;
                    case TRANSIT:
                        transit_cells[i]++;
                        break;
                    default:
                        differentiated_cells[i]++;
                        break;
                }

                ++cell_iterator;
            }
            times[i] = p_simulation_time->GetTime();
            i++;
        }
        TS_ASSERT_EQUALS(stem_cells[59], 1u);
        TS_ASSERT_EQUALS(transit_cells[59], 2u);
        TS_ASSERT_EQUALS(differentiated_cells[59], 8u);
    }


    void TestWithFixedDurationGenerationBasedCellCycleModel() throw(Exception)
    {
        // Simulation time is 6000 because we want to test that differentiated cells never divide.

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(6000.0, 1000);
        TissueConfig* p_params = TissueConfig::Instance();

        // This test needs particular cell cycle times
        TS_ASSERT_DELTA(p_params->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetSG2MDuration(), 10.0, 1e-12);

        p_simulation_time->IncrementTimeOneStep();
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        //  Creating different types of cells with different cell cycle models at SimulationTime = 6 hours
        FixedDurationGenerationBasedCellCycleModel* p_stem_model = new FixedDurationGenerationBasedCellCycleModel();
        p_stem_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_stem_cell(new TissueCell(p_healthy_state, p_stem_model));
        p_stem_cell->InitialiseCellCycleModel();

        StochasticDurationGenerationBasedCellCycleModel* p_stoch_model = new StochasticDurationGenerationBasedCellCycleModel();
        p_stoch_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_stochastic_stem_cell(new TissueCell(p_healthy_state, p_stoch_model));
        p_stochastic_stem_cell->InitialiseCellCycleModel();

        FixedDurationGenerationBasedCellCycleModel* p_diff_model = new FixedDurationGenerationBasedCellCycleModel();
        p_diff_model->SetCellProliferativeType(DIFFERENTIATED);
        p_diff_model->SetGeneration(6);
        TissueCellPtr p_differentiated_cell(new TissueCell(p_healthy_state, p_diff_model));
        p_differentiated_cell->InitialiseCellCycleModel();

        StochasticDurationGenerationBasedCellCycleModel* p_stoch_diff_model = new StochasticDurationGenerationBasedCellCycleModel();
        p_stoch_diff_model->SetCellProliferativeType(DIFFERENTIATED);
        p_stoch_diff_model->SetGeneration(6);
        TissueCellPtr p_stochastic_differentiated_cell(new TissueCell(p_healthy_state, p_stoch_diff_model));
        p_stochastic_differentiated_cell->InitialiseCellCycleModel();

        FixedDurationGenerationBasedCellCycleModel* p_transit_model = new FixedDurationGenerationBasedCellCycleModel();
        p_transit_model->SetCellProliferativeType(TRANSIT);
        p_transit_model->SetGeneration(2);
        TissueCellPtr p_transit_cell(new TissueCell(p_healthy_state, p_transit_model));
        p_transit_cell->InitialiseCellCycleModel();

        // SimulationTime = 6 hours
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        // SimulationTime = 18 hours
        TS_ASSERT_EQUALS(p_stochastic_stem_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_transit_cell->ReadyToDivide(), true);

        p_simulation_time->IncrementTimeOneStep();

        // SimulationTime = 24 hours
        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_stochastic_stem_cell->ReadyToDivide(), true);

        p_simulation_time->IncrementTimeOneStep();

        // SimulationTime = 30 hours
        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_stochastic_stem_cell->ReadyToDivide(), true);

        TissueCellPtr p_daughter_cell1 = p_stem_cell->Divide();
        TS_ASSERT(typeid(p_daughter_cell1->GetCellCycleModel()) == typeid(p_stem_cell->GetCellCycleModel()));

        // Go to large time to ensure that differentiated cells can not divide
        for (unsigned i=0; i<990; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }
        TS_ASSERT_EQUALS(p_differentiated_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_stochastic_differentiated_cell->ReadyToDivide(), false);
    }


    void TestStochasticCycleModel() throw(Exception)
    {
        // Go up in steps of 0.01 to test stochasticity in cell cycle models
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(54.0, 5400);
        TissueConfig* p_params = TissueConfig::Instance();

        // This test needs particular cell cycle times
        TS_ASSERT_DELTA(p_params->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetSG2MDuration(), 10.0, 1e-12);

        for (unsigned i=0; i<600; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }

        // Now at t=6.00
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        FixedDurationGenerationBasedCellCycleModel* p_transit_model = new FixedDurationGenerationBasedCellCycleModel();
        p_transit_model->SetCellProliferativeType(TRANSIT);
        p_transit_model->SetGeneration(2);
        TissueCellPtr p_transit_cell(new TissueCell(p_healthy_state, p_transit_model));
        p_transit_cell->InitialiseCellCycleModel();

        for (unsigned i=0; i<1199; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }

        // Now at t = 17.99, cell is 11.99 old
        TS_ASSERT_EQUALS(p_transit_cell->ReadyToDivide(), false);

        StochasticDurationGenerationBasedCellCycleModel* p_cell_cycle_model = new StochasticDurationGenerationBasedCellCycleModel;
        p_cell_cycle_model->SetCellProliferativeType(TRANSIT);

        // This now resets the age of the cell to 0.0 so more time added in underneath
        p_transit_cell->SetCellCycleModel(p_cell_cycle_model);
        p_transit_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_transit_cell->GetCellCycleModel(), p_cell_cycle_model);
        for (unsigned i=0; i<1399; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }
        TS_ASSERT_EQUALS(p_transit_cell->ReadyToDivide(), true);

        // Ensure transit cell divides
        while (!p_transit_cell->ReadyToDivide())
        {
            p_simulation_time->IncrementTimeOneStep();
        }
        TissueCellPtr p_daughter_cell2 = p_transit_cell->Divide();
        TS_ASSERT(typeid(p_daughter_cell2->GetCellCycleModel()) == typeid(p_transit_cell->GetCellCycleModel()));
    }


    void Test0DBucketStochastic()
    {
        TissueConfig* p_params = TissueConfig::Instance();

        // This test needs particular cell cycle times
        TS_ASSERT_DELTA(p_params->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetSG2MDuration(), 10.0, 1e-12);

        const double end_time = 70.0;
        const unsigned number_of_simulations = 1000;

        std::vector<TissueCellPtr> cells;
        std::vector<TissueCellPtr> newly_born;

        std::vector<unsigned> stem_cells(number_of_simulations);
        std::vector<unsigned> transit_cells(number_of_simulations);
        std::vector<unsigned> differentiated_cells(number_of_simulations);
        double stem_cell_mean = 0.0;
        double transit_cell_mean = 0.0;
        double differentiated_cell_mean = 0.0;

        for (unsigned simulation_number=0; simulation_number<number_of_simulations; simulation_number++)
        {
            SimulationTime::Destroy();
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(70.0, 70);
            boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

            StochasticDurationGenerationBasedCellCycleModel* p_model = new StochasticDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);
            TissueCellPtr p_stem_cell(new TissueCell(p_healthy_state, p_model));
            p_stem_cell->InitialiseCellCycleModel();
            cells.push_back(p_stem_cell);

            // Produce the offspring of the cells
            std::vector<TissueCellPtr>::iterator cell_iterator = cells.begin();

            while (p_simulation_time->GetTime()< end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                cell_iterator = cells.begin();
                while (cell_iterator < cells.end())
                {
                    if ((*cell_iterator)->ReadyToDivide())
                    {
                        newly_born.push_back((*cell_iterator)->Divide());
                    }
                    ++cell_iterator;
                }

                // Copy offspring in newly_born vector to cells vector
                cell_iterator = newly_born.begin();
                while (cell_iterator < newly_born.end())
                {
                    cells.push_back(*cell_iterator);
                    ++cell_iterator;
                }

                newly_born.clear();
            }
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                switch ((*cell_iterator)->GetCellCycleModel()->GetCellProliferativeType())
                {
                    case STEM:
                        stem_cells[simulation_number]++;
                        stem_cell_mean++;
                        break;
                    case TRANSIT:
                        transit_cells[simulation_number]++;
                        transit_cell_mean++;
                        break;
                    default:
                        differentiated_cells[simulation_number]++;
                        differentiated_cell_mean++;
                        break;
                }

                ++cell_iterator;
            }
            cells.clear();
        }
        stem_cell_mean=stem_cell_mean/(double) number_of_simulations;
        transit_cell_mean=transit_cell_mean/(double) number_of_simulations;
        differentiated_cell_mean=differentiated_cell_mean/(double) number_of_simulations;

        TS_ASSERT_DELTA(stem_cell_mean, 1.0, 1e-12);
        TS_ASSERT_DELTA(transit_cell_mean, 6.84, 1.0);
        TS_ASSERT_DELTA(differentiated_cell_mean, 16.0, 1.0);

        TS_ASSERT_DELTA(p_params->GetTransitCellG1Duration(), 2.0, 1e-12);
    }


    void TestInitialise0DBucket()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(60.0, 60);

        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        FixedDurationGenerationBasedCellCycleModel* p_stem_model = new FixedDurationGenerationBasedCellCycleModel();
        p_stem_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_stem_cell(new TissueCell(p_healthy_state, p_stem_model));
        p_stem_cell->InitialiseCellCycleModel();
        cells.push_back(p_stem_cell);

        FixedDurationGenerationBasedCellCycleModel* p_transit_model = new FixedDurationGenerationBasedCellCycleModel();
        p_transit_model->SetCellProliferativeType(TRANSIT);
        p_transit_model->SetGeneration(1);
        TissueCellPtr p_transit_cell_1(new TissueCell(p_healthy_state, p_transit_model));
        p_transit_cell_1->InitialiseCellCycleModel();
        cells.push_back(p_transit_cell_1);

        FixedDurationGenerationBasedCellCycleModel* p_transit_model2 = new FixedDurationGenerationBasedCellCycleModel();
        p_transit_model2->SetCellProliferativeType(TRANSIT);
        p_transit_model2->SetGeneration(2);
        TissueCellPtr p_transit_cell_2(new TissueCell(p_healthy_state, p_transit_model2));
        p_transit_cell_2->InitialiseCellCycleModel();
        cells.push_back(p_transit_cell_2);

        FixedDurationGenerationBasedCellCycleModel* p_transit_model3 = new FixedDurationGenerationBasedCellCycleModel();
        p_transit_model3->SetCellProliferativeType(TRANSIT);
        p_transit_model3->SetGeneration(3);
        TissueCellPtr p_transit_cell_3(new TissueCell(p_healthy_state, p_transit_model3));
        p_transit_cell_3->InitialiseCellCycleModel();
        cells.push_back(p_transit_cell_3);

        FixedDurationGenerationBasedCellCycleModel* p_diff_model = new FixedDurationGenerationBasedCellCycleModel();
        p_diff_model->SetCellProliferativeType(DIFFERENTIATED);
        p_diff_model->SetGeneration(4);
        TissueCellPtr p_differentiated_cell(new TissueCell(p_healthy_state, p_diff_model));
        p_differentiated_cell->InitialiseCellCycleModel();
        cells.push_back(p_differentiated_cell);

        std::vector<TissueCellPtr> newly_born;
        std::vector<unsigned> stem_cells(p_simulation_time->GetTotalNumberOfTimeSteps());
        std::vector<unsigned> transit_cells(p_simulation_time->GetTotalNumberOfTimeSteps());
        std::vector<unsigned> differentiated_cells(p_simulation_time->GetTotalNumberOfTimeSteps());
        std::vector<double> times(p_simulation_time->GetTotalNumberOfTimeSteps());

        std::vector<TissueCellPtr>::iterator cell_iterator;

        unsigned i = 0;
        while (!p_simulation_time->IsFinished())
        {
            p_simulation_time->IncrementTimeOneStep();

            // Produce the offspring of the cells
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                TS_ASSERT_EQUALS((*cell_iterator)->GetMutationState()->IsType<WildTypeCellMutationState>(), true);
                if ((*cell_iterator)->ReadyToDivide())
                {
                    newly_born.push_back((*cell_iterator)->Divide());
                }
                ++cell_iterator;
            }

            // Copy offspring in newly_born vector to cells vector
            cell_iterator = newly_born.begin();
            while (cell_iterator < newly_born.end())
            {
                cells.push_back(*cell_iterator);
                ++cell_iterator;
            }
            newly_born.clear();

            // Count number of cells of each type
            cell_iterator = cells.begin();
            stem_cells[i] = 0;
            transit_cells[i] = 0;
            differentiated_cells[i] = 0;
            while (cell_iterator < cells.end())
            {
                switch ((*cell_iterator)->GetCellCycleModel()->GetCellProliferativeType())
                {
                    case STEM:
                        stem_cells[i]++;
                        break;
                    case TRANSIT:
                        transit_cells[i]++;
                        break;
                    default:
                        differentiated_cells[i]++;
                        break;
                }
                ++cell_iterator;
            }

            times[i] = p_simulation_time->GetTime();
            i++;
        }

        TS_ASSERT_EQUALS(stem_cells[59], 1u);
        TS_ASSERT_EQUALS(transit_cells[59], 2u);
        TS_ASSERT_EQUALS(differentiated_cells[59], 23u);
    }


    /*
     * We are checking that the TissueCellPtrs work with the Wnt cell cycle models here
     * That division of wnt cells and stuff works OK.
     *
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModel() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        TissueConfig* p_parameters = TissueConfig::Instance();

        double SG2MDuration = p_parameters->GetSG2MDuration();

        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);

        double wnt_stimulus = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        p_cell_cycle_model1->SetCellProliferativeType(TRANSIT);
        TissueCellPtr p_wnt_cell(new TissueCell(p_healthy_state, p_cell_cycle_model1));
        p_wnt_cell->InitialiseCellCycleModel();

#ifdef CHASTE_CVODE
        const double expected_g1_duration = 5.96441;
#else
        const double expected_g1_duration = 5.971;
#endif //CHASTE_CVODE

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            if (time >= expected_g1_duration+SG2MDuration)
            {
                TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), false);
            }
        }

        p_simulation_time->IncrementTimeOneStep();
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);

        TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), true);

        TissueCellPtr p_wnt_cell2 = p_wnt_cell->Divide();

        double time_of_birth = p_wnt_cell->GetBirthTime();
        double time_of_birth2 = p_wnt_cell2->GetBirthTime();

        TS_ASSERT_DELTA(time_of_birth, time_of_birth2, 1e-9);

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            bool result1 = p_wnt_cell->ReadyToDivide();
            bool result2 = p_wnt_cell2->ReadyToDivide();

            if (time >= expected_g1_duration + SG2MDuration + time_of_birth)
            {
                TS_ASSERT_EQUALS(result1, true);
                TS_ASSERT_EQUALS(result2, true);
            }
            else
            {
                TS_ASSERT_EQUALS(result1, false);
                TS_ASSERT_EQUALS(result2, false);
            }
        }

        WntConcentration<2>::Destroy();
    }

    /*
     * We are checking that the TissueCellPtrs work with the StochasticWnt cell cycle models here
     * That division of wnt cells and stuff works OK.
     *
     * It checks that the cell division thing works nicely too.
     */
    void TestWithStochasticWntCellCycleModel() throw(Exception)
    {
        TissueConfig* p_parameters = TissueConfig::Instance();

        // These are the first three normal random with mean 10, s.d. 1 and this seed (0)
        double SG2MDuration1 = p_parameters->GetSDuration() + 3.16084 + p_parameters->GetMDuration();
        double SG2MDuration2 = p_parameters->GetSDuration() + 5.0468  + p_parameters->GetMDuration();
        double SG2MDuration3 = p_parameters->GetSDuration() + 3.34408 + p_parameters->GetMDuration();
        double g1_duration = 5.971;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);

        double wnt_stimulus = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        StochasticWntCellCycleModel* p_cell_model = new StochasticWntCellCycleModel();
        p_cell_model->SetDimension(2);
        p_cell_model->SetCellProliferativeType(TRANSIT);
        TissueCellPtr p_wnt_cell(new TissueCell(p_healthy_state, p_cell_model));
        p_wnt_cell->InitialiseCellCycleModel();

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            if (time >= g1_duration+SG2MDuration1)
            {
                TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), false);
            }
        }

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), true);

        TissueCellPtr p_wnt_cell2 = p_wnt_cell->Divide();

        double time_of_birth = p_wnt_cell->GetBirthTime();
        double time_of_birth2 = p_wnt_cell2->GetBirthTime();

        TS_ASSERT_DELTA(time_of_birth, time_of_birth2, 1e-9);

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            bool parent_ready = p_wnt_cell->ReadyToDivide();
            bool daughter_ready = p_wnt_cell2->ReadyToDivide();

            if (time >= g1_duration+SG2MDuration2+time_of_birth)
            {
                TS_ASSERT_EQUALS(parent_ready, true);
            }
            else
            {
                TS_ASSERT_EQUALS(parent_ready, false);
            }
            if (time >= g1_duration+SG2MDuration3+time_of_birth2)
            {
                TS_ASSERT_EQUALS(daughter_ready, true);
            }
            else
            {
                TS_ASSERT_EQUALS(daughter_ready, false);
            }
        }

        WntConcentration<2>::Destroy();
    }

    /*
     * We are checking that the TissueCellPtrs work with the T&N cell cycle models here
     * That division of wnt cells and stuff works OK.
     *
     * It checks that the cell division thing works nicely too.
     */
    void TestWithTysonNovakCellCycleModel() throw(Exception)
    {
        double standard_tyson_duration = 1.242;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(200.0/60.0, num_steps+1);

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        TysonNovakCellCycleModel* p_cell_model = new TysonNovakCellCycleModel();
        p_cell_model->SetCellProliferativeType(TRANSIT);
        TissueCellPtr p_tn_cell(new TissueCell(p_healthy_state, p_cell_model));
        p_tn_cell->InitialiseCellCycleModel();

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();
            if (time>standard_tyson_duration)
            {
                TS_ASSERT_EQUALS(p_tn_cell->ReadyToDivide(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(p_tn_cell->ReadyToDivide(), false);
            }
        }

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_tn_cell->ReadyToDivide(), true);

        TissueCellPtr p_tn_cell2 = p_tn_cell->Divide();

        double time_of_birth = p_tn_cell->GetBirthTime();
        double time_of_birth2 = p_tn_cell2->GetBirthTime();

        TS_ASSERT_DELTA(time_of_birth, time_of_birth2, 1e-9);

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();
            bool result1 = p_tn_cell->ReadyToDivide();
            bool result2 = p_tn_cell2->ReadyToDivide();

            if (time >= standard_tyson_duration + time_of_birth)
            {
                TS_ASSERT_EQUALS(result1, true);
                TS_ASSERT_EQUALS(result2, true);
            }
            else
            {
                TS_ASSERT_EQUALS(result1, false);
                TS_ASSERT_EQUALS(result2, false);
            }
        }
    }

    void TestTysonNovakSteadyState()
    {
        // Keep dividing until we reach steady-state
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps=100000;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(20000.0/60.0, num_steps+1);

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        TysonNovakCellCycleModel* p_cell_model = new TysonNovakCellCycleModel();
        p_cell_model->SetCellProliferativeType(TRANSIT);
        TissueCellPtr p_tn_cell(new TissueCell(p_healthy_state, p_cell_model));
        p_tn_cell->InitialiseCellCycleModel();

        unsigned num_divisions = 0;

        while (!p_simulation_time->IsFinished())
        {
            while (!p_simulation_time->IsFinished() && !p_tn_cell->ReadyToDivide())
            {
                p_simulation_time->IncrementTimeOneStep();
            }
            if (p_tn_cell->ReadyToDivide())
            {
                TissueCellPtr p_tn_cell2 = p_tn_cell->Divide();
                ++num_divisions;
            }
        }
        TS_ASSERT_EQUALS(num_divisions, 268u);
    }

    void TestApoptosisAndDeath()
    {
        // We are going to start at t=0 and jump up in steps of 0.2
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(0.6, 3);

        // This test needs particular apoptosis time
        TissueConfig* p_params = TissueConfig::Instance();
        TS_ASSERT_EQUALS(p_params->GetApoptosisTime(), 0.25);

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
        p_cell_model->SetCellProliferativeType(TRANSIT);
        TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_cell_model));
        p_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell->HasApoptosisBegun(), false);
        TS_ASSERT_EQUALS(p_cell->IsDead(), false);
        TS_ASSERT_THROWS_THIS(p_cell->GetTimeUntilDeath(),"Shouldn\'t be checking time until apoptosis as it isn\'t set");

        p_simulation_time->IncrementTimeOneStep(); // t=0.2

        p_cell->StartApoptosis();
        TS_ASSERT_THROWS_THIS(p_cell->StartApoptosis(),"StartApoptosis() called when already undergoing apoptosis");

        TS_ASSERT_EQUALS(p_cell->HasApoptosisBegun(), true);
        TS_ASSERT_EQUALS(p_cell->IsDead(), false);
        TS_ASSERT_DELTA(p_cell->GetTimeUntilDeath(),0.25,1e-12);

        // Check that we can copy a cell that has started apoptosis
        TissueCellPtr p_cell2 = p_cell;

        p_simulation_time->IncrementTimeOneStep(); // t=0.4
        TS_ASSERT_EQUALS(p_cell->HasApoptosisBegun(), true);
        TS_ASSERT_EQUALS(p_cell->IsDead(), false);
        TS_ASSERT_DELTA(p_cell->GetTimeUntilDeath(),0.05,1e-12);

        TS_ASSERT_EQUALS(p_cell2->HasApoptosisBegun(), true);
        TS_ASSERT_EQUALS(p_cell2->IsDead(), false);
        TS_ASSERT_DELTA(p_cell2->GetTimeUntilDeath(),0.05,1e-12);

        p_simulation_time->IncrementTimeOneStep(); // t=0.6
        TS_ASSERT_EQUALS(p_cell->HasApoptosisBegun(), true);
        TS_ASSERT_EQUALS(p_cell->IsDead(), true);
    }


    void TestCantDivideIfUndergoingApoptosis()
    {
        // We are going to start at t=0 and jump up to t=25
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 1);

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
        p_cell_model->SetCellProliferativeType(TRANSIT);
        TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_cell_model));
        p_cell->InitialiseCellCycleModel();
        p_simulation_time->IncrementTimeOneStep(); // t=25

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);
        p_cell->StartApoptosis();
        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), false);
    }


    void Test0DBucketWithDeath()
    {
        double end_time = 92.0;
        int time_steps = 92;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, time_steps);

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
        p_cell_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_stem_cell(new TissueCell(p_healthy_state, p_cell_model));
        p_stem_cell->InitialiseCellCycleModel();

        std::vector<TissueCellPtr> cells;
        std::vector<TissueCellPtr> newly_born;
        std::vector<unsigned> stem_cells(time_steps);
        std::vector<unsigned> transit_cells(time_steps);
        std::vector<unsigned> differentiated_cells(time_steps);
        std::vector<unsigned> dead_cells(time_steps);
        std::vector<double> times(time_steps);

        cells.push_back(p_stem_cell);
        std::vector<TissueCellPtr>::iterator cell_iterator;

        unsigned i=0;
        while (p_simulation_time->GetTime()< end_time)
        {
            // Produce the offspring of the cells

            p_simulation_time->IncrementTimeOneStep();
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                if (!(*cell_iterator)->IsDead())
                {
                    if ((*cell_iterator)->ReadyToDivide())
                    {
                        newly_born.push_back((*cell_iterator)->Divide());
                    }

                    if (((*cell_iterator)->GetAge() > 30))
                    {
                        (*cell_iterator)->StartApoptosis();
                    }
                }
                ++cell_iterator;
            }

            // Copy offspring in newly_born vector to cells vector
            cell_iterator = newly_born.begin();
            while (cell_iterator < newly_born.end())
            {
                cells.push_back(*cell_iterator);
                ++cell_iterator;
            }
            newly_born.clear();

            // Update cell counts
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                if (!(*cell_iterator)->IsDead())
                {
                    switch ((*cell_iterator)->GetCellCycleModel()->GetCellProliferativeType())
                    {
                        case STEM:
                            stem_cells[i]++;
                            break;
                        case TRANSIT:
                            transit_cells[i]++;
                            break;
                        default:
                            differentiated_cells[i]++;
                            break;
                    }
                }
                else
                {
                    dead_cells[i]++;
                }

                ++cell_iterator;
            }
            times[i] = p_simulation_time->GetTime();
            i++;
        }

        TS_ASSERT_EQUALS(stem_cells[time_steps-1], 1u);
        TS_ASSERT_EQUALS(transit_cells[time_steps-1], 2u);
        TS_ASSERT_EQUALS(differentiated_cells[time_steps-1], 8u);
        TS_ASSERT_EQUALS(dead_cells[time_steps-1], 8u);
    }


    void TestArchiveCell() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "cell->arch";

        // Archive a cell
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

            FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
            p_cell_model->SetCellProliferativeType(STEM);
            TissueCellPtr p_stem_cell(new TissueCell(p_healthy_state, p_cell_model));
            p_stem_cell->InitialiseCellCycleModel();
            p_simulation_time->IncrementTimeOneStep();

            TS_ASSERT_EQUALS(p_stem_cell->GetAge(), 0.5);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            TissueCellPtr const p_cell = p_stem_cell;

            // Write the cell to the archive
            output_arch << static_cast<const SimulationTime&> (*p_simulation_time);
            output_arch << p_cell;

            SimulationTime::Destroy();
        }

        // Restore TissueCellPtr
        {
            // Need to set up time to initialize a cell
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(1.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 1); // will be restored

            // Initialize a cell

            TissueCellPtr p_stem_cell;

            // Restore the cell
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> *p_simulation_time;
            input_arch >> p_stem_cell;

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.5);
            TS_ASSERT_EQUALS(p_simulation_time->GetTimeStep(), 0.5);

            TS_ASSERT_EQUALS(p_stem_cell->GetAge(), 0.5);
            TS_ASSERT_EQUALS(static_cast<FixedDurationGenerationBasedCellCycleModel*>(p_stem_cell->GetCellCycleModel())->GetGeneration(), 0u);
            TS_ASSERT_EQUALS(p_stem_cell->GetCellCycleModel()->GetCellProliferativeType(), STEM);

            AbstractCellCycleModel* p_model = p_stem_cell->GetCellCycleModel();

            TS_ASSERT_EQUALS(p_model->GetCell(), p_stem_cell);
        }
    }

    /*
     * We are checking that the TissueCellPtrs work with the Wnt
     * cell cycle models here. This just tests the set-up and checks that
     * the functions can all be called (not what they return).
     *
     * For more in depth tests see TestNightlyTissueCellPtr.hpp
     * (these test that the cell cycle times are correct for the
     * various mutant cells)
     */
    void TestWntMutantVariantsAndLabelling() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 10;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

        double wnt_stimulus = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);

        boost::shared_ptr<AbstractCellMutationState> p_apc_one_hit_state(CellMutationStateRegistry::Instance()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellMutationState> p_apc_two_hit_state(CellMutationStateRegistry::Instance()->Get<ApcTwoHitCellMutationState>());
        boost::shared_ptr<AbstractCellMutationState> p_bcat_one_hit_state(CellMutationStateRegistry::Instance()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellMutationState> p_labelled_state(CellMutationStateRegistry::Instance()->Get<LabelledCellMutationState>());

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        p_cell_cycle_model1->SetCellProliferativeType(TRANSIT);
        TissueCellPtr p_wnt_cell(new TissueCell(p_apc_one_hit_state, p_cell_cycle_model1));
        p_wnt_cell->InitialiseCellCycleModel();

        WntCellCycleModel* p_cell_cycle_model2 = new WntCellCycleModel();
        p_cell_cycle_model2->SetDimension(2);
        p_cell_cycle_model2->SetCellProliferativeType(TRANSIT);
        TissueCellPtr p_wnt_cell2(new TissueCell(p_bcat_one_hit_state, p_cell_cycle_model2));
        p_wnt_cell2->InitialiseCellCycleModel();

        WntCellCycleModel* p_cell_cycle_model3 = new WntCellCycleModel();
        p_cell_cycle_model3->SetDimension(2);
        p_cell_cycle_model3->SetCellProliferativeType(TRANSIT);
        TissueCellPtr p_wnt_cell3(new TissueCell(p_apc_two_hit_state, p_cell_cycle_model3));
        p_wnt_cell3->InitialiseCellCycleModel();

        WntCellCycleModel* p_cell_cycle_model4 = new WntCellCycleModel();
        p_cell_cycle_model4->SetDimension(2);
        p_cell_cycle_model4->SetCellProliferativeType(TRANSIT);
        TissueCellPtr p_wnt_cell4(new TissueCell(p_labelled_state, p_cell_cycle_model4));
        p_wnt_cell4->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_wnt_cell2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_wnt_cell3->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_wnt_cell4->ReadyToDivide(), false);

        WntConcentration<2>::Destroy();
    }

    void TestIsLogged()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 1);

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
        p_cell_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_cell_model));
        p_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell->IsLogged(), false);

        p_cell->SetLogged();

        TS_ASSERT_EQUALS(p_cell->IsLogged(), true);

        TissueCellPtr p_copied_cell = p_cell;

        TS_ASSERT_EQUALS(p_copied_cell->IsLogged(), true);

        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

        TissueCellPtr p_daughter_cell = p_cell->Divide();

        TS_ASSERT_EQUALS(p_cell->IsLogged(), true);
        TS_ASSERT_EQUALS(p_daughter_cell->IsLogged(), false);
    }

    void TestAncestors() throw (Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
        p_cell_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_cell_model));
        p_cell->InitialiseCellCycleModel();
        p_cell->SetAncestor(2u);

        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

        TissueCellPtr p_cell2 = p_cell->Divide();

        TS_ASSERT_EQUALS(p_cell->GetAncestor(), 2u);
        TS_ASSERT_EQUALS(p_cell2->GetAncestor(), 2u);
    }

    void TestCellId() throw (Exception)
    {
        // Resetting the Maximum cell Id to zero (to account for previous tests)
        TissueCell::ResetMaxCellId();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
        p_cell_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_cell(new TissueCell(p_healthy_state, p_cell_model));
        p_cell->InitialiseCellCycleModel();

        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

        TissueCellPtr p_cell2 = p_cell->Divide();

        TS_ASSERT_EQUALS(p_cell->GetCellId(), 0u);
        TS_ASSERT_EQUALS(p_cell2->GetCellId(), 1u);
    }

    void ThisTestBreaksSomeBuildsTestCellDivisionStops()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(54.0, 9);
        TissueConfig* p_params = TissueConfig::Instance();

        // If the value of GetStemCellG1Duration() changes in p_params the simulation time
        // step and end time will need to be changed accordingly so that
        // IncrementTimeOneStep() gets the cell to correct division times

        // This test needs particular cell cycle times
        TS_ASSERT_DELTA(p_params->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_params->GetSG2MDuration(), 10.0, 1e-12);

        // SimulationTime returns 0 hours
        p_simulation_time->IncrementTimeOneStep();
        // SimulationTime returns 6 hours

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());

        FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
        p_cell_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_stem_cell(new TissueCell(p_healthy_state, p_cell_model));
        p_stem_cell->InitialiseCellCycleModel();

        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        // SimulationTime returns 30 hours

        // Create transit progeny of stem
        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);
        TissueCellPtr p_daughter_cell = p_stem_cell->Divide();

        std::vector<TissueCellPtr> cells;
        std::vector<TissueCellPtr> newly_born;

        // Track all the offspring of the daughter cell
        // after 3 generations they should become differentiated
        // and stop dividing
        cells.push_back(p_daughter_cell);

        std::vector<TissueCellPtr>::iterator cell_iterator;

        TS_ASSERT_EQUALS(p_params->GetMaxTransitGenerations(), 3u);
        unsigned int expected_num_cells[6];
        expected_num_cells[0] = 0;
        expected_num_cells[1] = 1;
        expected_num_cells[2] = 2;
        expected_num_cells[3] = 4;
        expected_num_cells[4] = 8;
        expected_num_cells[5] = 8;

        TS_ASSERT_EQUALS(expected_num_cells[1], cells.size());

        for (unsigned generation=2; generation<6; generation++)
        {
            // Produce the offspring of the cells
            cell_iterator = cells.begin();

            p_simulation_time->IncrementTimeOneStep();
            p_simulation_time->IncrementTimeOneStep();

            while (cell_iterator < cells.end())
            {
                if ((*cell_iterator)->ReadyToDivide())
                {
                    TissueCellPtr p_new_cell = (*cell_iterator)->Divide();
                    TS_ASSERT_DELTA((*cell_iterator)->GetAge(), 0, 1e-9);
                    if ((*cell_iterator)->GetCellCycleModel()->GetCellProliferativeType()==STEM)
                    {
                        TS_ASSERT_EQUALS(p_new_cell->GetCellCycleModel()->GetCellProliferativeType(), TRANSIT);
                        TS_ASSERT_EQUALS(static_cast<FixedDurationGenerationBasedCellCycleModel*>((*cell_iterator)->GetCellCycleModel())->GetGeneration(), 0u);
                    }
                    else if ((*cell_iterator)->GetCellCycleModel()->GetCellProliferativeType()==TRANSIT)
                    {
                        TS_ASSERT_EQUALS(p_new_cell->GetCellCycleModel()->GetCellProliferativeType(), TRANSIT);
                        TS_ASSERT_EQUALS(static_cast<FixedDurationGenerationBasedCellCycleModel*>((*cell_iterator)->GetCellCycleModel())->GetGeneration(), generation);
                    }
                    else
                    {
                        TS_ASSERT_EQUALS(p_new_cell->GetCellCycleModel()->GetCellProliferativeType(), DIFFERENTIATED);
                        TS_ASSERT_EQUALS(static_cast<FixedDurationGenerationBasedCellCycleModel*>((*cell_iterator)->GetCellCycleModel())->GetGeneration(), generation);
                    }

                    TS_ASSERT_DELTA(p_new_cell->GetAge(), 0, 1e-9);
                    TS_ASSERT_EQUALS(static_cast<FixedDurationGenerationBasedCellCycleModel*>(p_new_cell->GetCellCycleModel())->GetGeneration(), generation);

                    newly_born.push_back(p_new_cell);
                }
                ++cell_iterator;
            }

            // Copy offspring in newly_born vector to cells vector
            cell_iterator = newly_born.begin();
            while (cell_iterator < newly_born.end())
            {
                cells.push_back(*cell_iterator);
                ++cell_iterator;
            }
            newly_born.clear();

            // Check cell counts
            TS_ASSERT_EQUALS(expected_num_cells[generation], cells.size());
        }

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_EQUALS(cells[i]->GetCellCycleModel()->GetCellProliferativeType(), DIFFERENTIATED);
        }
    }

};

#endif /*TESTTISSUECELL_HPP_*/
