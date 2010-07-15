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
#ifndef TESTTISSUECELLNIGHTLY_HPP_
#define TESTTISSUECELLNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>
#include <iostream>

#include "TissueCellPtr.hpp"
#include "WntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "AbstractCellBasedTestSuite.hpp"

class TestTissueCellPtrNightly: public AbstractCellBasedTestSuite
{
public:

    void TestTysonNovakImmortalStemCell()
    {
        double end_time = 100.0; // A good load of divisions to make sure nothing mucks up..
        // one division = 1.26 hours.
        int time_steps = 1000;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, time_steps);
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        TysonNovakCellCycleModel* p_model = new TysonNovakCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_stem_cell(new TissueCell(p_healthy_state, p_model));
        p_stem_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);

        unsigned divisions = 0;
        while (p_simulation_time->GetTime() < end_time)
        {
            p_simulation_time->IncrementTimeOneStep();

            if (p_stem_cell->ReadyToDivide())
            {
                p_stem_cell->Divide();
                TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);
                divisions++;
            }
        }
        TS_ASSERT_DELTA(divisions, (unsigned)(end_time/1.26), 1);
    }

    void Test0DBucketWithTysonNovak()
    {
        double end_time = 7.0; // not very long because cell cycle time is only 1.2
        // (75 mins) because Tyson Novaks is for yeast
        int time_steps = 100;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, time_steps);
        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);

        TysonNovakCellCycleModel* p_model = new TysonNovakCellCycleModel();
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

        unsigned i = 0;
        while (p_simulation_time->GetTime()< end_time)
        {
            // Produce the offspring of the cells
            p_simulation_time->IncrementTimeOneStep();
            times[i] = p_simulation_time->GetTime();
            cell_iterator = cells.begin();
            unsigned j = 0;
            while (cell_iterator < cells.end())
            {
                if ((*cell_iterator)->ReadyToDivide())
                {
                    newly_born.push_back(cell_iterator->Divide());
                }
                ++cell_iterator;
                j++;
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

            i++;
        }

        TS_ASSERT_EQUALS(stem_cells[time_steps-1], 1u);
        TS_ASSERT_EQUALS(transit_cells[time_steps-1], 31u);
        TS_ASSERT_EQUALS(differentiated_cells[time_steps-1], 0u);
    }

    /*
     * We are checking that the TissueCellPtrs work with the Wnt cell cycle models here
     * That division of wnt cells and stuff works OK.
     *
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModelAndMutationAPCONEHIT() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        TissueConfig* p_parameters = TissueConfig::Instance();
        double s_g2_duration = p_parameters->GetSG2MDuration();

        unsigned num_steps = 200;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);

        double wnt_stimulus = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        p_cell_cycle_model1->SetCellProliferativeType(TRANSIT);
        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);
        TissueCellPtr p_wnt_cell(new TissueCell(p_apc1, p_cell_cycle_model1));
        p_wnt_cell->InitialiseCellCycleModel();

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            if (time >= 4.804 + s_g2_duration)
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

            bool result1 = p_wnt_cell->ReadyToDivide();
            bool result2 = p_wnt_cell2->ReadyToDivide();

            if (time >= 4.804 + s_g2_duration + time_of_birth)
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
     * We are checking that the TissueCellPtrs work with the Wnt cell cycle models here
     * That division of wnt cells and stuff works OK.
     *
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModelAndMutationBetaCat() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        TissueConfig* p_parameters = TissueConfig::Instance();
        double s_g2_duration = p_parameters->GetSG2MDuration();

        unsigned num_steps = 200;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);

        double wnt_stimulus = 0.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        p_cell_cycle_model1->SetCellProliferativeType(TRANSIT);
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);
        TissueCellPtr p_wnt_cell(new TissueCell(p_bcat1, p_cell_cycle_model1));
        p_wnt_cell->InitialiseCellCycleModel();

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            if (time >= 7.82 + s_g2_duration)
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

            bool result1 = p_wnt_cell->ReadyToDivide();
            bool result2 = p_wnt_cell2->ReadyToDivide();

            if (time >= 7.82 + s_g2_duration + time_of_birth)
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

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    /*
     * We are checking that the TissueCellPtrs work with the Wnt cell cycle models here
     * That division of wnt cells and stuff works OK.
     *
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModelAndMutationAPC2() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        TissueConfig* p_parameters = TissueConfig::Instance();

        double s_g2_duration = p_parameters->GetSG2MDuration();

        unsigned num_steps = 200;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);

        double wnt_stimulus = 0.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        p_cell_cycle_model1->SetCellProliferativeType(TRANSIT);
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        TissueCellPtr p_wnt_cell(new TissueCell(p_apc2, p_cell_cycle_model1));
        p_wnt_cell->InitialiseCellCycleModel();

        boost::shared_ptr<AbstractCellMutationState> p_this_state = p_wnt_cell->GetMutationState();

        TS_ASSERT_EQUALS(p_this_state->IsType<ApcTwoHitCellMutationState>(), true);

        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);

        p_wnt_cell->SetMutationState(p_apc1);

        p_this_state = p_wnt_cell->GetMutationState();

        TS_ASSERT_EQUALS(p_this_state->IsType<ApcTwoHitCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_this_state->IsType<ApcOneHitCellMutationState>(), true);

        p_wnt_cell->SetMutationState(p_apc2);

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            if (time >= 3.9435 + s_g2_duration)
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

            bool result1 = p_wnt_cell->ReadyToDivide();
            bool result2 = p_wnt_cell2->ReadyToDivide();

            if (time >= 3.9435 + s_g2_duration + time_of_birth)
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

        // Tidy up
        WntConcentration<2>::Destroy();
    }
};


#endif /*TESTTISSUECELLNIGHTLY_HPP_*/
