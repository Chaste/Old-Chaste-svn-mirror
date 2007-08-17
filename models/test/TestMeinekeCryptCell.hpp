#ifndef TESTMEINEKECRYPTCELL_HPP_
#define TESTMEINEKECRYPTCELL_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "OutputFileHandler.hpp"
#include "MeinekeCryptCellTypes.hpp"
#include "CryptCellMutationStates.hpp"
#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "SimulationTime.hpp"
#include "WntGradient.hpp"
#include <iostream>

class TestMeinekeCryptCell: public CxxTest::TestSuite
{
public:

    void TestCellsAgeingCorrectly() throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        
        // These lines are added to cover the exception case that a cell is
        // created without simulation time being set up...
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        FixedCellCycleModel fixed_model;
        SimulationTime::Destroy();
        
        TS_ASSERT_THROWS_ANYTHING(MeinekeCryptCell bad_cell(STEM, // type
                                                            HEALTHY,//Mutation State
                                                            0,    // generation
                                                            &fixed_model));
                                                            
        // Proper test again
        
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);
        
        MeinekeCryptCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,    // generation
                                   new FixedCellCycleModel());
                                   
        p_simulation_time->IncrementTimeOneStep();
        
        TS_ASSERT_EQUALS(stem_cell.GetAge(), 0.5);
        
        //for coverage
        stem_cell.SetNodeIndex(3);
        TS_ASSERT_EQUALS((int)(stem_cell.GetNodeIndex()), 3);
        
        p_simulation_time->IncrementTimeOneStep();
        stem_cell.SetBirthTime(p_simulation_time->GetDimensionalisedTime());
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(stem_cell.GetAge(), 1.0);
        
        TS_ASSERT_EQUALS(stem_cell.IsDead(), false);
        stem_cell.Kill();
        TS_ASSERT_EQUALS(stem_cell.IsDead(), true);
        
        SimulationTime::Destroy();
    }
    
    void TestCellDivision()
    {
        CancerParameters::Instance()->Reset();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(54.0, 9);
        // We are going to start at t=0 and jump up in steps of 6.0
        CancerParameters *p_params = CancerParameters::Instance();
        
        // this test needs particular cell cycle times
        TS_ASSERT_EQUALS(p_params->GetStemCellCycleTime(), 24.0);
        TS_ASSERT_EQUALS(p_params->GetTransitCellCycleTime(), 12.0);
        
        p_simulation_time->IncrementTimeOneStep();//t=6
        MeinekeCryptCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,    // generation
                                   new FixedCellCycleModel());
        p_simulation_time->IncrementTimeOneStep();//t=12
        p_simulation_time->IncrementTimeOneStep();//t=18
        p_simulation_time->IncrementTimeOneStep();//t=24
        TS_ASSERT(!stem_cell.ReadyToDivide());
        
        
        p_simulation_time->IncrementTimeOneStep();//t=30
        TS_ASSERT(stem_cell.ReadyToDivide());
        
        // create transit progeny of stem
        MeinekeCryptCell daughter_cell = stem_cell.Divide();
        
        TS_ASSERT(!stem_cell.ReadyToDivide());
        
        TS_ASSERT(daughter_cell.GetGeneration() == 1);
        TS_ASSERT(daughter_cell.GetCellType() == TRANSIT);
        TS_ASSERT_DELTA(daughter_cell.GetAge(), 0 , 1e-9);
        
        p_simulation_time->IncrementTimeOneStep();//t=36
        TS_ASSERT(!daughter_cell.ReadyToDivide());
        p_simulation_time->IncrementTimeOneStep();//t=42
        TS_ASSERT(daughter_cell.ReadyToDivide());
        
        // create transit progeny of transit
        MeinekeCryptCell grandaughter_cell = daughter_cell.Divide();
        
        p_simulation_time->IncrementTimeOneStep();//t=48
        TS_ASSERT(!stem_cell.ReadyToDivide());
        
        TS_ASSERT(!grandaughter_cell.ReadyToDivide());
        TS_ASSERT(!daughter_cell.ReadyToDivide());
        
        // stem cell ready to divide again
        p_simulation_time->IncrementTimeOneStep();//t=54
        TS_ASSERT(stem_cell.ReadyToDivide());
        
        // both grandaughter and daughter cells should be ready to
        // divide
        TS_ASSERT(grandaughter_cell.ReadyToDivide());
        TS_ASSERT(daughter_cell.ReadyToDivide());
        
        SimulationTime::Destroy();
    }
    
    void TestCellDivisionStops()
    {
        CancerParameters::Instance()->Reset();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(54.0, 9);
        CancerParameters *p_params = CancerParameters::Instance();
        
        // If the value of GetStemCellCycleTime() changes in p_params the simulation time
        // step and end time will need to be changed accordingly so that
        // IncrementTimeOneStep() gets the cell to correct division times
        
        TS_ASSERT_EQUALS(p_params->GetStemCellCycleTime(), 24.0);
        TS_ASSERT_EQUALS(p_params->GetTransitCellCycleTime(), 12.0);
        
        // SimulationTime returns 0 hours
        p_simulation_time->IncrementTimeOneStep();
        // SimulationTime returns 6 hours
        
        MeinekeCryptCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,    // generation
                                   new FixedCellCycleModel());
                                   
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();
        
        // SimulationTime returns 30 hours
        
        // create transit progeny of stem
        TS_ASSERT(stem_cell.ReadyToDivide());
        MeinekeCryptCell daughter_cell = stem_cell.Divide();
        
        std::vector<MeinekeCryptCell> cells;
        std::vector<MeinekeCryptCell> newly_born;
        
        // track all the offspring of the daughter cell
        // after 3 generations they should become differentiated
        // and stop dividing
        cells.push_back(daughter_cell);
        
        std::vector<MeinekeCryptCell>::iterator cell_iterator;
        
        TS_ASSERT_EQUALS(p_params->GetMaxTransitGenerations(), 3u);
        unsigned int expected_num_cells[6];
        expected_num_cells[0]=0;
        expected_num_cells[1]=1;
        expected_num_cells[2]=2;
        expected_num_cells[3]=4;
        expected_num_cells[4]=8;
        expected_num_cells[5]=8;
        
        TS_ASSERT_EQUALS(expected_num_cells[1], cells.size());
        
        for (int generation=2; generation<6; generation++)
        {
            // produce the offspring of the cells
            cell_iterator = cells.begin();
            
            
            p_simulation_time->IncrementTimeOneStep();
            p_simulation_time->IncrementTimeOneStep();
            
            
            while (cell_iterator < cells.end())
            {
                if (cell_iterator->ReadyToDivide())
                {
                    newly_born.push_back(cell_iterator->Divide());
                }
                cell_iterator++;
            }
            
            // copy offspring in newly_born vector to cells vector
            cell_iterator = newly_born.begin();
            while (cell_iterator < newly_born.end())
            {
                cells.push_back(*cell_iterator);
                cell_iterator++;
            }
            newly_born.clear();
            
            // check cell counts
            TS_ASSERT_EQUALS(expected_num_cells[generation], cells.size());
        }
        
        for (unsigned i=0; i<cells.size() ; i++)
        {
            TS_ASSERT_EQUALS(cells[i].GetCellType(), DIFFERENTIATED);
        }
        SimulationTime::Destroy();
    }
    
    void TestUpdateCellTypes() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(200, 20);
        
        MeinekeCryptCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,  // generation
                                   new FixedCellCycleModel());
                                   
        stem_cell.UpdateCellType();
        
        TS_ASSERT_EQUALS(stem_cell.GetCellType(),STEM);
        
        stem_cell.SetCellType(TRANSIT);
        
        stem_cell.UpdateCellType();
        
        TS_ASSERT_EQUALS(stem_cell.GetCellType(),TRANSIT);
        
        // Test a Wnt dependent cell
        WntGradient wnt_gradient;
        
        MeinekeCryptCell wnt_cell(TRANSIT, // type
                                   HEALTHY,//Mutation State
                                   0,  // generation
                                   new WntCellCycleModel(0.0, wnt_gradient));
                                   
        wnt_cell.UpdateCellType();
        
        TS_ASSERT_EQUALS(wnt_cell.GetCellType(),DIFFERENTIATED);
        
        // Push up the Wnt stimulus and go forward through time
        std::vector<double> cell_cycle_influences;
        cell_cycle_influences.push_back(1.0);
        for (unsigned i=0 ; i<20 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }
        wnt_cell.ReadyToDivide(cell_cycle_influences);
        wnt_cell.UpdateCellType();
        
        TS_ASSERT_EQUALS(wnt_cell.GetCellType(),TRANSIT);
          
        SimulationTime::Destroy();
    }
    
    void Test0DBucket()
    {
        CancerParameters::Instance()->Reset();
        
        double end_time=61.0;
        int time_steps=61;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, time_steps);
        
        MeinekeCryptCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,  // generation
                                   new FixedCellCycleModel());
                                   
                                   
                                   
                                   
        std::vector<MeinekeCryptCell> cells;
        std::vector<MeinekeCryptCell> newly_born;
        std::vector<unsigned> stem_cells(time_steps);
        std::vector<unsigned> transit_cells(time_steps);
        std::vector<unsigned> differentiated_cells(time_steps);
        std::vector<double> times(time_steps);
        
        cells.push_back(stem_cell);
        std::vector<MeinekeCryptCell>::iterator cell_iterator;
        
        unsigned i=0;
        while (p_simulation_time->GetDimensionalisedTime()< end_time)
        {
            // produce the offspring of the cells
            
            p_simulation_time->IncrementTimeOneStep();
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                if (cell_iterator->ReadyToDivide())
                {
                    newly_born.push_back(cell_iterator->Divide());
                }
                cell_iterator++;
            }
            
            // copy offspring in newly_born vector to cells vector
            cell_iterator = newly_born.begin();
            while (cell_iterator < newly_born.end())
            {
                cells.push_back(*cell_iterator);
                cell_iterator++;
            }
            newly_born.clear();
            
            // update cell counts
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                switch (cell_iterator->GetCellType())
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
                
                cell_iterator++;
            }
            times[i]=p_simulation_time->GetDimensionalisedTime();
            i++;
        }
        TS_ASSERT_EQUALS(stem_cells[59], 1u);
        TS_ASSERT_EQUALS(transit_cells[59], 2u);
        TS_ASSERT_EQUALS(differentiated_cells[59], 8u);
        SimulationTime::Destroy();
    }
    
    void TestWithFixedCellCycleModel() throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        // Simulation time is 6000 because we want to test that differentiated cells never divide.
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(6000.0, 1000);
        CancerParameters *p_params = CancerParameters::Instance();
        
        // this test needs particular cell cycle times
        TS_ASSERT_EQUALS(p_params->GetStemCellCycleTime(), 24.0);
        TS_ASSERT_EQUALS(p_params->GetTransitCellCycleTime(), 12.0);
        
        RandomNumberGenerator::Instance();
        
        p_simulation_time->IncrementTimeOneStep();
        
        
        //  Creating different types of Meineke crypt cells with different cell cycle models at SImulationTime = 6 hours.
        MeinekeCryptCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,    // generation
                                   new FixedCellCycleModel());
                                   
        MeinekeCryptCell stochastic_stem_cell(STEM, // type
                                              HEALTHY,//Mutation State
                                              0,    // generation
                                              new StochasticCellCycleModel);
        MeinekeCryptCell differentiated_cell(DIFFERENTIATED, // type
                                             HEALTHY,//Mutation State
                                             6,    // generation
                                             new FixedCellCycleModel());
        MeinekeCryptCell stochastic_differentiated_cell(DIFFERENTIATED, // type
                                                        HEALTHY,//Mutation State
                                                        6,    // generation
                                                        new StochasticCellCycleModel);
        MeinekeCryptCell transit_cell(TRANSIT, // type
                                      HEALTHY,//Mutation State
                                      2,    // generation
                                      new FixedCellCycleModel());
                                      
                                      
        // SimulationTime = 6 hours
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();
        
        // SimulationTime = 18 hours
        
        TS_ASSERT(transit_cell.ReadyToDivide());
        
        p_simulation_time->IncrementTimeOneStep();
        
        // SimulationTime = 24 hours
        
        TS_ASSERT(!stem_cell.ReadyToDivide());
        TS_ASSERT(!stochastic_stem_cell.ReadyToDivide());
        
        p_simulation_time->IncrementTimeOneStep();
        
        // SimulationTime = 30 hours
        
        TS_ASSERT(stem_cell.ReadyToDivide());
        TS_ASSERT(stochastic_stem_cell.ReadyToDivide());
        
        MeinekeCryptCell daughter_cell1 = stem_cell.Divide();
        TS_ASSERT(typeid(daughter_cell1.GetCellCycleModel()) == typeid(stem_cell.GetCellCycleModel()));
        
        // Go to large time to ensure that differentiated cells can not divide
        
        for (int i=0; i<990; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }
        TS_ASSERT(!differentiated_cell.ReadyToDivide());
        TS_ASSERT(!stochastic_differentiated_cell.ReadyToDivide());
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        
    }
    
    
    
    void TestStochasticCycleModel() throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        
        // Go up in steps of 0.01 to test stochasticity in cell cycle models
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(54.0, 5400);
        CancerParameters *p_params = CancerParameters::Instance();
        
        // this test needs particular cell cycle times
        TS_ASSERT_EQUALS(p_params->GetStemCellCycleTime(), 24.0);
        TS_ASSERT_EQUALS(p_params->GetTransitCellCycleTime(), 12.0);
        
        RandomNumberGenerator::Instance();
        
        for (int i=0; i<600; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }
        // now at t=6.00
        MeinekeCryptCell transit_cell(TRANSIT, // type
                                      HEALTHY,//Mutation State
                                      2,    // generation
                                      new FixedCellCycleModel());
                                      
        for (int i=0; i<1199; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }
        // now at t = 17.99, cell is 11.99 old
        TS_ASSERT(!transit_cell.ReadyToDivide());
        
        StochasticCellCycleModel *cell_cycle_model = new StochasticCellCycleModel;
        // this now re-sets the age of the cell to 0.0 so more time added in underneath
        transit_cell.SetCellCycleModel(cell_cycle_model);
        TS_ASSERT_EQUALS(transit_cell.GetCellCycleModel(), cell_cycle_model);
        for (int i=0; i<1199; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }
        int ready_count=0;
        for (int i=0; i<100; i++)
        {
            //std::cout << "time = " << time << " transit cell age = " << transit_cell.GetAge() << "\n";
            if (transit_cell.ReadyToDivide())
            {
                ready_count++;
            }
        }
        TS_ASSERT(ready_count>0);
        
        // Ensure transit cell divides
        while (!transit_cell.ReadyToDivide())
        {
            p_simulation_time->IncrementTimeOneStep();
        }
        MeinekeCryptCell daughter_cell2 = transit_cell.Divide();
        TS_ASSERT(typeid(daughter_cell2.GetCellCycleModel()) == typeid(transit_cell.GetCellCycleModel()));
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    
    
    void Test0DBucketStochastic()
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        
        // this test needs particular cell cycle times
        TS_ASSERT_EQUALS(p_params->GetStemCellCycleTime(), 24.0);
        TS_ASSERT_EQUALS(p_params->GetTransitCellCycleTime(), 12.0);
        
        RandomNumberGenerator::Instance();
        
        const double end_time = 70.0;
        //const int time_steps = 70;
        const int number_of_simulations = 1000;
        //const double time_step= end_time/(double) time_steps;
        
        std::vector<MeinekeCryptCell> cells;
        std::vector<MeinekeCryptCell> newly_born;
        
        std::vector<unsigned> stem_cells(number_of_simulations);
        std::vector<unsigned> transit_cells(number_of_simulations);
        std::vector<unsigned> differentiated_cells(number_of_simulations);
        double stem_cell_mean = 0.0;
        double transit_cell_mean = 0.0;
        double differentiated_cell_mean = 0.0;
        
        for (int simulation_number=0; simulation_number<number_of_simulations; simulation_number++)
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(70.0, 70);
            
            MeinekeCryptCell stem_cell(STEM, // type
                                       HEALTHY,//Mutation State
                                       0,  // generation
                                       new StochasticCellCycleModel);
            cells.push_back(stem_cell);
            // produce the offspring of the cells
            std::vector<MeinekeCryptCell>::iterator cell_iterator = cells.begin();
            
            
            
            while (p_simulation_time->GetDimensionalisedTime()< end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                cell_iterator = cells.begin();
                while (cell_iterator < cells.end())
                {
                    if (cell_iterator->ReadyToDivide())
                    {
                        newly_born.push_back(cell_iterator->Divide());
                    }
                    cell_iterator++;
                }
                
                // copy offspring in newly_born vector to cells vector
                cell_iterator = newly_born.begin();
                while (cell_iterator < newly_born.end())
                {
                    cells.push_back(*cell_iterator);
                    cell_iterator++;
                }
                
                newly_born.clear();
            }
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                switch (cell_iterator->GetCellType())
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
                
                cell_iterator++;
            }
            cells.clear();
            SimulationTime::Destroy();
        }
        stem_cell_mean=stem_cell_mean/(double) number_of_simulations;
        transit_cell_mean=transit_cell_mean/(double) number_of_simulations;
        differentiated_cell_mean=differentiated_cell_mean/(double) number_of_simulations;
        
        TS_ASSERT_DELTA(stem_cell_mean, 1.0, 1e-12);
        TS_ASSERT_DELTA(transit_cell_mean, 2.0, 1.0);
        TS_ASSERT_DELTA(differentiated_cell_mean, 8.0, 1.0);
        
        
        TS_ASSERT_DELTA(p_params->GetTransitCellCycleTime(), 12.0, 1e-12);
        
        RandomNumberGenerator::Destroy();
    }
    
    /* We are setting up a 0d bucket with some initial cell population
     * This is deterministic so we can test it
     */
    void TestInitialise0DBucket()
    {
        CancerParameters::Instance()->Reset();
        
        //double end_time=60.0;
        //int time_steps=60;
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(60.0, 60);
        
        std::vector<MeinekeCryptCell> cells;
        
        MeinekeCryptCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,  // generation
                                   new FixedCellCycleModel());
                                   
        cells.push_back(stem_cell);
        
        MeinekeCryptCell transit_cell_1(TRANSIT, // type
                                        HEALTHY,//Mutation State
                                        1,  // generation
                                        new FixedCellCycleModel());
                                        
        cells.push_back(transit_cell_1);
        
        MeinekeCryptCell transit_cell_2(TRANSIT, // type
                                        HEALTHY,//Mutation State
                                        2,  // generation
                                        new FixedCellCycleModel());
                                        
        cells.push_back(transit_cell_2);
        
        MeinekeCryptCell transit_cell_3(TRANSIT, // type
                                        HEALTHY,//Mutation State
                                        3,  // generation
                                        new FixedCellCycleModel());
                                        
        cells.push_back(transit_cell_3);
        
        MeinekeCryptCell differentiated_cell(DIFFERENTIATED, // type
                                             HEALTHY,//Mutation State
                                             4,  // generation
                                             new FixedCellCycleModel());
                                             
        cells.push_back(differentiated_cell);
        
        //double time=0.0;
        
        //double time_step= end_time/(double) time_steps;
        
        std::vector<MeinekeCryptCell> newly_born;
        std::vector<unsigned> stem_cells(p_simulation_time->GetTotalNumberOfTimeSteps());
        std::vector<unsigned> transit_cells(p_simulation_time->GetTotalNumberOfTimeSteps());
        std::vector<unsigned> differentiated_cells(p_simulation_time->GetTotalNumberOfTimeSteps());
        std::vector<double> times(p_simulation_time->GetTotalNumberOfTimeSteps());
        
        std::vector<MeinekeCryptCell>::iterator cell_iterator;
        
        unsigned i=0;
        while (!p_simulation_time->IsFinished())
        {
            p_simulation_time->IncrementTimeOneStep();
            // produce the offspring of the cells
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                CryptCellMutationState this_cell_state;
                this_cell_state = cell_iterator->GetMutationState();
                TS_ASSERT(this_cell_state==HEALTHY);
                if (cell_iterator->ReadyToDivide())
                {
                    newly_born.push_back(cell_iterator->Divide());
                }
                cell_iterator++;
            }
            
            // copy offspring in newly_born vector to cells vector
            cell_iterator = newly_born.begin();
            while (cell_iterator < newly_born.end())
            {
                cells.push_back(*cell_iterator);
                cell_iterator++;
            }
            newly_born.clear();
            
            // count # cells of each type
            cell_iterator = cells.begin();
            stem_cells[i] = 0;
            transit_cells[i] = 0;
            differentiated_cells[i] = 0;
            while (cell_iterator < cells.end())
            {
                switch (cell_iterator->GetCellType())
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
                cell_iterator++;
            }
            
            times[i]=p_simulation_time->GetDimensionalisedTime();
            i++;
        }
        
        TS_ASSERT_EQUALS(stem_cells[59], 1u);
        TS_ASSERT_EQUALS(transit_cells[59], 2u);
        TS_ASSERT_EQUALS(differentiated_cells[59], 23u);
        
        SimulationTime::Destroy();
    }
    
    
    /*
     * We are checking that the MeinekeCryptCells work with the Wnt cell cycle models here
     * That division of wnt cells and stuff works OK.
     * 
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModel() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        CancerParameters *p_parameters = CancerParameters::Instance();
        p_parameters->Reset();
        
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        unsigned num_steps=100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);
        
        WntGradient wnt_gradient;
        
        double wnt_stimulus = 1.0;
        std::cout << "Hello Joe 1\n" << std::flush;
        MeinekeCryptCell wnt_cell(TRANSIT, // type
                                  HEALTHY,//Mutation State
                                  1,    // generation
                                  new WntCellCycleModel(wnt_stimulus, wnt_gradient));
        std::cout << "Hello Joe 2\n" << std::flush;
                 
        for (unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            std::vector<double> wnt;
            wnt.push_back(wnt_stimulus);
            if (time>=5.971+SG2MDuration)
            {
                TS_ASSERT(wnt_cell.ReadyToDivide(wnt)==true);
            }
            else
            {
                TS_ASSERT(wnt_cell.ReadyToDivide(wnt)==false);
            }
            //std::cout << "Time = " << time << " ready = " << wnt_cell.ReadyToDivide(wnt) << "\n" << std::endl;
        }
        
        std::vector<double> wnt;
        wnt.push_back(wnt_stimulus);
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT(wnt_cell.ReadyToDivide(wnt)==true);
        TS_ASSERT(wnt_cell.GetGeneration()==1);
        
        std::cout << "Hello Joe 3\n" << std::flush;
        MeinekeCryptCell wnt_cell2 = wnt_cell.Divide();
        std::cout << "Hello Joe 4\n" << std::flush;
        TS_ASSERT(wnt_cell.GetGeneration()==2);
        TS_ASSERT(wnt_cell2.GetGeneration()==2);
        
        //std::cout << "time now = " << p_simulation_time->GetDimensionalisedTime() << "\n" <<std::endl;
        
        double timeOfBirth = wnt_cell.GetBirthTime();
        double timeOfBirth2 = wnt_cell2.GetBirthTime();
        
        TS_ASSERT_DELTA(timeOfBirth, timeOfBirth2, 1e-9);
        
        //std::cout << "time of cell divisions = " << timeOfBirth << "\tand\t" << timeOfBirth2 << "\n" << std::endl;
        
        for (unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            std::vector<double> wnt;
            wnt.push_back(wnt_stimulus);
            bool result1=wnt_cell.ReadyToDivide(wnt);
            bool result2=wnt_cell2.ReadyToDivide(wnt);
            //std::cout << "Time = " << time << ", ready1 = " << result1 << ", ready2 = " << result2<< "\n" << std::endl;
            if (time>=5.971+SG2MDuration+timeOfBirth)
            {
                TS_ASSERT(result1==true);
                TS_ASSERT(result2==true);
            }
            else
            {
                TS_ASSERT(result1==false);
                TS_ASSERT(result2==false);
            }
        }
        
        SimulationTime::Destroy();
    }
    
    /*
     * We are checking that the MeinekeCryptCells work with the StochasticWnt cell cycle models here
     * That division of wnt cells and stuff works OK.
     * 
     * It checks that the cell division thing works nicely too.
     */
    void TestWithStochasticWntCellCycleModel() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        RandomNumberGenerator::Instance()->Reseed(0);
        CancerParameters *p_parameters = CancerParameters::Instance();
        p_parameters->Reset();
        
        // these are the first three normal random with mean 10, s.d. 1 and this seed (0)
        double SG2MDuration1 = 9.0676;
        double SG2MDuration2 = 11.1632;
        double SG2MDuration3 = 9.2712;
        
        unsigned num_steps=100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);
        
        double wnt_stimulus = 1.0;
        WntGradient wnt_gradient;
        MeinekeCryptCell wnt_cell(TRANSIT, // type
                                  HEALTHY,//Mutation State
                                  1,    // generation
                                  new StochasticWntCellCycleModel(wnt_stimulus, wnt_gradient));
                                  
        for (unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            std::vector<double> wnt;
            wnt.push_back(wnt_stimulus);
            if (time>=5.971+SG2MDuration1)
            {
                TS_ASSERT(wnt_cell.ReadyToDivide(wnt)==true);
            }
            else
            {
                TS_ASSERT(wnt_cell.ReadyToDivide(wnt)==false);
            }
            //std::cout << "Time = " << time << " ready = " << wnt_cell.ReadyToDivide(wnt) << "\n" << std::endl;
        }
        
        std::vector<double> wnt;
        wnt.push_back(wnt_stimulus);
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT(wnt_cell.ReadyToDivide(wnt)==true);
        TS_ASSERT(wnt_cell.GetGeneration()==1);
        
        MeinekeCryptCell wnt_cell2 = wnt_cell.Divide();
        
        TS_ASSERT(wnt_cell.GetGeneration()==2);
        TS_ASSERT(wnt_cell2.GetGeneration()==2);
        
        //std::cout << "time now = " << p_simulation_time->GetDimensionalisedTime() << "\n" <<std::endl;
        
        double timeOfBirth = wnt_cell.GetBirthTime();
        double timeOfBirth2 = wnt_cell2.GetBirthTime();
        
        TS_ASSERT_DELTA(timeOfBirth, timeOfBirth2, 1e-9);
        
        //std::cout << "time of cell divisions = " << timeOfBirth << "\tand\t" << timeOfBirth2 << "\n" << std::endl;
        
        for (unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            std::vector<double> wnt;
            wnt.push_back(wnt_stimulus);
            bool result1=wnt_cell.ReadyToDivide(wnt);
            bool result2=wnt_cell2.ReadyToDivide(wnt);
            //std::cout << "Time = " << time << ", ready1 = " << result1 << ", ready2 = " << result2<< "\n" << std::endl;
            if (time>=5.971+SG2MDuration2+timeOfBirth)
            {
                TS_ASSERT(result1==true);
            }
            else
            {
                TS_ASSERT(result1==false);
            }
            if (time>=5.971+SG2MDuration3+timeOfBirth)
            {
                TS_ASSERT(result2==true);
            }
            else
            {
                TS_ASSERT(result2==false);
            }
            
        }
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    /*
     * We are checking that the MeinekeCryptCells work with the T&N cell cycle models here
     * That division of wnt cells and stuff works OK.
     * 
     * It checks that the cell division thing works nicely too.
     */
    void TestWithTysonNovakCellCycleModel() throw(Exception)
    {
        CancerParameters::Instance()->Reset();

        double standard_tyson_duration = 75.19/60.0;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        unsigned num_steps=100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(200.0/60.0, num_steps+1);
        
        MeinekeCryptCell tn_cell(TRANSIT, // type
                                 HEALTHY,//Mutation State
                                 1,    // generation
                                 new TysonNovakCellCycleModel());
                                 
        for (unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            if (time>=standard_tyson_duration)
            {
                TS_ASSERT(tn_cell.ReadyToDivide()==true);
            }
            else
            {
                TS_ASSERT(tn_cell.ReadyToDivide()==false);
            }
            //std::cout << "Time = " << time << " ready = " << tn_cell.ReadyToDivide(wnt) << "\n" << std::endl;
        }
        
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT(tn_cell.ReadyToDivide()==true);
        TS_ASSERT(tn_cell.GetGeneration()==1);
        
        MeinekeCryptCell tn_cell2 = tn_cell.Divide();
        
        TS_ASSERT(tn_cell.GetGeneration()==2);
        TS_ASSERT(tn_cell2.GetGeneration()==2);
        
        //std::cout << "time now = " << p_simulation_time->GetDimensionalisedTime() << "\n" <<std::endl;
        
        double timeOfBirth = tn_cell.GetBirthTime();
        double timeOfBirth2 = tn_cell2.GetBirthTime();
        
        TS_ASSERT_DELTA(timeOfBirth, timeOfBirth2, 1e-9);
        
        //std::cout << "time of cell divisions = " << timeOfBirth << "\tand\t" << timeOfBirth2 << "\n" << std::endl;
        
        for (unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            bool result1=tn_cell.ReadyToDivide();
            bool result2=tn_cell2.ReadyToDivide();
            //std::cout << "Time = " << time << ", ready1 = " << result1 << ", ready2 = " << result2<< "\n" << std::endl;
            if (time>=standard_tyson_duration+timeOfBirth)
            {
                TS_ASSERT(result1==true);
                TS_ASSERT(result2==true);
            }
            else
            {
                TS_ASSERT(result1==false);
                TS_ASSERT(result2==false);
            }
        }
        
        SimulationTime::Destroy();
    }
    
    void TestApoptosisAndDeath()
    {
        CancerParameters::Instance()->Reset();
        
        // We are going to start at t=0 and jump up in steps of 0.2
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(0.6, 3);
        
        // this test needs particular apoptosis time
        CancerParameters *p_params = CancerParameters::Instance();
        TS_ASSERT_EQUALS(p_params->GetApoptosisTime(), 0.25);
        
        
        MeinekeCryptCell cell(TRANSIT, // type
                              HEALTHY,//Mutation State
                              0,    // generation
                              new FixedCellCycleModel());
                              
        TS_ASSERT_EQUALS(cell.HasApoptosisBegun(),false);
        TS_ASSERT_EQUALS(cell.IsDead(),false);
        TS_ASSERT_THROWS_ANYTHING(cell.TimeUntilDeath());
        
        p_simulation_time->IncrementTimeOneStep();//t=0.2
        
        cell.StartApoptosis();
        TS_ASSERT_THROWS_ANYTHING(cell.StartApoptosis());
        
        TS_ASSERT_EQUALS(cell.HasApoptosisBegun(),true);
        TS_ASSERT_EQUALS(cell.IsDead(),false);
        TS_ASSERT_DELTA(cell.TimeUntilDeath(),0.25,1e-12);
        
        // check that we can copy a cell that has started apoptosis
        MeinekeCryptCell cell2(cell);
        
        p_simulation_time->IncrementTimeOneStep();//t=0.4
        TS_ASSERT_EQUALS(cell.HasApoptosisBegun(),true);
        TS_ASSERT_EQUALS(cell.IsDead(),false);
        TS_ASSERT_DELTA(cell.TimeUntilDeath(),0.05,1e-12);
        
        TS_ASSERT_EQUALS(cell2.HasApoptosisBegun(),true);
        TS_ASSERT_EQUALS(cell2.IsDead(),false);
        TS_ASSERT_DELTA(cell2.TimeUntilDeath(),0.05,1e-12);
        
        p_simulation_time->IncrementTimeOneStep();//t=0.6
        TS_ASSERT_EQUALS(cell.HasApoptosisBegun(),true);
        TS_ASSERT_EQUALS(cell.IsDead(),true);
        
        SimulationTime::Destroy();
    }
    
    
    void TestCantDivideIfUndergoingApoptosis()
    {
        CancerParameters::Instance()->Reset();
        
        // We are going to start at t=0 and jump up to t=25
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 1);
        
        MeinekeCryptCell cell(TRANSIT, // type
                              HEALTHY,//Mutation State
                              0,    // generation
                              new FixedCellCycleModel());
                              
        p_simulation_time->IncrementTimeOneStep();//t=25
        
        TS_ASSERT_EQUALS(cell.ReadyToDivide(), true);
        cell.StartApoptosis();
        TS_ASSERT_EQUALS(cell.ReadyToDivide(), false);
        SimulationTime::Destroy();
    }
    
    
    
    void Test0DBucketWithDeath()
    {
        CancerParameters::Instance()->Reset();

        double end_time=92.0;
        int time_steps=92;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, time_steps);
        
        MeinekeCryptCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,  // generation
                                   new FixedCellCycleModel());
                                   
                                   
                                   
        std::vector<MeinekeCryptCell> cells;
        std::vector<MeinekeCryptCell> newly_born;
        std::vector<unsigned> stem_cells(time_steps);
        std::vector<unsigned> transit_cells(time_steps);
        std::vector<unsigned> differentiated_cells(time_steps);
        std::vector<unsigned> dead_cells(time_steps);
        std::vector<double> times(time_steps);
        
        cells.push_back(stem_cell);
        std::vector<MeinekeCryptCell>::iterator cell_iterator;
        
        unsigned i=0;
        while (p_simulation_time->GetDimensionalisedTime()< end_time)
        {
            // produce the offspring of the cells
            
            p_simulation_time->IncrementTimeOneStep();
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                if (!cell_iterator->IsDead())
                {
                    if (cell_iterator->ReadyToDivide())
                    {
                        newly_born.push_back(cell_iterator->Divide());
                    }
                    
                    if ((cell_iterator->GetAge() > 30))
                    {
                        cell_iterator->StartApoptosis();
                    }
                }
                
                
                cell_iterator++;
            }
            
            // copy offspring in newly_born vector to cells vector
            cell_iterator = newly_born.begin();
            while (cell_iterator < newly_born.end())
            {
                cells.push_back(*cell_iterator);
                cell_iterator++;
            }
            newly_born.clear();
            
            // update cell counts
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                if (!cell_iterator->IsDead())
                {
                    switch (cell_iterator->GetCellType())
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
                
                cell_iterator++;
            }
            times[i]=p_simulation_time->GetDimensionalisedTime();
            i++;
        }
        
        
        TS_ASSERT_EQUALS(stem_cells[time_steps-1], 1u);
        TS_ASSERT_EQUALS(transit_cells[time_steps-1], 2u);
        TS_ASSERT_EQUALS(differentiated_cells[time_steps-1], 8u);
        TS_ASSERT_EQUALS(dead_cells[time_steps-1], 8u);
        
        SimulationTime::Destroy();
    }
    
    
    
    void TestArchiveCell() throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "cell.arch";
        
        // Archive a Meineke Crypt cell
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);
            
            MeinekeCryptCell stem_cell(STEM, // type
                                       HEALTHY,//Mutation State
                                       0,    // generation
                                       new FixedCellCycleModel());
                                       
            p_simulation_time->IncrementTimeOneStep();
            
            TS_ASSERT_EQUALS(stem_cell.GetAge(), 0.5);
            
            stem_cell.SetNodeIndex(3);
            
            // Create an ouput archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            // and write the cell to the archive
            output_arch << static_cast<const SimulationTime&> (*p_simulation_time);
            output_arch << static_cast<const MeinekeCryptCell&>(stem_cell);
            SimulationTime::Destroy();
        }
        
        // Restore Meineke Crypt Cell
        {
            // need to set up time to initialise a cell
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(1.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 1); // will be restored
            
            // Initialise a cell
            AbstractCellCycleModel *temp_model = new FixedCellCycleModel();
            MeinekeCryptCell stem_cell(TRANSIT, // the type will be restored soon
                                       HEALTHY,//Mutation State
                                       1,    // generation
                                       temp_model);
                                       
            // restore the cell
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            input_arch >> *p_simulation_time;
            input_arch >> stem_cell;
            
            // Free memory
            delete temp_model;
            
            // check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetDimensionalisedTime(), 0.5);
            TS_ASSERT_EQUALS(p_simulation_time->GetTimeStep(), 0.5);
            
            TS_ASSERT_EQUALS(stem_cell.GetNodeIndex(), 3u);
            TS_ASSERT_EQUALS(stem_cell.GetAge(), 0.5);
            TS_ASSERT_EQUALS(stem_cell.GetGeneration(), 0u);
            TS_ASSERT_EQUALS(stem_cell.GetCellType(), STEM);
            
            SimulationTime::Destroy();
        }
    }
    
    /*
     * We are checking that the MeinekeCryptCells work with the Wnt 
     * cell cycle models here. This just tests the set-up and checks that
     * the functions can all be called (not what they return). 
     * 
     * For more in depth tests see TestNightlyMeinekeCryptCell.hpp
     * (these test that the cell cycle times are correct for the 
     * various mutant cells)
     */
    void TestWntMutantVariantsAndLabelling() throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        unsigned num_steps=10;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
        
        double wnt_stimulus = 1.0;
        WntGradient wnt_gradient;
        MeinekeCryptCell wnt_cell(TRANSIT, // type
                                  APC_ONE_HIT,//Mutation State
                                  1,    // generation
                                  new WntCellCycleModel(wnt_stimulus, wnt_gradient));
                                  
        MeinekeCryptCell wnt_cell2(TRANSIT, // type
                                  BETA_CATENIN_ONE_HIT,//Mutation State
                                  1,    // generation
                                  new WntCellCycleModel(wnt_stimulus, wnt_gradient));                          
                                  
        MeinekeCryptCell wnt_cell3(TRANSIT, // type
                                  APC_TWO_HIT,//Mutation State
                                  1,    // generation
                                  new WntCellCycleModel(wnt_stimulus, wnt_gradient)); 
                                  
        MeinekeCryptCell wnt_cell4(TRANSIT, // type
                                  LABELLED,//Mutation State
                                  1,    // generation
                                  new WntCellCycleModel(wnt_stimulus, wnt_gradient));                               

        std::vector<double> wnt;
        wnt.push_back(wnt_stimulus);

        TS_ASSERT(wnt_cell.ReadyToDivide(wnt)==false);

        TS_ASSERT(wnt_cell2.ReadyToDivide(wnt)==false);
        
        TS_ASSERT(wnt_cell3.ReadyToDivide(wnt)==false);
        
        TS_ASSERT(wnt_cell4.ReadyToDivide(wnt)==false);
        
        SimulationTime::Destroy();
    }
};



#endif /*TESTMEINEKECRYPTCELL_HPP_*/
