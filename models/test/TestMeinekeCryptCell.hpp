#ifndef TESTMEINEKECRYPTCELL_HPP_
#define TESTMEINEKECRYPTCELL_HPP_

#include <cxxtest/TestSuite.h>

#include "MeinekeCryptCellTypes.hpp"
#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"


class TestMeinekeCryptCell: public CxxTest::TestSuite
{
public:

    void TestMeinekeCryptCellClass()
    {
        MeinekeCryptCell stem_cell(STEM, // type
                                   0.1,  // birth time (hours)
                                   0,    // generation
                                   new FixedCellCycleModel());
                                   
        TS_ASSERT_EQUALS(stem_cell.GetAge(1.0), 0.9);
        
        //for coverage
        stem_cell.SetNodeIndex(3);
        TS_ASSERT_EQUALS((int)stem_cell.GetNodeIndex, (int)3);
    }
    
    void TestCellDivision()
    {
        CancerParameters *p_params = CancerParameters::Instance();
        const double birth_time = 0.1;
        
        MeinekeCryptCell stem_cell(STEM, // type
                                   birth_time,  // birth time (hours)
                                   0,    // generation
                                   new FixedCellCycleModel());
                                   
        TS_ASSERT(!stem_cell.ReadyToDivide(p_params->GetStemCellCycleTime()));
        
        TS_ASSERT(stem_cell.ReadyToDivide(p_params->GetStemCellCycleTime()+birth_time));
        
        // create transit progeny of stem
        MeinekeCryptCell daughter_cell = stem_cell.Divide(p_params->GetStemCellCycleTime()+birth_time);
        
        TS_ASSERT(!stem_cell.ReadyToDivide(p_params->GetStemCellCycleTime()+birth_time));
        
        TS_ASSERT(daughter_cell.GetGeneration() == 1);
        TS_ASSERT(daughter_cell.GetCellType() == TRANSIT);
        TS_ASSERT(daughter_cell.GetAge(p_params->GetStemCellCycleTime()+birth_time) == 0);
        
        TS_ASSERT(!daughter_cell.ReadyToDivide(p_params->GetStemCellCycleTime()+
                                               p_params->GetTransitCellCycleTime()));
        TS_ASSERT(daughter_cell.ReadyToDivide(p_params->GetStemCellCycleTime()+
                                              p_params->GetTransitCellCycleTime()+birth_time));
                                              
        // create transit progeny of transit
        MeinekeCryptCell grandaughter_cell = daughter_cell.Divide(p_params->GetStemCellCycleTime()+
                                                                  p_params->GetTransitCellCycleTime()+birth_time);
                                                                  
        TS_ASSERT(!stem_cell.ReadyToDivide(2*p_params->GetStemCellCycleTime()));
        // stem cell ready to divide again
        TS_ASSERT(stem_cell.ReadyToDivide(2*p_params->GetStemCellCycleTime()+birth_time));
        
        // both grandaughter and daughter cells should be ready to
        // divide at p_params->GetStemCellCycleTime()+2*p_params->GetTransitCellCycleTime()+0.1
        
        TS_ASSERT(!grandaughter_cell.ReadyToDivide(p_params->GetStemCellCycleTime()+
                                                   2*p_params->GetTransitCellCycleTime()));
        TS_ASSERT(!daughter_cell.ReadyToDivide(p_params->GetStemCellCycleTime()+
                                               2*p_params->GetTransitCellCycleTime()));
                                               
                                               
        TS_ASSERT(grandaughter_cell.ReadyToDivide(p_params->GetStemCellCycleTime()+
                                                  2*p_params->GetTransitCellCycleTime()+birth_time));
        TS_ASSERT(daughter_cell.ReadyToDivide(p_params->GetStemCellCycleTime()+
                                              2*p_params->GetTransitCellCycleTime()+birth_time));
    }
    
    void TestCellDivisionStops()
    {
        CancerParameters *p_params = CancerParameters::Instance();
        const double birth_time = 0.1;
        
        MeinekeCryptCell stem_cell(STEM, // type
                                   birth_time,  // birth time (hours)
                                   0,    // generation
                                   new FixedCellCycleModel());
                                   
        double time = p_params->GetStemCellCycleTime() + birth_time;
        // create transit progeny of stem
        TS_ASSERT(stem_cell.ReadyToDivide(time));
        MeinekeCryptCell daughter_cell = stem_cell.Divide(time);
        
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
            time += p_params->GetTransitCellCycleTime();
            while (cell_iterator < cells.end())
            {
                if (cell_iterator->ReadyToDivide(time))
                {
                    newly_born.push_back(cell_iterator->Divide(time));
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
    }
    
    void Test0DBucket()
    {
        MeinekeCryptCell stem_cell(STEM, // type
                                   0,  // birth time (hours)
                                   0,  // generation
                                   new FixedCellCycleModel());
                                   
        double time=0.0;
        double end_time=60.0;
        int time_steps=60;
        double time_step= end_time/(double) time_steps;
        
        std::vector<MeinekeCryptCell> cells;
        std::vector<MeinekeCryptCell> newly_born;
        std::vector<unsigned> stem_cells(time_steps);
        std::vector<unsigned> transit_cells(time_steps);
        std::vector<unsigned> differentiated_cells(time_steps);
        std::vector<double> times(time_steps);
        
        cells.push_back(stem_cell);
        std::vector<MeinekeCryptCell>::iterator cell_iterator;
        
        unsigned i=0;
        while (time< end_time)
        {
            // produce the offspring of the cells
            time+=time_step;
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                if (cell_iterator->ReadyToDivide(time))
                {
                    newly_born.push_back(cell_iterator->Divide(time));
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
            times[i]=time;
            i++;
        }
        TS_ASSERT_EQUALS(stem_cells[59], 1u);
        TS_ASSERT_EQUALS(transit_cells[59], 2u);
        TS_ASSERT_EQUALS(differentiated_cells[59], 8u);
    }
    
    void TestWithCellCycleModel() throw(Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        const double birth_time = 0.1;
        // Test Stem cell
        MeinekeCryptCell stem_cell(STEM, // type
                                   birth_time,  // birth time (hours)
                                   0,    // generation
                                   new FixedCellCycleModel());
        TS_ASSERT(!stem_cell.ReadyToDivide(p_params->GetStemCellCycleTime()));
        TS_ASSERT(stem_cell.ReadyToDivide(p_params->GetStemCellCycleTime()+birth_time));
        MeinekeCryptCell daughter_cell1 = stem_cell.Divide(p_params->GetStemCellCycleTime()+birth_time);
        TS_ASSERT(typeid(daughter_cell1.GetCellCycleModel()) == typeid(stem_cell.GetCellCycleModel()));
        
        MeinekeCryptCell stochastic_stem_cell(STEM, // type
                                              birth_time,  // birth time (hours)
                                              0,    // generation
                                              new StochasticCellCycleModel());
        TS_ASSERT(!stochastic_stem_cell.ReadyToDivide(p_params->GetStemCellCycleTime()));
        TS_ASSERT(stochastic_stem_cell.ReadyToDivide(p_params->GetStemCellCycleTime()+birth_time));
        
        // Test Fully-differentiated cell
        MeinekeCryptCell differentiated_cell(DIFFERENTIATED, // type
                                             birth_time,  // birth time (hours)
                                             6,    // generation
                                             new FixedCellCycleModel());
        TS_ASSERT(!differentiated_cell.ReadyToDivide(1e5));
        TS_ASSERT(!differentiated_cell.ReadyToDivide(1e50));
        TS_ASSERT(!differentiated_cell.ReadyToDivide(1e150));
        
        MeinekeCryptCell stochastic_differentiated_cell(DIFFERENTIATED, // type
                                                        birth_time,  // birth time (hours)
                                                        6,    // generation
                                                        new StochasticCellCycleModel());
        TS_ASSERT(!stochastic_differentiated_cell.ReadyToDivide(1e5));
        TS_ASSERT(!stochastic_differentiated_cell.ReadyToDivide(1e50));
        TS_ASSERT(!stochastic_differentiated_cell.ReadyToDivide(1e150));
        
        // Test transit cell
        MeinekeCryptCell transit_cell(TRANSIT, // type
                                      birth_time,  // birth time (hours)
                                      2,    // generation
                                      new FixedCellCycleModel());
        TS_ASSERT(!transit_cell.ReadyToDivide(p_params->GetTransitCellCycleTime()+birth_time-0.01));
        TS_ASSERT(transit_cell.ReadyToDivide(p_params->GetTransitCellCycleTime()+birth_time));
        
        StochasticCellCycleModel *cell_cycle_model = new StochasticCellCycleModel();
        transit_cell.SetCellCycleModel(cell_cycle_model);
        TS_ASSERT_EQUALS(transit_cell.GetCellCycleModel(), cell_cycle_model);
        int ready_count=0;
        for (int i=0; i<100; i++)
        {
            if (transit_cell.ReadyToDivide(p_params->GetTransitCellCycleTime()+birth_time-0.01))
            {
                ready_count++;
            }
        }
        TS_ASSERT(ready_count>0);
        
        double divide_time = p_params->GetTransitCellCycleTime();
        while (!transit_cell.ReadyToDivide(divide_time))
        {
            divide_time += 0.1;
        }
        MeinekeCryptCell daughter_cell2 = transit_cell.Divide(divide_time);
        TS_ASSERT(typeid(daughter_cell2.GetCellCycleModel()) == typeid(transit_cell.GetCellCycleModel()));
    }
    
    void Test0DBucketStochastic()
    {
        const double end_time = 70.0;
        const int time_steps = 70;
        const int number_of_simulations = 1000;
        const double time_step= end_time/(double) time_steps;
        
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
            MeinekeCryptCell stem_cell(STEM, // type
                                       0,  // birth time (hours)
                                       0,  // generation
                                       new StochasticCellCycleModel());
            cells.push_back(stem_cell);
            // produce the offspring of the cells
            std::vector<MeinekeCryptCell>::iterator cell_iterator = cells.begin();
            
            double time = 0.0;
            while (time< end_time)
            {
                time+=time_step;
                cell_iterator = cells.begin();
                while (cell_iterator < cells.end())
                {
                    if (cell_iterator->ReadyToDivide(time))
                    {
                        newly_born.push_back(cell_iterator->Divide(time));
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
        }
        stem_cell_mean=stem_cell_mean/(double) number_of_simulations;
        transit_cell_mean=transit_cell_mean/(double) number_of_simulations;
        differentiated_cell_mean=differentiated_cell_mean/(double) number_of_simulations;
        
        TS_ASSERT_DELTA(stem_cell_mean, 1.0, 1e-12);
        TS_ASSERT_DELTA(transit_cell_mean, 2.0, 1.0);
        TS_ASSERT_DELTA(differentiated_cell_mean, 8.0, 1.0);
        
        CancerParameters *p_params = CancerParameters::Instance();
        TS_ASSERT_DELTA(p_params->GetTransitCellCycleTime(), 12.0, 1e-12);
    }
    
    /* We are setting up a 0d bucket with some initial cell population
     * This is deterministic so we can test it
     */
    void TestInitialise0DBucket()
    {
        std::vector<MeinekeCryptCell> cells;
        
        MeinekeCryptCell stem_cell(STEM, // type
                                   0,  // birth time (hours)
                                   0,  // generation
                                   new FixedCellCycleModel());
                                   
        cells.push_back(stem_cell);
        
        MeinekeCryptCell transit_cell_1(TRANSIT, // type
                                        0,  // birth time (hours)
                                        1,  // generation
                                        new FixedCellCycleModel());
                                        
        cells.push_back(transit_cell_1);
        
        MeinekeCryptCell transit_cell_2(TRANSIT, // type
                                        0,  // birth time (hours)
                                        2,  // generation
                                        new FixedCellCycleModel());
                                        
        cells.push_back(transit_cell_2);
        
        MeinekeCryptCell transit_cell_3(TRANSIT, // type
                                        0,  // birth time (hours)
                                        3,  // generation
                                        new FixedCellCycleModel());
                                        
        cells.push_back(transit_cell_3);
        
        MeinekeCryptCell differentiated_cell(DIFFERENTIATED, // type
                                             0,  // birth time (hours)
                                             4,  // generation
                                             new FixedCellCycleModel());
                                             
        cells.push_back(differentiated_cell);
        
        double time=0.0;
        double end_time=60.0;
        int time_steps=60;
        double time_step= end_time/(double) time_steps;
        
        std::vector<MeinekeCryptCell> newly_born;
        std::vector<unsigned> stem_cells(time_steps);
        std::vector<unsigned> transit_cells(time_steps);
        std::vector<unsigned> differentiated_cells(time_steps);
        std::vector<double> times(time_steps);
        
        std::vector<MeinekeCryptCell>::iterator cell_iterator;
        
        unsigned i=0;
        while (time< end_time)
        {
            time+=time_step;
            // produce the offspring of the cells
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                if (cell_iterator->ReadyToDivide(time))
                {
                    newly_born.push_back(cell_iterator->Divide(time));
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
            
            times[i]=time;
            i++;
        }
        
        TS_ASSERT_EQUALS(stem_cells[59], 1u);
        TS_ASSERT_EQUALS(transit_cells[59], 2u);
        TS_ASSERT_EQUALS(differentiated_cells[59], 23u);
    }
};



#endif /*TESTMEINEKECRYPTCELL_HPP_*/
