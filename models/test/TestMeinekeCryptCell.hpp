#ifndef TESTMEINEKECRYPTCELL_HPP_
#define TESTMEINEKECRYPTCELL_HPP_

#include <cxxtest/TestSuite.h>

#include "MeinekeCryptCellTypes.hpp"
#include "CryptCellMutationStates.hpp"
#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "SimulationTime.hpp"

class TestMeinekeCryptCell: public CxxTest::TestSuite
{
public:

    void TestCellsAgeingCorrectly() throw(Exception)
    {
        // These lines are added to cover the exception case that a cell is 
        // created without simulation time being set up...
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);
    	FixedCellCycleModel fixed_model;
        SimulationTime::Destroy();
        
        TS_ASSERT_THROWS_ANYTHING(MeinekeCryptCell bad_cell(STEM, // type
        						   HEALTHY,//Mutation State
                                   0,    // generation
                                   &fixed_model));
                                   
        // Proper test again
    	
        p_simulation_time = SimulationTime::Instance();
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
        
        SimulationTime::Destroy();
    }
    
    void TestCellDivision()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
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
        SimulationTime* p_simulation_time = SimulationTime::Instance();
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
    
    void Test0DBucket()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(60.0, 60);
        
        MeinekeCryptCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,  // generation
                                   new FixedCellCycleModel());
                                   
                                   
        double end_time=60.0;
        int time_steps=60;
        
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
    
    void TestWithCellCycleModel() throw(Exception)
    {
        // Simulation time is 6000 because we want to test that differentiated cells never divide.
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(6000.0, 1000);
        CancerParameters *p_params = CancerParameters::Instance();
        
        // this test needs particular cell cycle times
        TS_ASSERT_EQUALS(p_params->GetStemCellCycleTime(), 24.0);
        TS_ASSERT_EQUALS(p_params->GetTransitCellCycleTime(), 12.0);
        
        RandomNumberGenerator rand_gen;
        
        p_simulation_time->IncrementTimeOneStep();
        
        
        //  Creating different types of Meineke crypt cells with different cell cycle models at SImulationTime = 6 hours.
        MeinekeCryptCell stem_cell(STEM, // type
								   HEALTHY,//Mutation State
                                   0,    // generation
                                   new FixedCellCycleModel());
                                   
        MeinekeCryptCell stochastic_stem_cell(STEM, // type
                                              HEALTHY,//Mutation State
                                   			  0,    // generation
                                              new StochasticCellCycleModel(&rand_gen));
        MeinekeCryptCell differentiated_cell(DIFFERENTIATED, // type
                                             HEALTHY,//Mutation State
                                   			 6,    // generation
                                             new FixedCellCycleModel());
        MeinekeCryptCell stochastic_differentiated_cell(DIFFERENTIATED, // type
                                                        HEALTHY,//Mutation State
                                   			 			6,    // generation
                                                        new StochasticCellCycleModel(&rand_gen));
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
        
    }
    
 
    
    void TestStochasticCycleModel() throw(Exception)
    {
    
        // Go up in steps of 0.01 to test stochasticity in cell cycle models
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(54.0, 5400);
        CancerParameters *p_params = CancerParameters::Instance();
        
        // this test needs particular cell cycle times
        TS_ASSERT_EQUALS(p_params->GetStemCellCycleTime(), 24.0);
        TS_ASSERT_EQUALS(p_params->GetTransitCellCycleTime(), 12.0);
        
        RandomNumberGenerator rand_gen;
        
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
        
        StochasticCellCycleModel *cell_cycle_model = new StochasticCellCycleModel(&rand_gen);
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
    }
    
    
    
    void Test0DBucketStochastic()
    {
    
        CancerParameters *p_params = CancerParameters::Instance();
        
        // this test needs particular cell cycle times
        TS_ASSERT_EQUALS(p_params->GetStemCellCycleTime(), 24.0);
        TS_ASSERT_EQUALS(p_params->GetTransitCellCycleTime(), 12.0);
        
        RandomNumberGenerator rand_gen;
        
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
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(70.0, 70);
            
            MeinekeCryptCell stem_cell(STEM, // type
                                       HEALTHY,//Mutation State
                                   	   0,  // generation
                                       new StochasticCellCycleModel(&rand_gen));
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
        
        
    }
    
    /* We are setting up a 0d bucket with some initial cell population
     * This is deterministic so we can test it
     */
    void TestInitialise0DBucket()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
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
        double end_time=60.0;
        int time_steps=60;
        //double time_step= end_time/(double) time_steps;
        
        std::vector<MeinekeCryptCell> newly_born;
        std::vector<unsigned> stem_cells(time_steps);
        std::vector<unsigned> transit_cells(time_steps);
        std::vector<unsigned> differentiated_cells(time_steps);
        std::vector<double> times(time_steps);
        
        std::vector<MeinekeCryptCell>::iterator cell_iterator;
        
        unsigned i=0;
        while (p_simulation_time->GetDimensionalisedTime()< end_time)
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
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        unsigned num_steps=100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);
        
        double wnt_stimulus = 1.0;
        MeinekeCryptCell wnt_cell(TRANSIT, // type
                                   HEALTHY,//Mutation State
                                   1,    // generation
                                   new WntCellCycleModel(wnt_stimulus));
        
        for(unsigned i=0 ; i<num_steps/2 ; i++)
        {
        	p_simulation_time->IncrementTimeOneStep();
        	double time = p_simulation_time->GetDimensionalisedTime();
        	std::vector<double> wnt;
        	wnt.push_back(wnt_stimulus);
        	if(time>=5.971+SG2MDuration)
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
        
        for(unsigned i=0 ; i<num_steps/2 ; i++)
        {
        	p_simulation_time->IncrementTimeOneStep();
        	double time = p_simulation_time->GetDimensionalisedTime();
        	std::vector<double> wnt;
        	wnt.push_back(wnt_stimulus);
        	bool result1=wnt_cell.ReadyToDivide(wnt);
        	bool result2=wnt_cell2.ReadyToDivide(wnt);
        	//std::cout << "Time = " << time << ", ready1 = " << result1 << ", ready2 = " << result2<< "\n" << std::endl;
        	if(time>=5.971+SG2MDuration+timeOfBirth)
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
     * We are checking that the MeinekeCryptCells work with the T&N cell cycle models here
     * That division of wnt cells and stuff works OK.
     * 
     * It checks that the cell division thing works nicely too.
     */
    void TestWithTysonNovakCellCycleModel() throw(Exception)
    {
        
        double standard_tyson_duration = 75.19/60.0;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        //CancerParameters *p_parameters = CancerParameters::Instance();
        //double SG2MDuration = p_parameters->GetSG2MDuration();
        
        unsigned num_steps=100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(200.0/60.0, num_steps+1);
        
        MeinekeCryptCell tn_cell(TRANSIT, // type
                                   HEALTHY,//Mutation State
                                   1,    // generation
                                   new TysonNovakCellCycleModel());
        
        for(unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            if(time>=standard_tyson_duration)
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
        
        for(unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            bool result1=tn_cell.ReadyToDivide();
            bool result2=tn_cell2.ReadyToDivide();
            //std::cout << "Time = " << time << ", ready1 = " << result1 << ", ready2 = " << result2<< "\n" << std::endl;
            if(time>=standard_tyson_duration+timeOfBirth)
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
    
    
    
    void Test0DBucketWithTysonNovak()
    {
        double end_time=5.0; // not very long because cell cycle time is only 1.2 
                             // (75 mins) because Tyson Novaks is for yeast
        int time_steps=60;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, time_steps);
        
        MeinekeCryptCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,  // generation
                                   new TysonNovakCellCycleModel());
                                   
        
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
        
        
//        for(int i=0; i<60; i++)
//        {
//            std::cout << i << ": " << stem_cells[i] << " " << transit_cells[i] 
//                      << " " << differentiated_cells[i] << "\n";
//        }
        
        TS_ASSERT_EQUALS(stem_cells[59], 1u);
        TS_ASSERT_EQUALS(transit_cells[59], 7u);
        TS_ASSERT_EQUALS(differentiated_cells[59], 8u);
        SimulationTime::Destroy();
    }
    
    /*
     * We are checking that the MeinekeCryptCells work with the Wnt cell cycle models here
     * That division of wnt cells and stuff works OK.
     * 
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModelAndMutationAPCONEHIT() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        CancerParameters *p_parameters = CancerParameters::Instance();
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        unsigned num_steps=200;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);
        
        double wnt_stimulus = 1.0;
        MeinekeCryptCell wnt_cell(TRANSIT, // type
                                   APC_ONE_HIT,//Mutation State
                                   1,    // generation
                                   new WntCellCycleModel(wnt_stimulus,1));
       
        for(unsigned i=0 ; i<num_steps/2 ; i++)
        {
        	p_simulation_time->IncrementTimeOneStep();
        	double time = p_simulation_time->GetDimensionalisedTime();
        	std::vector<double> wnt;
        	wnt.push_back(wnt_stimulus);
        	//std::cout << "time = " << time << std::endl;
        	if(time>=4.804+SG2MDuration)
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
        
        for(unsigned i=0 ; i<num_steps/2 ; i++)
        {
        	p_simulation_time->IncrementTimeOneStep();
        	double time = p_simulation_time->GetDimensionalisedTime();
        	std::vector<double> wnt;
        	wnt.push_back(wnt_stimulus);
        	bool result1=wnt_cell.ReadyToDivide(wnt);
        	bool result2=wnt_cell2.ReadyToDivide(wnt);
        	//std::cout << "Time = " << time << ", ready1 = " << result1 << ", ready2 = " << result2<< "\n" << std::endl;
        	if(time>=4.804+SG2MDuration+timeOfBirth)
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
     * We are checking that the MeinekeCryptCells work with the Wnt cell cycle models here
     * That division of wnt cells and stuff works OK.
     * 
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModelAndMutationBetaCat() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        CancerParameters *p_parameters = CancerParameters::Instance();
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        unsigned num_steps=200;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);
        
        double wnt_stimulus = 0.0;
        MeinekeCryptCell wnt_cell(TRANSIT, // type
                                   BETA_CATENIN_ONE_HIT,//Mutation State
                                   1,    // generation
                                   new WntCellCycleModel(wnt_stimulus,2));
       
        for(unsigned i=0 ; i<num_steps/2 ; i++)
        {
        	p_simulation_time->IncrementTimeOneStep();
        	double time = p_simulation_time->GetDimensionalisedTime();
        	std::vector<double> wnt;
        	wnt.push_back(wnt_stimulus);
        	//std::cout << "time = " << time << std::endl;
        	if(time>=7.82+SG2MDuration)
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
        
        for(unsigned i=0 ; i<num_steps/2 ; i++)
        {
        	p_simulation_time->IncrementTimeOneStep();
        	double time = p_simulation_time->GetDimensionalisedTime();
        	std::vector<double> wnt;
        	wnt.push_back(wnt_stimulus);
        	bool result1=wnt_cell.ReadyToDivide(wnt);
        	bool result2=wnt_cell2.ReadyToDivide(wnt);
        	//std::cout << "Time = " << time << ", ready1 = " << result1 << ", ready2 = " << result2<< "\n" << std::endl;
        	if(time>=7.82+SG2MDuration+timeOfBirth)
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
     * We are checking that the MeinekeCryptCells work with the Wnt cell cycle models here
     * That division of wnt cells and stuff works OK.
     * 
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModelAndMutationAPC2() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        CancerParameters *p_parameters = CancerParameters::Instance();
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        unsigned num_steps=200;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);
        
        double wnt_stimulus = 0.0;
        MeinekeCryptCell wnt_cell(TRANSIT, // type
                                   APC_ONE_HIT,//Mutation State
                                   1,    // generation
                                   new WntCellCycleModel(wnt_stimulus,3));
       
        CryptCellMutationState this_state = wnt_cell.GetMutationState();
        
        TS_ASSERT(this_state==APC_ONE_HIT);
        
        this_state = APC_TWO_HIT;
        
        wnt_cell.SetMutationState(this_state);
        
        this_state = wnt_cell.GetMutationState();
        
        TS_ASSERT(this_state==APC_TWO_HIT);
       
        for(unsigned i=0 ; i<num_steps/2 ; i++)
        {
        	p_simulation_time->IncrementTimeOneStep();
        	double time = p_simulation_time->GetDimensionalisedTime();
        	std::vector<double> wnt;
        	wnt.push_back(wnt_stimulus);
        	//std::cout << "time = " << time << std::endl;
        	if(time>=3.9435+SG2MDuration)
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
        
        for(unsigned i=0 ; i<num_steps/2 ; i++)
        {
        	p_simulation_time->IncrementTimeOneStep();
        	double time = p_simulation_time->GetDimensionalisedTime();
        	std::vector<double> wnt;
        	wnt.push_back(wnt_stimulus);
        	bool result1=wnt_cell.ReadyToDivide(wnt);
        	bool result2=wnt_cell2.ReadyToDivide(wnt);
        	//std::cout << "Time = " << time << ", ready1 = " << result1 << ", ready2 = " << result2<< "\n" << std::endl;
        	if(time>=3.9435+SG2MDuration+timeOfBirth)
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
    
};



#endif /*TESTMEINEKECRYPTCELL_HPP_*/
