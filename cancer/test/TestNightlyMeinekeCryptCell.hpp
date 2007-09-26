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
#include "WntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "SimulationTime.hpp"
#include <iostream>

class TestNightlyMeinekeCryptCell: public CxxTest::TestSuite
{
public:

    void TestTysonNovakImmortalStemCell()
    {
        CancerParameters::Instance()->Reset();

        double end_time=100.0; // A good load of divisions to make sure nothing mucks up..
        // one division = 1.26 hours.
        int time_steps=1000;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
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
        unsigned divisions=0;
        while (p_simulation_time->GetDimensionalisedTime()< end_time)
        {
            // produce the offspring of the cells
            
            p_simulation_time->IncrementTimeOneStep();
            times[i]=p_simulation_time->GetDimensionalisedTime();
            cell_iterator = cells.begin();
            unsigned j=0;
            while (cell_iterator < cells.end())
            {
                //std::cout << "Cell Cycle Model called" << std::endl;
                if (cell_iterator->ReadyToDivide())
                {
                    MeinekeCryptCell new_cell = cell_iterator->Divide();
                    divisions++;
                    //std::cout << "Division of stem cell at time = " << times[i] << std::endl;
                }
                cell_iterator++;
                j++;
            }
            i++;
        }
        TS_ASSERT_DELTA(divisions,(unsigned)(end_time/1.26),1);
        SimulationTime::Destroy();
    }
    
    void Test0DBucketWithTysonNovak()
    {
        CancerParameters::Instance()->Reset();

        double end_time=7.0; // not very long because cell cycle time is only 1.2
        // (75 mins) because Tyson Novaks is for yeast
        int time_steps=100;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
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
            times[i]=p_simulation_time->GetDimensionalisedTime();
            cell_iterator = cells.begin();
            unsigned j=0;
            while (cell_iterator < cells.end())
            {
            
                //std::cout << "Cell " << j << " is generation " << cell_iterator->GetGeneration() << std::endl;
                if (cell_iterator->ReadyToDivide())
                {
                    newly_born.push_back(cell_iterator->Divide());
                    //std::cout << "Division of cell "<<j<<" at time = " << times[i] << std::endl;
                }
                cell_iterator++;
                j++;
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
            
            i++;
        }
        
        TS_ASSERT_EQUALS(stem_cells[time_steps-1], 1u);
        TS_ASSERT_EQUALS(transit_cells[time_steps-1], 7u);
        TS_ASSERT_EQUALS(differentiated_cells[time_steps-1], 24u);
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
        CancerParameters::Instance()->Reset();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        CancerParameters *p_parameters = CancerParameters::Instance();
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        unsigned num_steps=200;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);
        
        double wnt_stimulus = 1.0;
        SingletonWntGradient::Instance()->SetConstantWntValueForTesting(wnt_stimulus);
        
        MeinekeCryptCell wnt_cell(TRANSIT, // type
                                  APC_ONE_HIT,//Mutation State
                                  1,    // generation
                                  new WntCellCycleModel(wnt_stimulus));
                                  
        wnt_cell.InitialiseCellCycleModel();
        
        for (unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();

            if (time>=4.804+SG2MDuration)
            {
                TS_ASSERT(wnt_cell.ReadyToDivide()==true);
            }
            else
            {
                TS_ASSERT(wnt_cell.ReadyToDivide()==false);
            }
            //std::cout << "Time = " << time << " ready = " << wnt_cell.ReadyToDivide(wnt) << "\n" << std::endl;
        }

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT(wnt_cell.ReadyToDivide()==true);
        TS_ASSERT(wnt_cell.GetGeneration()==1);
        
        MeinekeCryptCell wnt_cell2 = wnt_cell.Divide();
        
        TS_ASSERT(wnt_cell.GetGeneration()==2);
        TS_ASSERT(wnt_cell2.GetGeneration()==2);
        
        //std::cout << "time now = " << p_simulation_time->GetDimensionalisedTime() << "\n" <<std::endl;
        
        double time_of_birth = wnt_cell.GetBirthTime();
        double time_of_birth2 = wnt_cell2.GetBirthTime();
        
        TS_ASSERT_DELTA(time_of_birth, time_of_birth2, 1e-9);
        
        //std::cout << "time of cell divisions = " << timeOfBirth << "\tand\t" << timeOfBirth2 << "\n" << std::endl;
        
        for (unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();

            bool result1 = wnt_cell.ReadyToDivide();
            bool result2 = wnt_cell2.ReadyToDivide();
            
            if ( time >= 4.804+SG2MDuration+time_of_birth )
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
        SingletonWntGradient::Destroy();
    }
    
    /*
     * We are checking that the MeinekeCryptCells work with the Wnt cell cycle models here
     * That division of wnt cells and stuff works OK.
     * 
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModelAndMutationBetaCat() throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        CancerParameters *p_parameters = CancerParameters::Instance();
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        unsigned num_steps=200;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);
        
        double wnt_stimulus = 0.0;
        SingletonWntGradient::Instance()->SetConstantWntValueForTesting(wnt_stimulus);
        
        MeinekeCryptCell wnt_cell(TRANSIT, // type
                                  BETA_CATENIN_ONE_HIT,//Mutation State
                                  1,    // generation
                                  new WntCellCycleModel(wnt_stimulus));
                                  
        wnt_cell.InitialiseCellCycleModel();
                                  
        for (unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();

            if (time>=7.82+SG2MDuration)
            {
                TS_ASSERT(wnt_cell.ReadyToDivide()==true);
            }
            else
            {
                TS_ASSERT(wnt_cell.ReadyToDivide()==false);
            }
            //std::cout << "Time = " << time << " ready = " << wnt_cell.ReadyToDivide(wnt) << "\n" << std::endl;
        }
        
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT(wnt_cell.ReadyToDivide()==true);
        TS_ASSERT(wnt_cell.GetGeneration()==1);
        
        MeinekeCryptCell wnt_cell2 = wnt_cell.Divide();
        
        TS_ASSERT(wnt_cell.GetGeneration()==2);
        TS_ASSERT(wnt_cell2.GetGeneration()==2);
        
        //std::cout << "time now = " << p_simulation_time->GetDimensionalisedTime() << "\n" <<std::endl;
        
        double time_of_birth = wnt_cell.GetBirthTime();
        double time_of_birth2 = wnt_cell2.GetBirthTime();
        
        TS_ASSERT_DELTA(time_of_birth, time_of_birth2, 1e-9);
        
        //std::cout << "time of cell divisions = " << timeOfBirth << "\tand\t" << timeOfBirth2 << "\n" << std::endl;
        
        for (unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();

            bool result1=wnt_cell.ReadyToDivide();
            bool result2=wnt_cell2.ReadyToDivide();
            if ( time >= 7.82+SG2MDuration+time_of_birth )
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
        SingletonWntGradient::Destroy();
    }
    
    /*
     * We are checking that the MeinekeCryptCells work with the Wnt cell cycle models here
     * That division of wnt cells and stuff works OK.
     * 
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModelAndMutationAPC2() throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        CancerParameters *p_parameters = CancerParameters::Instance();
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        unsigned num_steps=200;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);
        
        double wnt_stimulus = 0.0;
        SingletonWntGradient::Instance()->SetConstantWntValueForTesting(wnt_stimulus);
        
        MeinekeCryptCell wnt_cell(TRANSIT, // type
                                  APC_TWO_HIT,//Mutation State
                                  1,    // generation
                                  new WntCellCycleModel(wnt_stimulus));
        
        wnt_cell.InitialiseCellCycleModel();
                                          
        CryptCellMutationState this_state = wnt_cell.GetMutationState();
        
        TS_ASSERT(this_state==APC_TWO_HIT);
        
        this_state = APC_ONE_HIT;
        
        wnt_cell.SetMutationState(this_state);
        
        this_state = wnt_cell.GetMutationState();
        
        TS_ASSERT(this_state==APC_ONE_HIT);
        
        wnt_cell.SetMutationState(APC_TWO_HIT);
        
        for (unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();

            if (time>=3.9435+SG2MDuration)
            {
                TS_ASSERT(wnt_cell.ReadyToDivide()==true);
            }
            else
            {
                TS_ASSERT(wnt_cell.ReadyToDivide()==false);
            }
            //std::cout << "Time = " << time << " ready = " << wnt_cell.ReadyToDivide(wnt) << "\n" << std::endl;
        }
        
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT(wnt_cell.ReadyToDivide()==true);
        TS_ASSERT(wnt_cell.GetGeneration()==1);
        
        MeinekeCryptCell wnt_cell2 = wnt_cell.Divide();
        
        TS_ASSERT(wnt_cell.GetGeneration()==2);
        TS_ASSERT(wnt_cell2.GetGeneration()==2);
        
        //std::cout << "time now = " << p_simulation_time->GetDimensionalisedTime() << "\n" <<std::endl;
        
        double time_of_birth = wnt_cell.GetBirthTime();
        double time_of_birth2 = wnt_cell2.GetBirthTime();
        
        TS_ASSERT_DELTA(time_of_birth, time_of_birth2, 1e-9);
        
        //std::cout << "time of cell divisions = " << timeOfBirth << "\tand\t" << timeOfBirth2 << "\n" << std::endl;
        
        for (unsigned i=0 ; i<num_steps/2 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();

            bool result1=wnt_cell.ReadyToDivide();
            bool result2=wnt_cell2.ReadyToDivide();
            if ( time >= 3.9435+SG2MDuration+time_of_birth )
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
        SingletonWntGradient::Destroy();
    }
    
    
};



#endif /*TESTMEINEKECRYPTCELL_HPP_*/
