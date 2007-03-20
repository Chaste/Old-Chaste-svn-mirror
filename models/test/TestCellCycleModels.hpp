#ifndef TESTCELLCYCLEMODELS_HPP_
#define TESTCELLCYCLEMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "OutputFileHandler.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"

class TestCellCycleModels : public CxxTest::TestSuite
{
public:
    void TestFixedCellCycleModel(void) throw(Exception)
    {
        // Make cell cycle models protest if simulation time is not set up
        TS_ASSERT_THROWS_ANYTHING(FixedCellCycleModel model1);
        
        CancerParameters *p_params = CancerParameters::Instance();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        TS_ASSERT_THROWS_ANYTHING(FixedCellCycleModel model2);
        
        unsigned num_steps = 100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0*p_params->GetStemCellCycleTime(), num_steps);
        
        TS_ASSERT_THROWS_NOTHING(FixedCellCycleModel model3);
        
        FixedCellCycleModel our_fixed_stem_cell_cycle_model;
        our_fixed_stem_cell_cycle_model.SetCellType(STEM);
        TS_ASSERT_EQUALS(our_fixed_stem_cell_cycle_model.GetCellType(),STEM);
        
        FixedCellCycleModel our_fixed_transit_cell_cycle_model;
        our_fixed_transit_cell_cycle_model.SetCellType(TRANSIT);
        
        FixedCellCycleModel our_fixed_diff_cell_cycle_model;
        our_fixed_diff_cell_cycle_model.SetCellType(DIFFERENTIATED);
        
        for (unsigned i = 0 ; i< num_steps ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            //std::cout << "Time = " << time << " cell age = " << our_fixed_stem_cell_cycle_model.GetAge() << "\n";
            // Test STEM cells
            if (time<p_params->GetStemCellCycleTime())
            {
                TS_ASSERT(!our_fixed_stem_cell_cycle_model.ReadyToDivide());
                //std::cout << "No stem division.\n";
            }
            else
            {
                TS_ASSERT(our_fixed_stem_cell_cycle_model.ReadyToDivide());
                //std::cout << "Stem Division.\n";
            }
            // Test a Transit Cell
            if (time<p_params->GetTransitCellCycleTime())
            {
                TS_ASSERT(!our_fixed_transit_cell_cycle_model.ReadyToDivide());
                //std::cout << "No transit division.\n";
            }
            else
            {
                TS_ASSERT(our_fixed_transit_cell_cycle_model.ReadyToDivide());
                //std::cout << "Transit Division.\n";
            }
            // Test a DIFFERENTIATED cell
            TS_ASSERT(!our_fixed_diff_cell_cycle_model.ReadyToDivide());
        }
        
        TS_ASSERT_DELTA(our_fixed_stem_cell_cycle_model.GetAge(), p_simulation_time->GetDimensionalisedTime(), 1e-9);
        TS_ASSERT_DELTA(our_fixed_transit_cell_cycle_model.GetAge(), p_simulation_time->GetDimensionalisedTime(), 1e-9);
        TS_ASSERT_DELTA(our_fixed_diff_cell_cycle_model.GetAge(), p_simulation_time->GetDimensionalisedTime(), 1e-9);
        SimulationTime::Destroy();
    }
    
    void TestStochasticCellCycleModel(void) throw(Exception)
    {
        TS_ASSERT_THROWS_ANYTHING(StochasticCellCycleModel cell_model1);
        
        RandomNumberGenerator::Instance();
        CancerParameters *p_params = CancerParameters::Instance();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        TS_ASSERT_THROWS_ANYTHING(StochasticCellCycleModel cell_model2);
        
        unsigned num_steps = 100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0*p_params->GetStemCellCycleTime(), num_steps);
        
        TS_ASSERT_THROWS_NOTHING(StochasticCellCycleModel cell_model3);
        
        StochasticCellCycleModel cell_model;
        
        for (unsigned i = 0 ; i< num_steps ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            //std::cout << "Time = " << time << " cell age = " << our_fixed_stem_cell_cycle_model.GetAge() << "\n";
            // Test STEM cells
            cell_model.SetCellType(STEM);
            if (time<p_params->GetStemCellCycleTime())
            {
                TS_ASSERT(!cell_model.ReadyToDivide());
                //std::cout << "No stem division.\n";
            }
            else
            {
                TS_ASSERT(cell_model.ReadyToDivide());
                //std::cout << "Stem Division.\n";
            }
            // Test Transit cells - new random number each time they are asked to divide...
            // shouldn't it be done so that they are given a random time when created?
            // Otherwise division time depends upon how often they are asked!
            cell_model.SetCellType(TRANSIT);
            const unsigned TESTS = 100;
            unsigned ready_count = 0;
            for (unsigned i=0; i<TESTS; i++)
            {
                if (cell_model.ReadyToDivide())
                {
                    ready_count++;
                }
            }
            //std::cout << time << " ready count = " << ready_count << "\n";
            if (time < 9.0)
            {
                TS_ASSERT(ready_count==0)
            }
            if (time > 15.0)
            {
                TS_ASSERT(ready_count==100)
            }
            if (time>11.75 && time < 12.25)
            {
                TS_ASSERT(ready_count==54)
            }
        }
        RandomNumberGenerator::Destroy();
        SimulationTime::Destroy();
    }
    
    
    void TestTysonNovakCellCycleModel(void) throw(Exception)
    {
        TS_ASSERT_THROWS_ANYTHING(TysonNovakCellCycleModel bad_cell_model);
        
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        int num_timesteps = 100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, num_timesteps);// just choosing 5 hours for now - in the Tyson and Novak model cells are yeast and cycle in 75 mins
        
        // cover another exception: create a cell model, delete the time, then
        // try to create another cell model
        std::vector<double> some_proteins(1); // not used except in next line
        TysonNovakCellCycleModel cell_model_1;
        SimulationTime::Destroy();
        TS_ASSERT_THROWS_ANYTHING(TysonNovakCellCycleModel* p_another_cell_model = static_cast<TysonNovakCellCycleModel*>(cell_model_1.CreateCellCycleModel());delete p_another_cell_model;)
        
        //CancerParameters *p_params = CancerParameters::Instance();
        p_simulation_time = SimulationTime::Instance();
        
        double standard_divide_time = 75.19/60.0;
        
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, num_timesteps);// just choosing 5 hours for now - in the Tyson and Novak model cells are yeast and cycle in 75 mins
        TysonNovakCellCycleModel cell_model;
        
        
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            //std::cout << "Time = " << time << "\n";
            bool result = cell_model.ReadyToDivide();
            //std::cout << result << "\n";
            if (time>standard_divide_time)
            {
                TS_ASSERT(result==true);
            }
            else
            {
                TS_ASSERT(result==false);
            }
        }
        
        
        std::vector<double> proteins = cell_model.GetProteinConcentrations();
        
        TS_ASSERT(proteins.size()==6);
        
        TS_ASSERT_DELTA(proteins[0],0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(proteins[1],0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(proteins[2],1.54216806705641, 1e-2);
        TS_ASSERT_DELTA(proteins[3],1.40562614481544, 1e-2);
        TS_ASSERT_DELTA(proteins[4],0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(proteins[5],0.95328206604519, 1e-2);
        
        //double divide_time = p_simulation_time->GetDimensionalisedTime();
        cell_model.ResetModel();
        TysonNovakCellCycleModel *p_cell_model2 = static_cast <TysonNovakCellCycleModel*> (cell_model.CreateCellCycleModel());
        
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            //std::cout << "Time = " << time << "\n";
            bool result = cell_model.ReadyToDivide();
            bool result2 = p_cell_model2->ReadyToDivide();
            //std::cout << result << "\n";
            if (time> 2.0* standard_divide_time)
            {
                TS_ASSERT(result==true);
                TS_ASSERT(result2==true);
            }
            else
            {
                TS_ASSERT(result==false);
                TS_ASSERT(result2==false);
            }
        }
        
        proteins = cell_model.GetProteinConcentrations();
        
        TS_ASSERT(proteins.size()==6);
        
        TS_ASSERT_DELTA(proteins[0],0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(proteins[1],0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(proteins[2],1.54216806705641, 1e-2);
        TS_ASSERT_DELTA(proteins[3],1.40562614481544, 1e-1);
        TS_ASSERT_DELTA(proteins[4],0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(proteins[5],0.95328206604519, 1e-2);
        
        //coverage
        cell_model.SetBirthTime(1.0);
        
        delete p_cell_model2;
        SimulationTime::Destroy();
    }
    
    void TestWntCellCycleModelForVaryingWntStimulus(void) throw(Exception)
    {
        // Here we have a system at rest at Wnt = 1.0 - it would normally go into S phase at 5.971.
        // Instead we reduce Wnt linearly over 0<t<1 to zero and the cell doesn't divide.
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        double endTime = 10.0; //hours
        int numTimesteps = 1000*(int)endTime;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(endTime, numTimesteps);// 15.971 hours to go into S phase
        double wnt_level = 1.0;
        double mutation = 0.0;
        TS_ASSERT_THROWS_NOTHING(WntCellCycleModel cell_model(wnt_level));
        WntCellCycleModel cell_model(wnt_level);
        //TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel wntmodel);
        std::vector<double> cell_cycle_influences;
        cell_cycle_influences.push_back(wnt_level);
        cell_cycle_influences.push_back(mutation);
        for (int i=0; i<numTimesteps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime() ;
            cell_cycle_influences[0] = wnt_level;
            bool result = cell_model.ReadyToDivide(cell_cycle_influences);
            
            if (time <= 1.0)
            {
                wnt_level = 1.0-time;
            }
            else
            {
                wnt_level = 0.0;
            }
            TS_ASSERT(result==false)
        }
        std::vector<double> testResults = cell_model.GetProteinConcentrations();
        TS_ASSERT_DELTA(testResults[0] , 7.330036281693106e-01 , 1e-5);
        TS_ASSERT_DELTA(testResults[1] , 1.715690244022676e-01 , 1e-5);
        TS_ASSERT_DELTA(testResults[2] , 6.127460817296076e-02 , 1e-5);
        TS_ASSERT_DELTA(testResults[3] , 1.549402358669023e-07 , 1e-5);
        TS_ASSERT_DELTA(testResults[4] , 4.579067802591843e-08 , 1e-5);
        TS_ASSERT_DELTA(testResults[5] , 9.999999999999998e-01 , 1e-5);
        TS_ASSERT_DELTA(testResults[6] , 0.5*7.415537855270896e-03 , 1e-5);
        TS_ASSERT_DELTA(testResults[7] , 0.5*7.415537855270896e-03 , 1e-5);
        TS_ASSERT_DELTA(testResults[8] , 0.0 , 1e-6);
        TS_ASSERT_DELTA(testResults[9] , 0.0 , 1e-6);
        
        double diff = 1.0;
        testResults[6] = testResults[6] + diff;
        
        cell_model.SetProteinConcentrationsForTestsOnly(1.0, testResults);
        
        testResults = cell_model.GetProteinConcentrations();
        
        TS_ASSERT_DELTA(testResults[6] , diff + 0.5*7.415537855270896e-03 , 1e-5);
        TS_ASSERT_DELTA(testResults[5] , 9.999999999999998e-01 , 1e-5);
        SimulationTime::Destroy();
        
    }
    
    void TestWntCellCycleModelForAPCSingleHit(void) throw(Exception)
    {
        int num_timesteps = 500;
        double wnt_level = 1.0;
        double mutation = 1.0;
        
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel cell_model_15(wnt_level,1));
        
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase
        WntCellCycleModel cell_model_1(wnt_level);
        SimulationTime::Destroy();
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel *p_cell_model_13 = static_cast<WntCellCycleModel*> (cell_model_1.CreateCellCycleModel()); delete p_cell_model_13;);
        
        p_simulation_time = SimulationTime::Instance();
        
        CancerParameters *p_parameters = CancerParameters::Instance();
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel cell_model_2(wnt_level));
        
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase
        TS_ASSERT_THROWS_NOTHING(WntCellCycleModel cell_model_3(wnt_level));
        
        WntCellCycleModel cell_model(wnt_level,(unsigned)mutation);
        
        
        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 4.804 hrs and then finish dividing
        // 10 hours later at 14.804 hours.
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            //std::cout << "Time = " << time << "\n";
            std::vector <double> cell_cycle_params;
            cell_cycle_params.push_back(wnt_level);
            cell_cycle_params.push_back(mutation);
            bool result = cell_model.ReadyToDivide(cell_cycle_params);
            //std::cout << "divide = " << result << "\n";
            if (time < 4.804+SG2MDuration)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }
        cell_model.ResetModel();
        WntCellCycleModel cell_model_2 = cell_model;
        double second_cycle_start = cell_model_2.GetBirthTime();
        
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime() ;
            //std::cout << "Time = " << time << "\n";
            std::vector <double> cell_cycle_params;
            cell_cycle_params.push_back(wnt_level);
            cell_cycle_params.push_back(mutation);
            bool result = cell_model_2.ReadyToDivide(cell_cycle_params);
            //std::cout << "divide = " << result << "\n";
            if (time< second_cycle_start+4.804+SG2MDuration)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }
        
        SimulationTime::Destroy();
    }
    
    void TestWntCellCycleModelForBetaCatSingleHit(void) throw(Exception)
    {
        int num_timesteps = 500;
        double wnt_level = 0.0;
        double mutation = 2.0;
        
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel cell_model_15(wnt_level));
        
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase
        WntCellCycleModel cell_model_1(wnt_level);
        SimulationTime::Destroy();
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel *p_cell_model_13 = static_cast<WntCellCycleModel*> (cell_model_1.CreateCellCycleModel()); delete p_cell_model_13;);
        
        p_simulation_time = SimulationTime::Instance();
        
        CancerParameters *p_parameters = CancerParameters::Instance();
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel cell_model_2(wnt_level));
        
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase
        TS_ASSERT_THROWS_NOTHING(WntCellCycleModel cell_model_3(wnt_level));
        
        WntCellCycleModel cell_model(wnt_level,(unsigned)mutation);
        
        
        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 7.82 hrs and then finish dividing
        // 10 hours later at 17.82 hours.
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            //std::cout << "Time = " << time << "\n";
            std::vector <double> cell_cycle_params;
            cell_cycle_params.push_back(wnt_level);
            cell_cycle_params.push_back(mutation);
            bool result = cell_model.ReadyToDivide(cell_cycle_params);
            //std::cout << "divide = " << result << "\n";
            if (time < 7.82+SG2MDuration)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }
        cell_model.ResetModel();
        WntCellCycleModel cell_model_2 = cell_model;
        double second_cycle_start = cell_model_2.GetBirthTime();
        
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime() ;
            //std::cout << "Time = " << time << "\n";
            std::vector <double> cell_cycle_params;
            cell_cycle_params.push_back(wnt_level);
            cell_cycle_params.push_back(mutation);
            bool result = cell_model_2.ReadyToDivide(cell_cycle_params);
            //std::cout << "divide = " << result << "\n";
            if (time< second_cycle_start+7.82+SG2MDuration)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }
        
        SimulationTime::Destroy();
    }
    
    void TestWntCellCycleModelForAPCDoubleHit(void) throw(Exception)
    {
        int num_timesteps = 500;
        
        double wnt_level = 0.738;// This shouldn't matter for this kind of cell!
        double mutation = 3;
        
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel cell_model_15(wnt_level,(unsigned)3));
        
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase
        WntCellCycleModel cell_model_1(wnt_level);
        SimulationTime::Destroy();
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel *p_cell_model_13 = static_cast<WntCellCycleModel*> (cell_model_1.CreateCellCycleModel()); delete p_cell_model_13;);
        
        p_simulation_time = SimulationTime::Instance();
        
        CancerParameters *p_parameters = CancerParameters::Instance();
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel cell_model_2(wnt_level,(unsigned)3));
        
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase
        TS_ASSERT_THROWS_NOTHING(WntCellCycleModel cell_model_3(wnt_level));
        
        WntCellCycleModel cell_model(wnt_level,(unsigned) mutation);
        
        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 3.943 hrs and then finish dividing
        // 10 hours later at 13.9435 hours.
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            //std::cout << "Time = " << time << "\n";
            std::vector <double> cell_cycle_params;
            cell_cycle_params.push_back(wnt_level);
            cell_cycle_params.push_back(mutation);
            bool result = cell_model.ReadyToDivide(cell_cycle_params);
            //std::cout << "divide = " << result << "\n";
            if (time < 3.9435+SG2MDuration)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }
        cell_model.ResetModel();
        WntCellCycleModel cell_model_2 = cell_model;
        double second_cycle_start = cell_model_2.GetBirthTime();
        
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime() ;
            //std::cout << "Time = " << time << "\n";
            std::vector <double> cell_cycle_params;
            cell_cycle_params.push_back(wnt_level);
            cell_cycle_params.push_back(mutation);
            bool result = cell_model_2.ReadyToDivide(cell_cycle_params);
            //std::cout << "divide = " << result << "\n";
            if (time< second_cycle_start+3.9435+SG2MDuration)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }
        
        SimulationTime::Destroy();
    }
    
    void TestWntCellCycleModelForConstantWntStimulusHealthyCell(void) throw(Exception)
    {
        int num_timesteps = 500;
        double wnt_level = 1.0;
        double mutation = 0.0;
        
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel cell_model_15(wnt_level));
        
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase
        WntCellCycleModel cell_model_1(wnt_level);
        SimulationTime::Destroy();
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel *p_cell_model_13 = static_cast<WntCellCycleModel*> (cell_model_1.CreateCellCycleModel()); delete p_cell_model_13;);
        
        p_simulation_time = SimulationTime::Instance();
        
        CancerParameters *p_parameters = CancerParameters::Instance();
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel cell_model_2(wnt_level));
        
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase
        TS_ASSERT_THROWS_NOTHING(WntCellCycleModel cell_model_3(wnt_level));
        
        WntCellCycleModel cell_model(wnt_level, (unsigned) mutation);
        
        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 5.971 hrs and then finish dividing
        // 10 hours later at 15.971 hours.
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            //std::cout << "Time = " << time << "\n";
            std::vector <double> cell_cycle_params;
            cell_cycle_params.push_back(wnt_level);
            cell_cycle_params.push_back(mutation);
            bool result = cell_model.ReadyToDivide(cell_cycle_params);
            //std::cout << "divide = " << result << "\n";
            if (time < 5.971+SG2MDuration)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }
        cell_model.ResetModel();
        WntCellCycleModel cell_model_2 = cell_model;
        double second_cycle_start = cell_model_2.GetBirthTime();
        
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime() ;
            //std::cout << "Time = " << time << "\n";
            std::vector <double> cell_cycle_params;
            cell_cycle_params.push_back(wnt_level);
            cell_cycle_params.push_back(mutation);
            bool result = cell_model_2.ReadyToDivide(cell_cycle_params);
            //std::cout << "divide = " << result << "\n";
            if (time< second_cycle_start+5.971+SG2MDuration)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }
        
        SimulationTime::Destroy();
    }
    
    
    void TestArchiveFixedCellCycleModels()
    {
        OutputFileHandler handler("archive");
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "fixed.arch";
        
        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);
            FixedCellCycleModel model;
            p_simulation_time->IncrementTimeOneStep();
            
            model.SetCellType(TRANSIT);
            model.SetBirthTime(-1.0);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            output_arch << static_cast<const SimulationTime&>(*p_simulation_time);
            output_arch << static_cast<const FixedCellCycleModel&>(model);
            
            SimulationTime::Destroy();
        }
        
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            FixedCellCycleModel model;
            model.SetCellType(STEM);
            model.SetBirthTime(-2.0);
            
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore from the archive
            input_arch >> *p_simulation_time;
            input_arch >> model;
            
            // Check
            TS_ASSERT_DELTA(model.GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_EQUALS(model.GetCellType(),TRANSIT);
            TS_ASSERT_DELTA(model.GetAge(),1.5,1e-12);
            
            SimulationTime::Destroy();
        }
    }
    
    void TestArchiveStochasticCellCycleModels()
    {
        OutputFileHandler handler("archive");
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "stoch_cycle.arch";
        
        // Create an ouput archive
        {
            RandomNumberGenerator::Instance();
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);
            StochasticCellCycleModel model;
            p_simulation_time->IncrementTimeOneStep();
            
            model.SetCellType(TRANSIT);
            model.SetBirthTime(-1.0);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            output_arch << static_cast<const SimulationTime&>(*p_simulation_time);
            output_arch << static_cast<const CancerParameters&>(*CancerParameters::Instance());
            output_arch << static_cast<const StochasticCellCycleModel&>(model);
            
            RandomNumberGenerator::Destroy();
            SimulationTime::Destroy();
        }
        
        {
            RandomNumberGenerator::Instance();
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            StochasticCellCycleModel model;
            model.SetCellType(STEM);
            model.SetBirthTime(-2.0);
            
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            CancerParameters *inst1 = CancerParameters::Instance();
            
            inst1->SetSG2MDuration(101.0);
            
            // restore from the archive
            input_arch >> *p_simulation_time;
            input_arch >> *inst1;
            input_arch >> model;
            
            // Check
            TS_ASSERT_DELTA(model.GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_EQUALS(model.GetCellType(),TRANSIT);
            TS_ASSERT_DELTA(model.GetAge(),1.5,1e-12);
            
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
            RandomNumberGenerator::Destroy();
            SimulationTime::Destroy();
        }
    }
    
    void TestArchiveTysonNovakCellCycleModels()
    {
        OutputFileHandler handler("archive");
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "tyson_novak.arch";
        
        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(100.0, 1);
            TysonNovakCellCycleModel model;
            p_simulation_time->IncrementTimeOneStep();
            
            TS_ASSERT_EQUALS(model.ReadyToDivide(),true);
            
            model.SetCellType(TRANSIT);
            model.SetBirthTime(-1.0);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            output_arch << static_cast<const SimulationTime&>(*p_simulation_time);
            output_arch << static_cast<const TysonNovakCellCycleModel&>(model);
            
            SimulationTime::Destroy();
        }
        
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            TysonNovakCellCycleModel model;
            model.SetCellType(STEM);
            model.SetBirthTime(-2.0);
            
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore from the archive
            input_arch >> *p_simulation_time;
            input_arch >> model;
            
            // Check
            TS_ASSERT_EQUALS(model.ReadyToDivide(),true);
            TS_ASSERT_DELTA(model.GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_EQUALS(model.GetCellType(),TRANSIT);
            TS_ASSERT_DELTA(model.GetAge(),101.0,1e-12);
            SimulationTime::Destroy();
        }
    }
    
    void TestArchiveWntCellCycleModels()
    {
        OutputFileHandler handler("archive");
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "wnt.arch";
        
        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(16, 2);
            
            WntCellCycleModel model(1.0);
            p_simulation_time->IncrementTimeOneStep();
            std::vector<double> cell_cycle_influence;
            cell_cycle_influence.push_back(1.0);
            cell_cycle_influence.push_back(0.0);
            TS_ASSERT_EQUALS(model.ReadyToDivide(cell_cycle_influence),false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(model.ReadyToDivide(cell_cycle_influence),true);
            
            model.SetCellType(TRANSIT);
            model.SetBirthTime(-1.0);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            output_arch << static_cast<const SimulationTime&>(*p_simulation_time);
            output_arch << static_cast<const CancerParameters&>(*CancerParameters::Instance());
            output_arch << static_cast<const WntCellCycleModel&>(model);
            
            SimulationTime::Destroy();
        }
        
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            CancerParameters *inst1 = CancerParameters::Instance();
            
            inst1->SetSG2MDuration(101.0);
            
            WntCellCycleModel model(0.0);
            model.SetCellType(STEM);
            model.SetBirthTime(-2.0);
            
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore from the archive
            input_arch >> *p_simulation_time;
            input_arch >> *inst1;
            input_arch >> model;
            
            // Check
            std::vector<double> cell_cycle_influence;
            cell_cycle_influence.push_back(1.0);
            cell_cycle_influence.push_back(0.0);
            TS_ASSERT_EQUALS(model.ReadyToDivide(cell_cycle_influence),true);
            TS_ASSERT_DELTA(model.GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_EQUALS(model.GetCellType(),TRANSIT);
            TS_ASSERT_DELTA(model.GetAge(),17.0,1e-12);
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
            SimulationTime::Destroy();
        }
    }
    
};

#endif /*TESTCELLCYCLEMODELS_HPP_*/
