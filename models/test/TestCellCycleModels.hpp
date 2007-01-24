#ifndef TESTCELLCYCLEMODELS_HPP_
#define TESTCELLCYCLEMODELS_HPP_

#include <cxxtest/TestSuite.h>

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
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0*p_params->GetStemCellCycleTime(), num_steps);
        
        TS_ASSERT_THROWS_NOTHING(FixedCellCycleModel model3);
        
        FixedCellCycleModel our_fixed_stem_cell_cycle_model;
        our_fixed_stem_cell_cycle_model.SetCellType(STEM);
        
        FixedCellCycleModel our_fixed_transit_cell_cycle_model;
        our_fixed_transit_cell_cycle_model.SetCellType(TRANSIT);
        
        FixedCellCycleModel our_fixed_diff_cell_cycle_model;
        our_fixed_diff_cell_cycle_model.SetCellType(DIFFERENTIATED);
        
        for(unsigned i = 0 ; i< num_steps ; i++)
        {
        	p_simulation_time->IncrementTimeOneStep();
        	double time = p_simulation_time->GetDimensionalisedTime();
        	//std::cout << "Time = " << time << " cell age = " << our_fixed_stem_cell_cycle_model.GetAge() << "\n";
        	// Test STEM cells
        	if(time<p_params->GetStemCellCycleTime())
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
        	if(time<p_params->GetTransitCellCycleTime())
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
    	RandomNumberGenerator rand_gen;
    	TS_ASSERT_THROWS_ANYTHING(StochasticCellCycleModel cell_model1(&rand_gen));
    	
        CancerParameters *p_params = CancerParameters::Instance();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        TS_ASSERT_THROWS_ANYTHING(StochasticCellCycleModel cell_model2(&rand_gen));
        
        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0*p_params->GetStemCellCycleTime(), num_steps);
        		
        TS_ASSERT_THROWS_NOTHING(StochasticCellCycleModel cell_model3(&rand_gen));
        
		StochasticCellCycleModel cell_model(&rand_gen);
        
        for(unsigned i = 0 ; i< num_steps ; i++)
        {
        	p_simulation_time->IncrementTimeOneStep();
        	double time = p_simulation_time->GetDimensionalisedTime();
        	//std::cout << "Time = " << time << " cell age = " << our_fixed_stem_cell_cycle_model.GetAge() << "\n";
        	// Test STEM cells
        	cell_model.SetCellType(STEM);
        	if(time<p_params->GetStemCellCycleTime())
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
	        if(time>11.75 && time < 12.25)
	        {
	        	TS_ASSERT(ready_count==54)	
	        } 
        }   
        
		SimulationTime::Destroy();
    }
        
    // Backward Euler solver is still broken so this won't work...
    void no___TestTysonNovakCellCycleModel(void) throw(Exception)
    {
        //CancerParameters *p_params = CancerParameters::Instance();
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        
        int num_timesteps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(5, num_timesteps);// just choosing 5 hours for now - in the Tyson and Novak model cells are yeast and cycle in 75 mins
        TysonNovakCellCycleModel cell_model;
        
        for(int i=0; i<num_timesteps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            std::cout << "Time = " << time << "\n";
            bool result = cell_model.ReadyToDivide();
            std::cout << result << "\n";
            TS_ASSERT(result==false);
        }
        SimulationTime::Destroy();
    }
    
    void TestWntCellCycleModelForConstantWntStimulus(void) throw(Exception)
    {
    	double WntLevel = 1.0;
    	TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel cellModel1(WntLevel));
    	
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        CancerParameters *p_parameters = CancerParameters::Instance();
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel cellModel2(WntLevel));
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        int num_timesteps = 500;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase 
        TS_ASSERT_THROWS_NOTHING(WntCellCycleModel cellModel3(WntLevel));
        
        WntCellCycleModel cellModel(WntLevel);
        
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel wntmodel);
        
        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 5.971 hrs and then finish dividing
        // 10 hours later at 15.971 hours.
        for(int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            //std::cout << "Time = " << time << "\n";
            std::vector <double> cell_cycle_params;
            cell_cycle_params.push_back(WntLevel);
            bool result = cellModel.ReadyToDivide(cell_cycle_params);
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
        
        cellModel.ResetModel();
        
        WntCellCycleModel cellModel2 = cellModel;
        double second_cycle_start = cellModel2.GetBirthTime();
        
        for(int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime() ;
            //std::cout << "Time = " << time << "\n";
            std::vector <double> cell_cycle_params;
            cell_cycle_params.push_back(WntLevel);
            bool result = cellModel2.ReadyToDivide(cell_cycle_params);
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
    
    void TestWntCellCycleModelForVaryingWntStimulus(void) throw(Exception)
    {
        // Here we have a system at rest at Wnt = 1.0 - it would normally go into S phase at 5.971.
        // Instead we reduce Wnt linearly over 0<t<1 to zero and the cell doesn't divide.
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        
        double endTime = 10.0; //hours
        int numTimesteps = 1000*(int)endTime;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(endTime, numTimesteps);// 15.971 hours to go into S phase 
        double WntLevel = 1.0;
        TS_ASSERT_THROWS_NOTHING(WntCellCycleModel cellModel(WntLevel));
        WntCellCycleModel cellModel(WntLevel);
        TS_ASSERT_THROWS_ANYTHING(WntCellCycleModel wntmodel);
        std::vector<double> cellCycleInfluences;
        cellCycleInfluences.push_back(WntLevel);
        for(int i=0; i<numTimesteps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime() ;
			cellCycleInfluences[0] = WntLevel;
            bool result = cellModel.ReadyToDivide(cellCycleInfluences);
            
            if (time <= 1.0)
            {
            	WntLevel = 1.0-time;
            }
            else
            {
            	WntLevel = 0.0;
            }
            TS_ASSERT(result==false)
        }
        std::vector<double> testResults = cellModel.GetProteinConcentrations();
        TS_ASSERT_DELTA(testResults[0] , 7.330036281693106e-01 , 1e-5);
        TS_ASSERT_DELTA(testResults[1] , 1.715690244022676e-01 , 1e-5);
		TS_ASSERT_DELTA(testResults[2] , 6.127460817296076e-02 , 1e-5);
        TS_ASSERT_DELTA(testResults[3] , 1.549402358669023e-07 , 1e-5);
        TS_ASSERT_DELTA(testResults[4] , 4.579067802591843e-08 , 1e-5);
		TS_ASSERT_DELTA(testResults[5] , 9.999999999999998e-01 , 1e-5);
        TS_ASSERT_DELTA(testResults[6] , 7.415537855270896e-03 , 1e-5);
        TS_ASSERT_DELTA(testResults[7] , 0.0 , 1e-6);
        
        double diff = 1.0;
        testResults[6] = testResults[6] + diff;
        
        cellModel.SetProteinConcentrationsForTestsOnly(1.0, testResults);
        
        testResults = cellModel.GetProteinConcentrations();
        
        TS_ASSERT_DELTA(testResults[6] , diff + 7.415537855270896e-03 , 1e-5);
        
        SimulationTime::Destroy();
    }
    
    
};

#endif /*TESTCELLCYCLEMODELS_HPP_*/
