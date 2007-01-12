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
        CancerParameters *p_params = CancerParameters::Instance();
        FixedCellCycleModel our_fixed_cell;
        
        our_fixed_cell.SetCellType(TRANSIT);
        TS_ASSERT(!our_fixed_cell.ReadyToDivide(p_params->GetTransitCellCycleTime()-0.01));
        TS_ASSERT(our_fixed_cell.ReadyToDivide(p_params->GetTransitCellCycleTime()));
        TS_ASSERT(our_fixed_cell.ReadyToDivide(p_params->GetTransitCellCycleTime()+0.01));
        
        our_fixed_cell.SetCellType(STEM);
        TS_ASSERT(!our_fixed_cell.ReadyToDivide(p_params->GetStemCellCycleTime()-0.01));
        TS_ASSERT(our_fixed_cell.ReadyToDivide(p_params->GetStemCellCycleTime()));
        TS_ASSERT(our_fixed_cell.ReadyToDivide(p_params->GetStemCellCycleTime()+0.01));
        
        our_fixed_cell.SetCellType(DIFFERENTIATED);
        TS_ASSERT(!our_fixed_cell.ReadyToDivide(1.0));
        TS_ASSERT(!our_fixed_cell.ReadyToDivide(1e10));
        TS_ASSERT(!our_fixed_cell.ReadyToDivide(1e100));
    }
    
    void TestStochasticCellCycleModel(void) throw(Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator rand_gen;
        StochasticCellCycleModel cell_model(&rand_gen);
        
        cell_model.SetCellType(STEM);
        TS_ASSERT(!cell_model.ReadyToDivide(p_params->GetStemCellCycleTime()-0.01));
        TS_ASSERT(cell_model.ReadyToDivide(p_params->GetStemCellCycleTime()));
        TS_ASSERT(cell_model.ReadyToDivide(p_params->GetStemCellCycleTime()+0.01));
        cell_model.SetCellType(DIFFERENTIATED);
        TS_ASSERT(!cell_model.ReadyToDivide(1.0));
        TS_ASSERT(!cell_model.ReadyToDivide(1e10));
        TS_ASSERT(!cell_model.ReadyToDivide(1e100));
        
        // Testing a random generator is hard...
        cell_model.SetCellType(TRANSIT);
        const int TESTS = 100;
        int ready_count = 0;
        for (int i=0; i<TESTS; i++)
        {
            if (cell_model.ReadyToDivide(p_params->GetTransitCellCycleTime()-0.1))
            {
                ready_count++;
            }
        }
        TS_ASSERT(ready_count>0);
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
            bool result = cell_model.ReadyToDivide(0);
            std::cout << result << "\n";
            TS_ASSERT(result==false);
        }
        SimulationTime::Destroy();
    }
    
    void TestWntCellCycleModelForConstantWntStimulus(void) throw(Exception)
    {
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        
        int num_timesteps = 500;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(20, num_timesteps);// 15.971 hours to go into S phase 
        double WntLevel = 1.0;
        WntCellCycleModel cellModel(WntLevel);
        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 5.971 hrs and then finish dividing
        // 10 hours later at 15.971 hours.
        for(int i=0; i<num_timesteps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime() ;
            //std::cout << "Time = " << time << "\n";
            bool result = cellModel.ReadyToDivide(WntLevel);
            //std::cout << "divide = " << result << "\n";
            if (time<15.971)
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
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(endTime, numTimesteps);// 5.971 hours to go into S phase 
        double WntLevel = 1.0;
        WntCellCycleModel cellModel(WntLevel);
        for(int i=0; i<numTimesteps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime() ;

            bool result = cellModel.ReadyToDivide(WntLevel);
            
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
        SimulationTime::Destroy();
    }
    
    
};

#endif /*TESTCELLCYCLEMODELS_HPP_*/
