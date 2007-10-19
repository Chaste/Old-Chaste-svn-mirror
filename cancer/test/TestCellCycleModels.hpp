#ifndef TESTCELLCYCLEMODELS_HPP_
#define TESTCELLCYCLEMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "ConformingTetrahedralMesh.cpp"
#include "CellsGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "IngeWntSwatCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "OxygenBasedCellCycleModel.hpp"
#include "WntGradient.hpp"

class TestCellCycleModels : public CxxTest::TestSuite
{
public:
    void TestFixedCellCycleModel(void) throw(Exception)
    {   
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        unsigned num_steps = 100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            2.0*(p_params->GetStemCellG1Duration()
                  +p_params->GetSG2MDuration()     ), num_steps);
        
        TS_ASSERT_THROWS_NOTHING(FixedCellCycleModel model3);
        
        FixedCellCycleModel* p_our_fixed_stem_cell_cycle_model = new FixedCellCycleModel;
        TissueCell stem_cell(STEM, // type
                           HEALTHY,//Mutation State
                           0,  // generation
                           p_our_fixed_stem_cell_cycle_model);
        
        
        TS_ASSERT_EQUALS(stem_cell.GetCellType(),STEM);
        
        FixedCellCycleModel* p_our_fixed_transit_cell_cycle_model = new FixedCellCycleModel;
        TissueCell transit_cell(TRANSIT, // type
                           HEALTHY,//Mutation State
                           0,  // generation
                           p_our_fixed_transit_cell_cycle_model);
                           
        TS_ASSERT_EQUALS(transit_cell.GetCellType(),TRANSIT);
        
        FixedCellCycleModel* p_our_fixed_diff_cell_cycle_model = new FixedCellCycleModel;
        TissueCell diff_cell(DIFFERENTIATED, // type
                           HEALTHY,//Mutation State
                           0,  // generation
                           p_our_fixed_diff_cell_cycle_model);
                           
        TS_ASSERT_EQUALS(diff_cell.GetCellType(),DIFFERENTIATED);
        
        FixedCellCycleModel* p_our_fixed_hepa_one_cell_cycle_model = new FixedCellCycleModel;
        TissueCell hepa_one_cell(HEPA_ONE, // type
                           HEALTHY,//Mutation State
                           0,  // generation
                           p_our_fixed_hepa_one_cell_cycle_model);
                
        for (unsigned i = 0 ; i< num_steps ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();

            // Test STEM cells
            if (time<p_params->GetStemCellG1Duration() + p_params->GetSG2MDuration())
            {
                TS_ASSERT(!p_our_fixed_stem_cell_cycle_model->ReadyToDivide());
            }
            else
            {
                TS_ASSERT(p_our_fixed_stem_cell_cycle_model->ReadyToDivide());
            }
            // Test a Transit Cell
            if (time<p_params->GetTransitCellG1Duration() + p_params->GetSG2MDuration())
            {
                TS_ASSERT(!p_our_fixed_transit_cell_cycle_model->ReadyToDivide());
            }
            else
            {
                TS_ASSERT(p_our_fixed_transit_cell_cycle_model->ReadyToDivide());
            }
            // Test a DIFFERENTIATED cell
            TS_ASSERT(!p_our_fixed_diff_cell_cycle_model->ReadyToDivide());
            // Test a HEPA_ONE cell
            if (time<p_params->GetHepaOneCellG1Duration() + p_params->GetSG2MDuration())
            {
                TS_ASSERT(!p_our_fixed_hepa_one_cell_cycle_model->ReadyToDivide());
            }
            else
            {
                TS_ASSERT(p_our_fixed_hepa_one_cell_cycle_model->ReadyToDivide());
            }
        }
        
        TS_ASSERT_DELTA(p_our_fixed_stem_cell_cycle_model->GetAge(), p_simulation_time->GetDimensionalisedTime(), 1e-9);
        TS_ASSERT_DELTA(p_our_fixed_transit_cell_cycle_model->GetAge(), p_simulation_time->GetDimensionalisedTime(), 1e-9);
        TS_ASSERT_DELTA(p_our_fixed_diff_cell_cycle_model->GetAge(), p_simulation_time->GetDimensionalisedTime(), 1e-9);
        TS_ASSERT_DELTA(p_our_fixed_hepa_one_cell_cycle_model->GetAge(), p_simulation_time->GetDimensionalisedTime(), 1e-9);
        SimulationTime::Destroy();
    }
    
    void TestStochasticCellCycleModel(void) throw(Exception)
    {
        RandomNumberGenerator::Instance()->Reseed(0);
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();

        SimulationTime* p_simulation_time = SimulationTime::Instance();        
        unsigned num_steps = 100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            2.0*(p_params->GetStemCellG1Duration() + p_params->GetSG2MDuration()), num_steps);
        
        TS_ASSERT_THROWS_NOTHING(StochasticCellCycleModel cell_model3);
        
        StochasticCellCycleModel* p_stem_model = new StochasticCellCycleModel;
        StochasticCellCycleModel* p_transit_model = new StochasticCellCycleModel;
        
        TissueCell stem_cell(STEM, // type
                              HEALTHY,//Mutation State
                              0,  // generation
                              p_stem_model);
        
        TissueCell transit_cell(TRANSIT, // type
                              HEALTHY,//Mutation State
                              0,  // generation
                              p_transit_model);        
                              
        for (unsigned i = 0 ; i< num_steps ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
         
            // Test STEM cells
            if (time<p_params->GetStemCellG1Duration()+ p_params->GetSG2MDuration())
            {
                TS_ASSERT_EQUALS(p_stem_model->ReadyToDivide(), false);
            }
            else
            {
                TS_ASSERT_EQUALS(p_stem_model->ReadyToDivide(), true);
            }
            // Test Transit cells
            if (time < 11.0676) // this is first random normal deviate of transit cell cycle time.
            {
                TS_ASSERT_EQUALS(p_transit_model->ReadyToDivide(),false);
            }
            else
            {
                TS_ASSERT_EQUALS(p_transit_model->ReadyToDivide(),true);
            }
        }
        RandomNumberGenerator::Destroy();
        SimulationTime::Destroy();
    }
    
    
    void TestTysonNovakCellCycleModel(void) throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        int num_timesteps = 100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, num_timesteps);// just choosing 5 hours for now - in the Tyson and Novak model cells are yeast and cycle in 75 mins
                
        double standard_divide_time = 75.19/60.0;
        
        TysonNovakCellCycleModel* p_cell_model = new TysonNovakCellCycleModel;
        //coverage
        p_cell_model->SetBirthTime(p_simulation_time->GetDimensionalisedTime());           
        TissueCell cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,  // generation
                                   p_cell_model);
                                   
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();

            bool result = p_cell_model->ReadyToDivide();

            if (time>standard_divide_time)
            {
                TS_ASSERT(result==true);
            }
            else
            {
                TS_ASSERT(result==false);
            }
        }
        
        std::vector<double> proteins = p_cell_model->GetProteinConcentrations();
        
        TS_ASSERT(proteins.size()==6);
        
        TS_ASSERT_DELTA(proteins[0],0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(proteins[1],0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(proteins[2],1.54216806705641, 1e-2);
        TS_ASSERT_DELTA(proteins[3],1.40562614481544, 1e-2);
        TS_ASSERT_DELTA(proteins[4],0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(proteins[5],0.95328206604519, 1e-2);
        
        //double divide_time = p_simulation_time->GetDimensionalisedTime();
        p_cell_model->ResetModel();
        TysonNovakCellCycleModel *p_cell_model2 = static_cast <TysonNovakCellCycleModel*> (p_cell_model->CreateCellCycleModel());
        
        TissueCell stem_cell_2(STEM, // type
                                     APC_ONE_HIT,//Mutation State
                                     0,  // generation
                                     p_cell_model2);
        
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();

            bool result = p_cell_model->ReadyToDivide();
            bool result2 = p_cell_model2->ReadyToDivide();

            if (time> 2.0* standard_divide_time)
            {
                TS_ASSERT_EQUALS(result,true);
                TS_ASSERT_EQUALS(result2,true);
            }
            else
            {
                TS_ASSERT_EQUALS(result,false);
                TS_ASSERT_EQUALS(result2,false);
            }
        }
        
        proteins = p_cell_model->GetProteinConcentrations();
        
        TS_ASSERT_EQUALS(proteins.size(),6u);
        
        TS_ASSERT_DELTA(proteins[0],0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(proteins[1],0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(proteins[2],1.54216806705641, 1e-2);
        TS_ASSERT_DELTA(proteins[3],1.40562614481544, 1e-1);
        TS_ASSERT_DELTA(proteins[4],0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(proteins[5],0.95328206604519, 1e-2);
        
        SimulationTime::Destroy();
    }
    
    void TestWntCellCycleModelForVaryingWntStimulus() throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        
        // Here we have a system at rest at Wnt = 1.0 - it would normally go into S phase at 5.971.
        // Instead we reduce Wnt linearly over 0<t<1 to zero and the cell doesn't divide.
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        
        double end_time = 10.0; //hours
        int num_timesteps = 1000*(int)end_time;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);// 15.971 hours to go into S phase
                
        double wnt_level = 1.0;
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);

        WntCellCycleModel* p_cell_model = new WntCellCycleModel();
        
        TissueCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,  // generation
                                   p_cell_model);
                           
        stem_cell.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(stem_cell.GetCellType(),TRANSIT);
        
        for (int i=0; i<num_timesteps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();            
            bool result = p_cell_model->ReadyToDivide();
            
            if (time <= 1.0)
            {
                wnt_level = 1.0-time;
            }
            else
            {
                wnt_level = 0.0;
            }
            WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
            
            TS_ASSERT(result==false)
        }
        
        std::vector<double> test_results = p_cell_model->GetProteinConcentrations();
        TS_ASSERT_DELTA(test_results[0] , 7.330036281693106e-01 , 1e-5);
        TS_ASSERT_DELTA(test_results[1] , 1.715690244022676e-01 , 1e-5);
        TS_ASSERT_DELTA(test_results[2] , 6.127460817296076e-02 , 1e-5);
        TS_ASSERT_DELTA(test_results[3] , 1.549402358669023e-07 , 1e-5);
        TS_ASSERT_DELTA(test_results[4] , 4.579067802591843e-08 , 1e-5);
        TS_ASSERT_DELTA(test_results[5] , 9.999999999999998e-01 , 1e-5);
        TS_ASSERT_DELTA(test_results[6] , 0.5*7.415537855270896e-03 , 1e-5);
        TS_ASSERT_DELTA(test_results[7] , 0.5*7.415537855270896e-03 , 1e-5);
        TS_ASSERT_DELTA(test_results[8] , 0.0 , 1e-6);
        
        TS_ASSERT_EQUALS(stem_cell.GetCellType(), DIFFERENTIATED);
        
        double diff = 1.0;
        test_results[6] = test_results[6] + diff;
        
        p_cell_model->SetProteinConcentrationsForTestsOnly(1.0, test_results);
        
        test_results = p_cell_model->GetProteinConcentrations();
        
        TS_ASSERT_DELTA(test_results[6] , diff + 0.5*7.415537855270896e-03 , 1e-5);
        TS_ASSERT_DELTA(test_results[5] , 9.999999999999998e-01 , 1e-5);
                
        SimulationTime::Destroy();        
        WntGradient::Destroy();
    }
    
    void TestIngeWntSwatCellCycleModel() throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        
        // Here we have a system at rest at Wnt = 1.0 - it would normally go into S phase at 5.971.
        // Instead we reduce Wnt linearly over 0<t<1 to zero and the cell doesn't divide.
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        
        double end_time = 30.0; //hours
        int num_timesteps = 1000*(int)end_time;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);// 15.971 hours to go into S phase
                
        double wnt_level = 1.0;
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);

        // cover exception
        TS_ASSERT_THROWS_ANYTHING(IngeWntSwatCellCycleModel model(0));

        IngeWntSwatCellCycleModel* p_cell_model = new IngeWntSwatCellCycleModel(1);
        
        TS_ASSERT_EQUALS(p_cell_model->GetHypothesis(), 1u);
        
        TissueCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,  // generation
                                   p_cell_model);
                                   
        // Coverage of cell cycle model copying without an ODE system set up
        TissueCell stem_cell2 = stem_cell;
        TS_ASSERT_EQUALS(stem_cell2.GetMutationState(), HEALTHY);
                           
        stem_cell.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(stem_cell.GetCellType(),TRANSIT);
        WntGradient::Instance()->SetConstantWntValueForTesting(1.0);
        for (int i=0; i<2*num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            bool result = stem_cell.ReadyToDivide();
            if (SimulationTime::Instance()->GetDimensionalisedTime()
                    >=  6.1877 + CancerParameters::Instance()->GetSG2MDuration())
            {
                TS_ASSERT_EQUALS(result, true);
            }
            else
            {
                TS_ASSERT_EQUALS(result, false);
            }
        }
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetDimensionalisedTime(), 20.0, 1e-4);
        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), true);
        
        std::vector<double> test_results = p_cell_model->GetProteinConcentrations();

        TS_ASSERT_DELTA(test_results[0] , 2.937528298476724e-01, 1e-3);
        TS_ASSERT_DELTA(test_results[1] , 1.000000000000000e+00, 1e-2);
        TS_ASSERT_DELTA(test_results[2] , 2.400571020059542e+00, 1e-3);
        TS_ASSERT_DELTA(test_results[3] , 1.392891502782046e+00, 1e-3);
        TS_ASSERT_DELTA(test_results[4] , 1.358708742498684e-01, 1e-3);
        TS_ASSERT_DELTA(test_results[5] , 1.428571428571429e-01, 1e-4);
        TS_ASSERT_DELTA(test_results[6] , 2.857142857142857e-02, 1e-4);
        TS_ASSERT_DELTA(test_results[7] , 2.120643654085205e-01, 1e-4);
        TS_ASSERT_DELTA(test_results[8] , 1.439678172957377e+01, 1e-3);
        TS_ASSERT_DELTA(test_results[9] ,                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[10],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[11],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[12], 1.000000000000002e+01, 1e-4);
        TS_ASSERT_DELTA(test_results[13], 1.028341552112414e+02, 1e-2);
        TS_ASSERT_DELTA(test_results[14],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[15], 2.499999999999999e+01, 1e-4);
        TS_ASSERT_DELTA(test_results[16], 1.439678172957377e+01, 1e-3);
        TS_ASSERT_DELTA(test_results[17],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[18],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[19],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[20], 2.235636835087684e+00, 1e-4);
        TS_ASSERT_DELTA(test_results[21], 1.000000000000000e+00, 1e-4);
        
        // Acts as if it was divided at time = 16.1877... which is OK 
        // (cell cycle model dictates division time, not when the cell is manually
        // divided)
        TissueCell daughter_cell = stem_cell.Divide();
        AbstractCellCycleModel* p_cell_model2 = daughter_cell.GetCellCycleModel();
        
        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_cell_model2->ReadyToDivide(), false);
        
        // Check first 5 protein levels have been reset and the rest are the same.
        test_results = p_cell_model->GetProteinConcentrations();
        TS_ASSERT_DELTA(test_results[0] , 2.631865125420296e-01, 1e-3);
        TS_ASSERT_DELTA(test_results[1] , 2.678271949808561e-01, 1e-3);
        TS_ASSERT_DELTA(test_results[2] , 2.389956081120099e+00, 1e-3);
        TS_ASSERT_DELTA(test_results[3] , 1.390258620103223e+00, 1e-3);
        TS_ASSERT_DELTA(test_results[4] , 1.218603203963113e-01, 1e-3);
        TS_ASSERT_DELTA(test_results[5] , 1.428571428571429e-01, 1e-4);
        TS_ASSERT_DELTA(test_results[6] , 2.857142857142857e-02, 1e-4);
        TS_ASSERT_DELTA(test_results[7] , 2.120643654085205e-01, 1e-4);
        TS_ASSERT_DELTA(test_results[8] , 1.439678172957377e+01, 1e-3);
        TS_ASSERT_DELTA(test_results[9] ,                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[10],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[11],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[12], 1.000000000000002e+01, 1e-4);
        TS_ASSERT_DELTA(test_results[13], 1.028341552112414e+02, 1e-2);
        TS_ASSERT_DELTA(test_results[14],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[15], 2.499999999999999e+01, 1e-4);
        TS_ASSERT_DELTA(test_results[16], 1.439678172957377e+01, 1e-3);
        TS_ASSERT_DELTA(test_results[17],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[18],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[19],                     0, 1e-4);
        TS_ASSERT_DELTA(test_results[20], 2.235636835087684e+00, 1e-4);
        TS_ASSERT_DELTA(test_results[21], 1.000000000000000e+00, 1e-4);
        
        for (int i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();            
            bool result = p_cell_model->ReadyToDivide();
            bool result2 = p_cell_model2->ReadyToDivide();
            
            if (time <= 21.0)
            {
                wnt_level = 21.0 - time;
            }
            else
            {
                wnt_level = 0.0;
            }
            WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
            
            TS_ASSERT_EQUALS(result, false);
            TS_ASSERT_EQUALS(result2, false);
        }
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetDimensionalisedTime(), 30.0, 1e-4);
        
        test_results = p_cell_model->GetProteinConcentrations();       
        
        TS_ASSERT_DELTA(test_results[0] , 0.3636, 1e-3);
        TS_ASSERT_DELTA(test_results[1] , 0.9900, 1e-2);
        TS_ASSERT_DELTA(test_results[2] , 1.4996, 1e-3);
        TS_ASSERT_DELTA(test_results[3] , 0.8180, 1e-4);
        TS_ASSERT_DELTA(test_results[4] , 0.1000, 1e-4);
        TS_ASSERT_DELTA(test_results[5] , 0.6666, 1e-4);
        TS_ASSERT_DELTA(test_results[6] , 0.0666, 1e-4);
        TS_ASSERT_DELTA(test_results[7] , 0.7311, 1e-3);
        TS_ASSERT_DELTA(test_results[8] , 6.0481, 1e-3);
        TS_ASSERT_DELTA(test_results[9] , 0, 1e-4);
        TS_ASSERT_DELTA(test_results[10], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[11], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[12], 17.5432, 1e-4);
        TS_ASSERT_DELTA(test_results[13], 75.8316, 1e-2);
        TS_ASSERT_DELTA(test_results[14], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[15], 29.5727, 1e-4);
        TS_ASSERT_DELTA(test_results[16], 7.1564, 1e-3);
        TS_ASSERT_DELTA(test_results[17], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[18], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[19], 0, 1e-4);
        TS_ASSERT_DELTA(test_results[20], 1.6047, 1e-4);
        TS_ASSERT_DELTA(test_results[21], 0.0000, 1e-4);

        TS_ASSERT_EQUALS(stem_cell.GetCellType(), DIFFERENTIATED);

        // membrane_beta_cat = Ca + Ma
        double membrane_beta_cat = test_results[13]+test_results[14];

        // cytoplasmic_beta_cat = Cu + Co + Cc + Mo + Mc
        double cytoplasm_beta_cat = test_results[7] + test_results[8]
                          + test_results[9] + test_results[10]+test_results[11];

        // nuclear_beta_cat = Cot + Cct + Mot + Mct
        double nuclear_beta_cat = test_results[16] + test_results[17]
                                    + test_results[18] + test_results[19];

        TS_ASSERT_DELTA(p_cell_model->GetMembraneBoundBetaCateninLevel(), membrane_beta_cat, 1e-4);
        TS_ASSERT_DELTA(p_cell_model->GetCytoplasmicBetaCateninLevel(), cytoplasm_beta_cat, 1e-4);
        TS_ASSERT_DELTA(p_cell_model->GetNuclearBetaCateninLevel(), nuclear_beta_cat, 1e-4);
        
        SimulationTime::Destroy();        
        WntGradient::Destroy();
    }
    
    void TestWntCellCycleModelForAPCSingleHit(void) throw(Exception)
    {
        int num_timesteps = 500;

        CancerParameters::Instance()->Reset();
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase
        
        double wnt_level = 1.0;        
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);

        WntCellCycleModel* p_cell_model = new WntCellCycleModel();
                
        TissueCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,  // generation
                                   p_cell_model);

        stem_cell.InitialiseCellCycleModel();
                                                              
        double SG2MDuration = CancerParameters::Instance()->GetSG2MDuration();
        TS_ASSERT_THROWS_NOTHING(WntCellCycleModel cell_model_3());
        
        WntCellCycleModel* p_cell_model_1 = new WntCellCycleModel();
                
        TissueCell stem_cell_1(STEM, // type
                                     APC_ONE_HIT,//Mutation State
                                     0,  // generation
                                     p_cell_model_1);
        stem_cell_1.InitialiseCellCycleModel();
        
        // Wnt cells not set up to deal with unknown mutations...
        TissueCell cell(STEM, // type
                        ALARCON_NORMAL,//Mutation State
                        0,  // generation
                        new WntCellCycleModel());
        TS_ASSERT_THROWS_ANYTHING(cell.InitialiseCellCycleModel());
                
        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 4.804 hrs and then finish dividing
        // 10 hours later at 14.804 hours.
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();            
            bool result = p_cell_model_1->ReadyToDivide();
            if (time < 4.804+SG2MDuration)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }
        p_cell_model_1->ResetModel();
        double second_cycle_start = p_cell_model_1->GetBirthTime();
        
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();            
            bool result = p_cell_model_1->ReadyToDivide();
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
        WntGradient::Destroy();
    }
    
    void TestWntCellCycleModelForBetaCatSingleHit(void) throw(Exception)
    {
        CancerParameters::Instance()->Reset();

        int num_timesteps = 500;
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase

        double wnt_level = 0.0;
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
        
        WntCellCycleModel* p_cell_model = new WntCellCycleModel();
                
        TissueCell stem_cell(STEM, // type
                                   BETA_CATENIN_ONE_HIT,//Mutation State
                                   0,  // generation
                                   p_cell_model);
        stem_cell.InitialiseCellCycleModel();
                        
        CancerParameters *p_parameters = CancerParameters::Instance();
        
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        TS_ASSERT_THROWS_NOTHING(WntCellCycleModel cell_model_3());
        
        WntCellCycleModel* p_cell_model_1 = new WntCellCycleModel();
                
        TissueCell stem_cell_1(STEM, // type
                                     BETA_CATENIN_ONE_HIT,//Mutation State
                                     0,  // generation
                                     p_cell_model_1);        
        stem_cell_1.InitialiseCellCycleModel();
        
        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 7.82 hrs and then finish dividing
        // 10 hours later at 17.82 hours.
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            TS_ASSERT_THROWS_ANYTHING(p_cell_model_1->UpdateCellType());
            bool result = p_cell_model_1->ReadyToDivide();
            
            if (time < 7.82+SG2MDuration)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }
        
        p_cell_model_1->ResetModel();
        double second_cycle_start = p_cell_model_1->GetBirthTime();
        
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();            
            bool result = p_cell_model_1->ReadyToDivide();

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
        WntGradient::Destroy();
    }
    
    void TestWntCellCycleModelForAPCDoubleHit(void) throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        
        int num_timesteps = 500;
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase
        
        double wnt_level = 0.738;// This shouldn't matter for this kind of cell!
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
                
        WntCellCycleModel* p_cell_model_1 = new WntCellCycleModel();
        
        TissueCell stem_cell_1(STEM, // type
                                     APC_TWO_HIT,//Mutation State
                                     0,  // generation
                                     p_cell_model_1);   
        stem_cell_1.InitialiseCellCycleModel();
                
        CancerParameters *p_parameters = CancerParameters::Instance();
        
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        WntCellCycleModel* p_cell_model_2 = new WntCellCycleModel();
                
        TissueCell stem_cell_2(STEM, // type
                                     APC_TWO_HIT,//Mutation State
                                     0,  // generation
                                     p_cell_model_2);   
        stem_cell_2.InitialiseCellCycleModel();
        
        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 3.943 hrs and then finish dividing
        // 10 hours later at 13.9435 hours.
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();
            
            bool result = p_cell_model_2->ReadyToDivide();
            
            if (time < 3.9435+SG2MDuration)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }
        p_cell_model_2->ResetModel();
        double second_cycle_start = p_cell_model_2->GetBirthTime();
        
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime() ;
            bool result = p_cell_model_2->ReadyToDivide();

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
        WntGradient::Destroy();
    }
    
    void TestWntCellCycleModelForConstantWntStimulusHealthyCell(void) throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        int num_timesteps = 500;
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(40, num_timesteps);// 15.971 hours to go into S phase
        
        double wnt_level = 1.0;
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
                
        WntCellCycleModel* p_cell_model_1 = new WntCellCycleModel();
                
        TissueCell stem_cell_1(STEM, // type
                                     HEALTHY,//Mutation State
                                     0,  // generation
                                     p_cell_model_1);   
        stem_cell_1.InitialiseCellCycleModel();
        
        CancerParameters *p_parameters = CancerParameters::Instance();
         
        double SG2MDuration = p_parameters->GetSG2MDuration();
        
        WntCellCycleModel* p_cell_model_2 = new WntCellCycleModel();
                
        TissueCell stem_cell_2(STEM, // type
                                     HEALTHY,//Mutation State
                                     0,  // generation
                                     p_cell_model_2);   
        stem_cell_2.InitialiseCellCycleModel();
        
        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 5.971 hrs and then finish dividing
        // 10 hours later at 15.971 hours.
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();            
            bool result = p_cell_model_2->ReadyToDivide();
            
            if (time < 5.971+SG2MDuration)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }
        p_cell_model_2->ResetModel();
        double second_cycle_start = p_cell_model_2->GetBirthTime();
        
        for (int i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime() ;
            bool result = p_cell_model_2->ReadyToDivide();

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
        WntGradient::Destroy();
    }
    
    void TestStochasticWntCellCycleModel() throw (Exception)
    {
        RandomNumberGenerator::Instance()->Reseed(0);
        
        int num_timesteps = 100;
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(20, num_timesteps);// 15.971 hours to go into S phase

        double wnt_level = 1.0;
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
        
        StochasticWntCellCycleModel* p_cell_model = new StochasticWntCellCycleModel();
                
        TissueCell stem_cell(STEM, // type
                                   HEALTHY,//Mutation State
                                   0,  // generation
                                   p_cell_model);   
        stem_cell.InitialiseCellCycleModel();
               
        // A WntCellCycleModel does this:
        // Run the Wnt model for a full constant Wnt stimulus for 20 hours.
        // Model should enter S phase at 5.971 hrs and then finish dividing
        // 10 hours later at 15.971 hours.
        // 
        // A StochasticWntCellCycleModel does this:
        // divides at the same time with a random normal distribution 
        // for the SG2M time (default 10) in this case 9.0676
        
        for (int i=0; i<num_timesteps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();            
            bool result = p_cell_model->ReadyToDivide();
            
            if (time < 5.971 + 9.0676)
            {
                TS_ASSERT(result==false);
            }
            else
            {
                TS_ASSERT(result==true);
            }
        }
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
    }
    
    void TestSimpleWntCellCycleModel() throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        // Set up the simulation time
        SimulationTime *p_simulation_time = SimulationTime::Instance();        
        double end_time = 60.0; 
        unsigned num_timesteps = 1000*(unsigned)end_time;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);
        
        double wnt_level = 1.0;
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
        
        /* First three random cell cycle times are
         * 11.2712
         * 13.1632
         * 11.0676
         */
        double cycle_time_one = 11.0676;
        double cycle_time_two = 13.1632;
        double cycle_time_three = 11.2712;
        SimpleWntCellCycleModel* p_cycle_model = new SimpleWntCellCycleModel;
        TissueCell cell(STEM, HEALTHY, 0, p_cycle_model);
                
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();            
            bool result = p_cycle_model->ReadyToDivide();
            
            if (time < cycle_time_one)
            {
                TS_ASSERT_EQUALS(result, false);
            }
            else
            {
                TS_ASSERT_EQUALS(result, true);
            }
            TS_ASSERT_EQUALS(cell.GetCellType(), TRANSIT);
        }
        
        // divide the cell
        double division_time = SimulationTime::Instance()->GetDimensionalisedTime();
        p_cycle_model->ResetModel();
        SimpleWntCellCycleModel *p_cycle_model2 = static_cast <SimpleWntCellCycleModel*> (p_cycle_model->CreateCellCycleModel());
        
        TissueCell cell2(STEM, APC_TWO_HIT, 0, p_cycle_model2);
        cell.SetMutationState(LABELLED);
        // Reduce Wnt gradient
        wnt_level = 0.6;
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
        
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();            
            bool result = p_cycle_model->ReadyToDivide();
            bool result2 = p_cycle_model2->ReadyToDivide();
            
            if (time < division_time + cycle_time_two)
            {
                TS_ASSERT_EQUALS(result, false);
            }
            else
            {
                TS_ASSERT_EQUALS(result, true);
            }
            
            if (time < division_time + cycle_time_three)
            {
                TS_ASSERT_EQUALS(result2, false);
            }
            else
            {
                TS_ASSERT_EQUALS(result2, true);
            }
            TS_ASSERT_EQUALS(cell.GetCellType(), TRANSIT);
            TS_ASSERT_EQUALS(cell2.GetCellType(), TRANSIT);
        }
        
        p_cycle_model->ResetModel();
        p_cycle_model2->ResetModel();
        division_time = SimulationTime::Instance()->GetDimensionalisedTime();
        // double cycle_time_four = 11.2204
        double cycle_time_five = 10.747;
        // Reduce Wnt gradient so only beta-cat or APC2 hit cells divide.
        wnt_level = 0.15;
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
        
        cell.SetMutationState(APC_ONE_HIT);
        cell2.SetMutationState(BETA_CATENIN_ONE_HIT);
        
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetDimensionalisedTime();            
            bool result = p_cycle_model->ReadyToDivide();
            bool result2 = p_cycle_model2->ReadyToDivide();
            
            TS_ASSERT_EQUALS(result, false);    // these mutants under wnt threshold
            
            if (time < division_time + cycle_time_five)
            {
                TS_ASSERT_EQUALS(result2, false);
            }
            else
            {
                TS_ASSERT_EQUALS(result2, true);
            }
            TS_ASSERT_EQUALS(cell.GetCellType(), DIFFERENTIATED);
            TS_ASSERT_EQUALS(cell2.GetCellType(), TRANSIT);
        }
        
        RandomNumberGenerator::Destroy();
        SimulationTime::Destroy();
        WntGradient::Destroy();
    }
    
    void TestOxygenBasedCellCycleModel() throw(Exception)
    {        
        CancerParameters::Instance()->Reset();
       
        // set up SimulationTime         
        SimulationTime *p_simulation_time = SimulationTime::Instance();   
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(20.0, 2);
        
        // set up oxygen_concentration        
        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        // create model
        OxygenBasedCellCycleModel* p_cell_model = new OxygenBasedCellCycleModel();
        
        // create cell 
        TissueCell cell(HEPA_ONE, ALARCON_NORMAL, 0, p_cell_model);    
        
        // Coverage of cell cycle model copying without an ODE system set up
        TissueCell stem_cell2 = cell;
        TS_ASSERT_EQUALS(stem_cell2.GetMutationState(), ALARCON_NORMAL);
                               
        cell.InitialiseCellCycleModel();
        
        // check oxygen concentration is correct in cell cycle model
        TS_ASSERT_DELTA(p_cell_model->GetProteinConcentrations()[5], 1.0, 1e-5);
        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), false);        
        
        // divide a cell    
        OxygenBasedCellCycleModel *p_cell_model2 = static_cast <OxygenBasedCellCycleModel*> (p_cell_model->CreateCellCycleModel());
        
        TissueCell cell2(HEPA_ONE, ALARCON_NORMAL, 0, p_cell_model2);
        
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),false)
        TS_ASSERT_EQUALS(p_cell_model2->ReadyToDivide(),false);
        
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),true)
        TS_ASSERT_EQUALS(p_cell_model2->ReadyToDivide(),true);
        
        TS_ASSERT_THROWS_NOTHING(p_cell_model->ResetModel());     

        SimulationTime::Destroy();
        CellwiseData<2>::Destroy();
    }
        
    void TestArchiveFixedCellCycleModels() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "fixed.arch";
        
        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);
            FixedCellCycleModel model;
            p_simulation_time->IncrementTimeOneStep();
            
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
            model.SetBirthTime(-2.0);
            
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore from the archive
            input_arch >> *p_simulation_time;
            input_arch >> model;
            
            // Check
            TS_ASSERT_DELTA(model.GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_DELTA(model.GetAge(),1.5,1e-12);
            
            SimulationTime::Destroy();
        }
    }
    
    void TestArchiveStochasticCellCycleModels()
    {
        CancerParameters::Instance()->Reset();

        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "stoch_cycle.arch";
        
        double random_number_test = 0;
        
        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);
            
            StochasticCellCycleModel model;
            p_simulation_time->IncrementTimeOneStep();
            model.SetBirthTime(-1.0);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            output_arch << static_cast<const SimulationTime&>(*p_simulation_time);
            output_arch << static_cast<const StochasticCellCycleModel&>(model);
            
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            random_number_test = p_gen->ranf();
            
            RandomNumberGenerator::Destroy();
            SimulationTime::Destroy();
        }
        
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            p_gen->Reseed(36);
            
            StochasticCellCycleModel model;
            model.SetBirthTime(-2.0);
            
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            CancerParameters *inst1 = CancerParameters::Instance();
            
            inst1->SetSDuration(101.0);
            
            // restore from the archive
            input_arch >> *p_simulation_time;
            input_arch >> model;
            
            // Check
            TS_ASSERT_DELTA(model.GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_DELTA(model.GetAge(),1.5,1e-12);
            TS_ASSERT_DELTA(p_gen->ranf(),random_number_test,1e-7);
            
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
            
            RandomNumberGenerator::Destroy();
            SimulationTime::Destroy();
        }
    }
    
    void TestArchiveTysonNovakCellCycleModels()
    {
        CancerParameters::Instance()->Reset();

        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "tyson_novak.arch";
        
        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(100.0, 1);
            TysonNovakCellCycleModel* p_model = new TysonNovakCellCycleModel;
            p_simulation_time->IncrementTimeOneStep();
            
            TissueCell cell(TRANSIT, // type
                            HEALTHY,//Mutation State
                            0,  // generation
                            p_model);
            cell.InitialiseCellCycleModel();  
            
            TS_ASSERT_EQUALS(p_model->ReadyToDivide(),true);
            
            p_model->SetBirthTime(-1.0);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            TissueCell* const p_cell = &cell;
            
            output_arch << static_cast<const SimulationTime&>(*p_simulation_time);
            output_arch << p_cell;
            
            SimulationTime::Destroy();
        }
                
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            TissueCell* p_cell;            
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore from the archive
            input_arch >> *p_simulation_time;
            input_arch >> p_cell;
            
            AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();
            
            // Check
            TS_ASSERT_EQUALS(p_model->ReadyToDivide(),true);
            TS_ASSERT_DELTA(p_model->GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(),101.0,1e-12);
            SimulationTime::Destroy();
            delete p_cell;
        }
    }
    
    // NB - to archive a cell cycle model it has to be archived via the cell that owns it.
    void TestArchiveWntCellCycleModel()
    {
        CancerParameters::Instance()->Reset();

        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "wnt_cell_cycle.arch";
        WntGradient::Instance()->SetConstantWntValueForTesting(1.0);

        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(16, 2);
                       
            WntCellCycleModel* p_cell_model = new WntCellCycleModel();
            
            TissueCell stem_cell(STEM, // type
                                       HEALTHY,//Mutation State
                                       0,  // generation
                                       p_cell_model);
            stem_cell.InitialiseCellCycleModel();  
            
            p_simulation_time->IncrementTimeOneStep();            
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(),false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(),true);

            stem_cell.GetCellCycleModel()->SetBirthTime(-1.0);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            TissueCell* const p_cell = &stem_cell;
            
            output_arch << static_cast<const SimulationTime&>(*p_simulation_time);
            output_arch << static_cast<const CancerParameters&>(*CancerParameters::Instance());
            output_arch << p_cell;
            SimulationTime::Destroy();
        }
        
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            CancerParameters *inst1 = CancerParameters::Instance();
            
            inst1->SetSDuration(101.0);
            
            TissueCell* p_cell;
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore from the archive
            input_arch >> *p_simulation_time;
            input_arch >> *inst1;
            input_arch >> p_cell;
            
            // Check
            
            AbstractCellCycleModel* p_cell_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());            
                 
            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),true);
            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_DELTA(p_cell_model->GetAge(),17.0,1e-12);
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
            SimulationTime::Destroy();
            delete p_cell;
        }

        WntGradient::Destroy();
    }   
    
    // NB - to archive a cell cycle model it has to be archived via the cell that owns it.
    void TestArchiveIngeWntSwatCellCycleModel()
    {
        CancerParameters::Instance()->Reset();

        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "inge_wnt_swat_cell_cycle.arch";
        WntGradient::Instance()->SetConstantWntValueForTesting(1.0);

        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(17, 2);
                       
            IngeWntSwatCellCycleModel* p_cell_model = new IngeWntSwatCellCycleModel(1);
            
            TissueCell stem_cell(STEM, // type
                                 HEALTHY,//Mutation State
                                 0,  // generation
                                 p_cell_model);
            stem_cell.InitialiseCellCycleModel();  
            
            p_simulation_time->IncrementTimeOneStep();            
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(),false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(),true);

            stem_cell.GetCellCycleModel()->SetBirthTime(-1.0);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            TissueCell* const p_cell = &stem_cell;
            
            output_arch << static_cast<const SimulationTime&>(*p_simulation_time);
            output_arch << static_cast<const CancerParameters&>(*CancerParameters::Instance());
            output_arch << p_cell;
            SimulationTime::Destroy();
        }
        
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            CancerParameters *inst1 = CancerParameters::Instance();
            
            inst1->SetSDuration(101.0);
            
            TissueCell* p_cell;
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore from the archive
            input_arch >> *p_simulation_time;
            input_arch >> *inst1;
            input_arch >> p_cell;
            
            // Check
            
            AbstractCellCycleModel* p_cell_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());            
                 
            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),true);
            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_DELTA(p_cell_model->GetAge(),18.0,1e-12);
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
            SimulationTime::Destroy();
            delete p_cell;
        }

        WntGradient::Destroy();
    }   
    
    // NB - to archive a cell cycle model it has to be archived via the cell that owns it.
    void TestArchiveSimpleWntCellCycleModel()
    {
        CancerParameters::Instance()->Reset();

        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "simple_wnt_cell_cycle.arch";
        WntGradient::Instance()->SetConstantWntValueForTesting(1.0);

        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(15, 30);
                       
            SimpleWntCellCycleModel* p_cell_model = new SimpleWntCellCycleModel();
            
            TissueCell stem_cell(STEM, // type
                                       HEALTHY,//Mutation State
                                       0,  // generation
                                       p_cell_model);
                                       
            p_cell_model->SetBirthTime(-1.0);
            // cell divides at t=11.0676, so this loop takes it up to time = 10, age = 11.
            for (unsigned i=0 ; i<20 ; i++)
            {
                p_simulation_time->IncrementTimeOneStep();            
                TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(),false);
            }
                        
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            TissueCell* const p_cell = &stem_cell;
            
            output_arch << static_cast<const SimulationTime&>(*p_simulation_time);
            output_arch << static_cast<const CancerParameters&>(*CancerParameters::Instance());
            output_arch << p_cell;
            
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(),true);
            
            
            SimulationTime::Destroy();
        }
        
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            CancerParameters *inst1 = CancerParameters::Instance();
            
            inst1->SetSDuration(101.0);
            
            TissueCell* p_cell;
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore from the archive
            input_arch >> *p_simulation_time;
            input_arch >> *inst1;
            input_arch >> p_cell;
            
            // Check            
            AbstractCellCycleModel* p_cell_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());            
                 
            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),true);
            
            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_cell_model->GetAge(), 11.5, 1e-12);
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(), 10.0, 1e-12);
            SimulationTime::Destroy();
            delete p_cell;
        }

        WntGradient::Destroy();
    }    
    
    void TestArchiveStochasticWntCellCycleModels()
    {
        CancerParameters::Instance()->Reset();
        RandomNumberGenerator::Instance()->Reseed(0);   // reset at beginning of this test.
        // In this case the first cycle time will be 5.971+9.0676 = 15.0386
        // note that the S-G2-M time is assigned when the cell finishes G1
        //(i.e. at time 5.971 here so the model has to be archived BEFORE that.
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "stochastic_wnt_cell_cycle.arch";
        WntGradient::Instance()->SetConstantWntValueForTesting(1.0);
        
        // Create an ouput archive
        {   // In this test the RandomNumberGenerator in existence 
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(16.0, 1000);
            
            StochasticWntCellCycleModel* p_stoc_model = new StochasticWntCellCycleModel();                    
                                           
            TissueCell stoc_cell(STEM, // type
                                       HEALTHY,//Mutation State
                                       0,  // generation
                                       p_stoc_model); 
            stoc_cell.InitialiseCellCycleModel();                                       
            
            WntCellCycleModel* p_wnt_model = new WntCellCycleModel();
            
            TissueCell wnt_cell(STEM, // type
                                      HEALTHY,//Mutation State
                                      0,  // generation
                                      p_wnt_model); 
            wnt_cell.InitialiseCellCycleModel();                                       
                                       
            p_simulation_time->IncrementTimeOneStep(); // 5.5
            
            while (p_simulation_time->GetDimensionalisedTime() < 4.0)
            {
                p_simulation_time->IncrementTimeOneStep();   
            }
            TS_ASSERT_EQUALS(stoc_cell.GetCellCycleModel()->ReadyToDivide(),false);
            TS_ASSERT_EQUALS(wnt_cell.GetCellCycleModel()->ReadyToDivide(),false);
            // When these are included here they pass - so are moved down into 
            // after load to see if they still pass.
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            TissueCell* const p_wnt_cell = &wnt_cell;
            TissueCell* const p_stoc_cell = &stoc_cell;
            
            output_arch << static_cast<const SimulationTime&>(*p_simulation_time);
            output_arch << static_cast<const CancerParameters&>(*CancerParameters::Instance());
            output_arch << p_stoc_cell;
            output_arch << p_wnt_cell;
            SimulationTime::Destroy();
        }
        
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(16.0, 2);
            
            CancerParameters *inst1 = CancerParameters::Instance();
            
            inst1->SetSDuration(101.0);
            
            TissueCell* p_stoc_cell; 
            TissueCell* p_wnt_cell; 
                      
            std::vector<double> cell_cycle_influence1;
            cell_cycle_influence1.push_back(1.0);
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            input_arch >> *p_simulation_time;
            input_arch >> *inst1;
            input_arch >> p_stoc_cell;
            input_arch >> p_wnt_cell;
            
            // Check - stochastic should divide at 15.03
            // Wnt should divide at 15.971
                        
            while (p_simulation_time->GetDimensionalisedTime() < 15.0)
            {
                p_simulation_time->IncrementTimeOneStep();   
            }
            TS_ASSERT_EQUALS(p_stoc_cell->GetCellCycleModel()->ReadyToDivide(),false);
            TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->ReadyToDivide(),false);
            
            while (p_simulation_time->GetDimensionalisedTime() < 15.5)
            {
                p_simulation_time->IncrementTimeOneStep();   
            }
            TS_ASSERT_EQUALS(p_stoc_cell->GetCellCycleModel()->ReadyToDivide(),true);// only for stochastic
            TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->ReadyToDivide(),false);
            
            while (p_simulation_time->GetDimensionalisedTime() < 16.0)
            {
                p_simulation_time->IncrementTimeOneStep();   
            }
            TS_ASSERT_EQUALS(p_stoc_cell->GetCellCycleModel()->ReadyToDivide(),true);
            TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->ReadyToDivide(),true);
            
            TS_ASSERT_DELTA(p_stoc_cell->GetCellCycleModel()->GetBirthTime(),0.0,1e-12);
            TS_ASSERT_DELTA(p_stoc_cell->GetCellCycleModel()->GetAge(),16.0,1e-12);
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
            
            delete p_stoc_cell;
            delete p_wnt_cell;
            SimulationTime::Destroy();
        }
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
    }    
    
    // NB - to archive a cell cycle model it has to be archived via the cell that owns it.
    void TestArchiveOxygenBasedCellCycleModels()
    {
        CancerParameters::Instance()->Reset();

        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "oxygen_based_cell_cycle.arch";

        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 2);
            
            OxygenBasedCellCycleModel* p_cell_model = new OxygenBasedCellCycleModel();
            
            TissueCell cell(HEPA_ONE, ALARCON_NORMAL, 0, p_cell_model);
            cell.InitialiseCellCycleModel();  
            cell.GetCellCycleModel()->SetBirthTime(-10.0);
            
            p_simulation_time->IncrementTimeOneStep();            
            TS_ASSERT_EQUALS(cell.GetCellCycleModel()->ReadyToDivide(),false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(cell.GetCellCycleModel()->ReadyToDivide(),true);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            TissueCell* const p_cell = &cell;
            
            output_arch << static_cast<const SimulationTime&>(*p_simulation_time);
            output_arch << static_cast<const CancerParameters&>(*CancerParameters::Instance());
            output_arch << p_cell;
            SimulationTime::Destroy();            
        }
        
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            CancerParameters *inst1 = CancerParameters::Instance();
            
            inst1->SetSDuration(101.0);
            
            TissueCell* p_cell;
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore from the archive
            input_arch >> *p_simulation_time;
            input_arch >> *inst1;
            input_arch >> p_cell;
            
            // Check            
            AbstractCellCycleModel* p_cell_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());            
                 
            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),true);
            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(),-10.0,1e-12);
            TS_ASSERT_DELTA(p_cell_model->GetAge(),20.0,1e-12);
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
            SimulationTime::Destroy();
            delete p_cell;
        }

        CellwiseData<2>::Destroy();
    }    
    
};

#endif /*TESTCELLCYCLEMODELS_HPP_*/
