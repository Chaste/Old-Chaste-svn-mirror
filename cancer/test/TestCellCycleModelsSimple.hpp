#ifndef TESTCELLCYCLEMODELSSIMPLE_HPP_
#define TESTCELLCYCLEMODELSSIMPLE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "ConformingTetrahedralMesh.cpp"
#include "CellsGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "WntGradient.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestCellCycleModelsSimple : public AbstractCancerTestSuite
{               
public:

    void TestFixedCellCycleModel() throw(Exception)
    {   
        CancerParameters *p_params = CancerParameters::Instance();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        
        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            2.0*(p_params->GetStemCellG1Duration()
                  +p_params->GetSG2MDuration()     ), num_steps);
        
        TS_ASSERT_THROWS_NOTHING(FixedCellCycleModel model3);
        
        FixedCellCycleModel* p_stem_model = new FixedCellCycleModel;
        TS_ASSERT(!p_stem_model->UsesBetaCat());
        TissueCell stem_cell(STEM, HEALTHY, p_stem_model);
        stem_cell.InitialiseCellCycleModel();
        
        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(),M_PHASE);
        TS_ASSERT_EQUALS(p_stem_model->GetGeneration(), 0u);
        
        TS_ASSERT_EQUALS(stem_cell.GetCellType(),STEM);
        
        FixedCellCycleModel* p_transit_model = new FixedCellCycleModel;
        TissueCell transit_cell(TRANSIT, HEALTHY, p_transit_model);
        transit_cell.InitialiseCellCycleModel();
                           
        TS_ASSERT_EQUALS(transit_cell.GetCellType(),TRANSIT);
        TS_ASSERT_EQUALS(p_transit_model->GetGeneration(), 0u);
        
        FixedCellCycleModel* p_diff_model = new FixedCellCycleModel;
        TissueCell diff_cell(DIFFERENTIATED, HEALTHY, p_diff_model);
        diff_cell.InitialiseCellCycleModel();
                           
        TS_ASSERT_EQUALS(diff_cell.GetCellType(),DIFFERENTIATED);
        TS_ASSERT_EQUALS(p_diff_model->GetGeneration(), 0u);
                        
        for (unsigned i = 0 ; i< num_steps ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, p_params->GetStemCellG1Duration());
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, p_params->GetTransitCellG1Duration());
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 100);           
        }
        
        TS_ASSERT_DELTA(p_stem_model->GetAge(), p_simulation_time->GetDimensionalisedTime(), 1e-9);
        TS_ASSERT_DELTA(p_transit_model->GetAge(), p_simulation_time->GetDimensionalisedTime(), 1e-9);
        TS_ASSERT_DELTA(p_diff_model->GetAge(), p_simulation_time->GetDimensionalisedTime(), 1e-9);
        
        double hepa_one_cell_birth_time = p_simulation_time->GetDimensionalisedTime();
        
        p_params->SetHepaOneParameters();        
        FixedCellCycleModel* p_hepa_one_model = new FixedCellCycleModel;
        TissueCell hepa_one_cell(STEM, HEALTHY, p_hepa_one_model);
        hepa_one_cell.InitialiseCellCycleModel();
                
        for (unsigned i = 0 ; i< num_steps ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, p_params->GetHepaOneCellG1Duration());
        }        
        
        TS_ASSERT_DELTA(p_hepa_one_model->GetAge() + hepa_one_cell_birth_time, p_simulation_time->GetDimensionalisedTime(), 1e-9);
    }
    
    
    void TestStochasticCellCycleModel(void) throw(Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();

        SimulationTime* p_simulation_time = SimulationTime::Instance();        
        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            2.0*(p_params->GetStemCellG1Duration() + p_params->GetSG2MDuration()), num_steps);
        
        TS_ASSERT_THROWS_NOTHING(StochasticCellCycleModel cell_model3);
        
        StochasticCellCycleModel* p_stem_model = new StochasticCellCycleModel;
        StochasticCellCycleModel* p_transit_model = new StochasticCellCycleModel;
        StochasticCellCycleModel* p_diff_model = new StochasticCellCycleModel;
                
        TissueCell stem_cell(STEM, HEALTHY,  p_stem_model);
        stem_cell.InitialiseCellCycleModel();
                        
        TissueCell transit_cell(TRANSIT, HEALTHY, p_transit_model);
        transit_cell.InitialiseCellCycleModel();
        
        TissueCell diff_cell(DIFFERENTIATED, HEALTHY, p_diff_model); 
        diff_cell.InitialiseCellCycleModel();
                                  
        for (unsigned i = 0 ; i< num_steps ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            
            // The numbers for the G1 durations below are taken from the first three 
            // random numbers generated    
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, 4.36075);
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, 1.78877);            
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132);  // any old number            
        }
        
        p_params->SetHepaOneParameters();                
        StochasticCellCycleModel* p_hepa_one_model = new StochasticCellCycleModel;
        TissueCell hepa_one_cell(STEM, HEALTHY, p_hepa_one_model);
        hepa_one_cell.InitialiseCellCycleModel();
                
        for (unsigned i = 0 ; i< num_steps ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, 4.1324);
        }        
    }
    
    
    void TestSimpleWntCellCycleModel() throw(Exception)
    {
        CancerParameters* p_params = CancerParameters::Instance();
        
        // Set up the simulation time
        SimulationTime *p_simulation_time = SimulationTime::Instance();   
        
        double end_time = 60.0;         
        unsigned num_timesteps = 1000*(unsigned)end_time;     
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);
        
        // Set up the Wnt gradient
        double wnt_level = 1.0;
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
                
        SimpleWntCellCycleModel* p_cycle_model = new SimpleWntCellCycleModel;
        TissueCell cell(STEM, HEALTHY, p_cycle_model);
        cell.InitialiseCellCycleModel();
        
        for (unsigned i = 0 ; i< num_timesteps/3 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            
            // The number for the G1 duration is taken from 
            // the first random number generated    
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, 1.0676);
        }
        // Stem cell should have been changed into a transit cell by wnt cell cycle model
        TS_ASSERT_EQUALS(cell.GetCellType(), TRANSIT);
        
        // Divide the cell
        TS_ASSERT_EQUALS(cell.ReadyToDivide(), true);
        TissueCell cell2 = cell.Divide();
        cell.SetMutationState(LABELLED);
        
        SimpleWntCellCycleModel *p_cycle_model2 = static_cast <SimpleWntCellCycleModel*> (cell2.GetCellCycleModel());        
        
        // Now reduce the Wnt gradient
        wnt_level = 0.7;
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
        
        double division_time = SimulationTime::Instance()->GetDimensionalisedTime();
        
        // The numbers for the G1 durations are taken from 
        // the first two random numbers generated
        double new_g1_duration = 3.16316;
        double new_g1_duration2 = 1.2712;
        for (unsigned i = 0 ; i< num_timesteps/3 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();            
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, new_g1_duration);
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model2, new_g1_duration2);
        }
        
        TS_ASSERT_EQUALS(cell.GetCellType(), TRANSIT);
        TS_ASSERT_EQUALS(cell2.GetCellType(), TRANSIT);
                                    
        p_cycle_model->ResetForDivision();
        p_cycle_model2->ResetForDivision();
                 
        division_time = SimulationTime::Instance()->GetDimensionalisedTime();                
                
        // Now reduce the Wnt gradient so only beta-cat or APC2 hit cells divide.
        wnt_level = 0.15;
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
        
        cell.SetMutationState(APC_ONE_HIT);
        cell2.SetMutationState(BETA_CATENIN_ONE_HIT);
        
        TS_ASSERT(!p_cycle_model->ReadyToDivide());
        TS_ASSERT(!p_cycle_model2->ReadyToDivide());
        
        // Coverage...
        cell.SetMutationState(APC_TWO_HIT);
        TS_ASSERT(!p_cycle_model->ReadyToDivide());
        cell.SetMutationState(APC_ONE_HIT);
        TS_ASSERT(!p_cycle_model->ReadyToDivide());
        
        // The numbers for the G1 durations are taken from 
        // the first two random numbers generated
        new_g1_duration = 1.22037;
        new_g1_duration2 = 0.74699;
        
        for (unsigned i = 0 ; i< num_timesteps/3 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();            
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, new_g1_duration);
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model2, new_g1_duration2);
        }
                
        TS_ASSERT_EQUALS(cell.GetCellType(), DIFFERENTIATED);
        TS_ASSERT_EQUALS(cell2.GetCellType(), TRANSIT);
        
        // For coverage...
        SimpleWntCellCycleModel* p_cycle_model1 = new SimpleWntCellCycleModel;
        TissueCell cell1(DIFFERENTIATED, HEALTHY, p_cycle_model1);
        cell1.InitialiseCellCycleModel();
        // ...end of coverage
        
        
        /*
         * Test the case of a radial Wnt gradient
         */ 
        
        p_params->Reset(); 
        RandomNumberGenerator::Instance()->Reseed(0);
        
        // Set up the Wnt gradient
        wnt_level = p_params->GetWntStemThreshold() + 0.01;
        WntGradient::Destroy();
        WntGradient::Instance()->SetType(RADIAL);
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
                
        // Set up a cell cycle model and cell        
        SimpleWntCellCycleModel* p_cycle_model4 = new SimpleWntCellCycleModel;
        TissueCell cell4(STEM, HEALTHY,  p_cycle_model4);
        cell4.InitialiseCellCycleModel();
        
        // Test the GetCurrentCellCyclePhase() and ReadyToDivide() methods
        double first_g1_duration = 1.0676;
        for (unsigned i = 0 ; i< num_timesteps/3 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            
            // The number for the G1 duration is taken from 
            // the first random number generated    
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model4, first_g1_duration);
        }

        // We should still have a stem cell since the WntGradient exceeds mRadialWntThreshold
        TS_ASSERT_EQUALS(cell4.GetCellType(), STEM);
        
        // Divide the cell
        TS_ASSERT_EQUALS(cell4.ReadyToDivide(), true);
        TS_ASSERT_EQUALS(cell4.GetCellType(), STEM);
        TissueCell cell5 = cell4.Divide();
        TS_ASSERT_EQUALS(cell4.GetCellType(), STEM);
        TS_ASSERT_EQUALS(cell5.GetCellType(), TRANSIT);
        cell2.SetMutationState(LABELLED);
            
        // Now reduce the Wnt gradient
        wnt_level = p_params->GetWntStemThreshold() - 0.01;
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
              
        // The numbers for the G1 durations are taken from 
        // the first two random numbers generated
        new_g1_duration = 3.16316;
        for (unsigned i = 0 ; i< num_timesteps/3 ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();            
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model4, new_g1_duration);
        }
        
        TS_ASSERT_DELTA(WntGradient::Instance()->GetWntLevel(&cell4), wnt_level, 1e-12);
        TS_ASSERT_EQUALS(cell4.GetCellType(), TRANSIT);
        TS_ASSERT_EQUALS(cell5.GetCellType(), TRANSIT);
        
        WntGradient::Destroy();
    }
    
    
    void TestArchiveFixedCellCycleModels() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "fixed_cell_cycle.arch";
        
        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 4);
            FixedCellCycleModel* p_model = new FixedCellCycleModel;
            
            TissueCell cell(TRANSIT, HEALTHY, p_model);
            cell.InitialiseCellCycleModel();
            
            p_simulation_time->IncrementTimeOneStep();
            p_simulation_time->IncrementTimeOneStep();
            
            p_model->SetBirthTime(-1.0);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            // Update cell phase
            p_model->ReadyToDivide();
            
            TissueCell* const p_cell = &cell;
            
            output_arch << p_cell;
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(),G_ONE_PHASE);
            
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
            
            // Restore from the archive
            input_arch >> p_cell;
            
            AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();
            
            // Check private data has been restored correctly.
            TS_ASSERT_DELTA(p_model->GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(),2.5,1e-12);
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(),G_ONE_PHASE);
            
            delete p_cell;
        }
    }
    
        
    void TestArchiveStochasticCellCycleModels()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "stoch_cell_cycle.arch";
        
        double random_number_test = 0;
        
        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);
            
            StochasticCellCycleModel* p_model = new StochasticCellCycleModel;
            
            TissueCell cell(TRANSIT,  HEALTHY, p_model);
            cell.InitialiseCellCycleModel();
            cell.SetBirthTime(-1.1);
            p_simulation_time->IncrementTimeOneStep();
            p_simulation_time->IncrementTimeOneStep();
            
            cell.ReadyToDivide(); // updates phases
            
            TissueCell* const p_cell = &cell;
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            TS_ASSERT_DELTA(CancerParameters::Instance()->GetSDuration(),5.0,1e-12);
            
            output_arch << p_cell;
            
            TS_ASSERT_DELTA(p_model->GetBirthTime(),-1.1,1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(),2.1,1e-12);
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(),G_ONE_PHASE);
            
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            random_number_test = p_gen->ranf();
            SimulationTime::Destroy();
        }
        
        {           
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            p_gen->Reseed(128);
            
            TissueCell* p_cell;
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            CancerParameters *inst1 = CancerParameters::Instance();
            
            inst1->SetSDuration(101.0);
            
            // Restore from the archive
            input_arch >> p_cell;
            
            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(),random_number_test,1e-7);
            
            AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();
            
            // Check
            TS_ASSERT_DELTA(p_model->GetBirthTime(),-1.1,1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(),2.1,1e-12);
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(),G_ONE_PHASE);
            
            TS_ASSERT_DELTA(inst1->GetSDuration(),5.0,1e-12);
            
            delete p_cell;
        }
    }
    
    
    void TestArchiveSimpleOxygenBasedCycleModels() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "oxygen_based_cell_cycle.arch";

        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);
                        
        // Create an ouput archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);
            
            SimpleOxygenBasedCellCycleModel model;
            
            p_simulation_time->IncrementTimeOneStep();
            
            model.SetBirthTime(-1.0);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);            
                        
            output_arch << static_cast<const SimpleOxygenBasedCellCycleModel&>(model);
            
            SimulationTime::Destroy();     
        }
        
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            SimpleOxygenBasedCellCycleModel model;
            model.SetBirthTime(-2.0);            
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // Restore from the archive
            input_arch >> model;
            
            // Check that archiiving worked correctly
            TS_ASSERT_DELTA(model.GetBirthTime(),-1.0,1e-12);
            TS_ASSERT_DELTA(model.GetAge(),1.5,1e-12);
            TS_ASSERT_EQUALS(model.GetCurrentCellCyclePhase(),M_PHASE);     
        }
    }
       
    
    void TestArchiveSimpleWntCellCycleModel()
    {
        CancerParameters* p_params = CancerParameters::Instance();
        
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "simple_wnt_cell_cycle.arch";
        
        // Set up the Wnt gradient
        double wnt_level = 1.0;
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
        
        double random_number_test = 0;
        
        // Create an ouput archive
        {
            // Set up the simulation time
            SimulationTime *p_simulation_time = SimulationTime::Instance();   
            
            // The number for the G1 duration is taken from 
            // the first random number generated    
            double g1_duration = 1.0676;
            
            double end_time = g1_duration + p_params->GetSG2MDuration() + 5.0; 
                    
            unsigned num_timesteps = 50;   
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);
                           
            SimpleWntCellCycleModel* p_cell_model = new SimpleWntCellCycleModel();
            
            p_cell_model->SetBirthTime(-1.0);
            
            TissueCell stem_cell(STEM, HEALTHY, p_cell_model);
            stem_cell.InitialiseCellCycleModel();
                                       
            while (p_cell_model->GetAge() < 
                g1_duration + p_params->GetSG2MDuration() 
                - p_simulation_time->GetTimeStep()) // minus one to match birth time.
            {
                p_simulation_time->IncrementTimeOneStep();  
                CheckReadyToDivideAndPhaseIsUpdated(p_cell_model, g1_duration);
            }
            
            // Wnt should change this to a transit cell
            TS_ASSERT_EQUALS(stem_cell.GetCellType(), TRANSIT);
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->ReadyToDivide(), false);
            TS_ASSERT_EQUALS(stem_cell.GetCellCycleModel()->GetCurrentCellCyclePhase(), G_TWO_PHASE);                             
                       
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            TissueCell* const p_cell = &stem_cell;
            
            output_arch << p_cell;
            
            TS_ASSERT_EQUALS(p_cell->ReadyToDivide(),false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_cell->ReadyToDivide(),true);
            
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT(stem_cell.GetCellCycleModel()->ReadyToDivide());     
                   
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            random_number_test = p_gen->ranf();
            SimulationTime::Destroy();
        }        
        {             
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            CancerParameters *inst1 = CancerParameters::Instance();
            
            inst1->SetSDuration(101.0);
            
            TissueCell* p_cell;
                        
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            p_gen->Reseed(36);
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // Restore from the archive
            input_arch >> p_cell;
            
            // Check            
            AbstractCellCycleModel* p_cell_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());            
                 
            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),true);
            
            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(), 10.0, 1e-12);
            
            TS_ASSERT_DELTA(p_gen->ranf(),random_number_test,1e-7);
            
            SimulationTime::Destroy();
            delete p_cell;
        }
        
        /*
         * Test the case of a radial Wnt gradient
         */ 

        p_params->Reset();
        RandomNumberGenerator::Instance()->Reseed(0);
        
        OutputFileHandler handler2("archive", false);
        archive_filename = handler2.GetOutputDirectoryFullPath() + "crypt_projection_cell_cycle.arch";
        
        // Set up the Wnt gradient
        wnt_level = p_params->GetWntStemThreshold() - 0.01;
        WntGradient::Destroy();
        WntGradient::Instance()->SetConstantWntValueForTesting(wnt_level);
        
        random_number_test = 0;
        
        // Create an ouput archive
        {
            // Set up the simulation time
            SimulationTime *p_simulation_time = SimulationTime::Instance();   
            
            // The number for the G1 duration is taken from 
            // the first random number generated    
            double g1_duration = 1.0676;
            
            double end_time = g1_duration + p_params->GetSG2MDuration() + 5.0; 
                    
            unsigned num_timesteps = 50;   
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);
                           
            SimpleWntCellCycleModel* p_cell_model = new SimpleWntCellCycleModel();
            
            p_cell_model->SetBirthTime(-1.0);
            
            TissueCell cell(STEM, HEALTHY, p_cell_model);
            cell.InitialiseCellCycleModel();
            
            // Run to division age minus one time step to match birth time                                   
            while (p_cell_model->GetAge() < g1_duration + p_params->GetSG2MDuration() 
                                            - p_simulation_time->GetTimeStep()) 
            {                     
                p_simulation_time->IncrementTimeOneStep();  
                CheckReadyToDivideAndPhaseIsUpdated(p_cell_model, g1_duration);
            }
            
            // Wnt should change this to a transit cell
            TS_ASSERT_EQUALS(cell.GetCellType(), TRANSIT);
            TS_ASSERT_EQUALS(cell.GetCellCycleModel()->ReadyToDivide(), false);
            TS_ASSERT_EQUALS(cell.GetCellCycleModel()->GetCurrentCellCyclePhase(), G_TWO_PHASE);                             
                       
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            TissueCell* const p_cell = &cell;
            
            output_arch << p_cell;
            
            TS_ASSERT_EQUALS(p_cell->ReadyToDivide(),false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_cell->ReadyToDivide(),true);
            
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT(cell.GetCellCycleModel()->ReadyToDivide());     
                   
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            random_number_test = p_gen->ranf();
                       
            RandomNumberGenerator::Destroy();
            SimulationTime::Destroy();
        }        
        {             
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            
            CancerParameters *inst1 = CancerParameters::Instance();
            
            inst1->SetSDuration(101.0);
            
            TissueCell* p_cell;
                        
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            p_gen->Reseed(36);
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // Restore from the archive
            input_arch >> p_cell;
            
            // Check            
            AbstractCellCycleModel* p_cell_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());            
            
            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),false);
            TS_ASSERT_EQUALS(p_cell->GetCellType(), TRANSIT);     
            
            p_simulation_time->IncrementTimeOneStep();
            
            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),true);
            TS_ASSERT_EQUALS(p_cell->GetCellType(), TRANSIT);     
            
            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(), 10.0, 1e-12);
            
            TS_ASSERT_DELTA(p_gen->ranf(),random_number_test,1e-7);            
            
            delete p_cell;
        }
        
        WntGradient::Destroy();
    }   
    
    
    void TestSimpleOxygenBasedCellCycleModel(void) throw(Exception)
    {           
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->SetHepaOneParameters();
        
        // Check that mCurrentHypoxiaOnsetTime and mCurrentHypoxicDuration are 
        // updated correctly
        SimulationTime *p_simulation_time = SimulationTime::Instance();  
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);
        
        SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
        TissueCell cell(STEM, HEALTHY, p_model);
        cell.InitialiseCellCycleModel();
        
        // Set up constant oxygen_concentration     
        std::vector<double> low_oxygen_concentration;
        std::vector<double> high_oxygen_concentration;
        low_oxygen_concentration.push_back(0.0);
        high_oxygen_concentration.push_back(1.0);
        
        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);
        
        p_model->ReadyToDivide();                
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);
        
        p_simulation_time->IncrementTimeOneStep(); // t=1.0
        p_model->ReadyToDivide(); 
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);
        
        CellwiseData<2>::Instance()->SetConstantDataForTesting(high_oxygen_concentration);
        
        p_simulation_time->IncrementTimeOneStep(); // t=2.0
        p_model->ReadyToDivide(); 
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 2.0, 1e-12);
        
        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);
        p_simulation_time->IncrementTimeOneStep(); // t=3.0
        p_model->ReadyToDivide(); 
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 2.0, 1e-12);
        
        // Set up SimulationTime       
        SimulationTime::Destroy();  
        p_simulation_time = SimulationTime::Instance();   
                
        unsigned num_steps = 100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            4.0*(p_params->GetHepaOneCellG1Duration()
                  +p_params->GetSG2MDuration()     ), num_steps);
        
        // Set up constant oxygen_concentration     
        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        TS_ASSERT_THROWS_NOTHING(SimpleOxygenBasedCellCycleModel model);
        
        // Create cell cycle model
        SimpleOxygenBasedCellCycleModel* p_hepa_one_model = new SimpleOxygenBasedCellCycleModel();
        SimpleOxygenBasedCellCycleModel* p_diff_model = new SimpleOxygenBasedCellCycleModel();
        
        // Create cell 
        TissueCell hepa_one_cell(STEM, HEALTHY, p_hepa_one_model);
        hepa_one_cell.InitialiseCellCycleModel();
        
        TissueCell diff_cell(DIFFERENTIATED, HEALTHY, p_diff_model);
        diff_cell.InitialiseCellCycleModel();
        
        // Check that the cell cycle phase and ready to divide
        // are updated correctly        
        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(),false);        
        TS_ASSERT_EQUALS(p_hepa_one_model->GetCurrentCellCyclePhase(),M_PHASE);
        
        TS_ASSERT_EQUALS(p_diff_model->ReadyToDivide(),false);
        TS_ASSERT_EQUALS(p_diff_model->GetCurrentCellCyclePhase(),G_ZERO_PHASE);
                
        for (unsigned i = 0 ; i< num_steps ; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            
            // Note that we need to pass in the updated G1 duration            
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, p_hepa_one_model->GetG1Duration());         
        }
        
        TS_ASSERT_DELTA(p_hepa_one_model->GetAge(), p_simulation_time->GetDimensionalisedTime(), 1e-9);
        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(),true);  

        // Check that cell division correctly resets the cell cycle phase
        SimpleOxygenBasedCellCycleModel *p_hepa_one_model2 = static_cast <SimpleOxygenBasedCellCycleModel*> (p_hepa_one_model->CreateCellCycleModel());
        
        TissueCell hepa_one_cell2(STEM, HEALTHY, p_hepa_one_model2);
        TS_ASSERT_EQUALS(p_hepa_one_model2->ReadyToDivide(), false);        
        TS_ASSERT_EQUALS(p_hepa_one_model2->GetCurrentCellCyclePhase(), M_PHASE);
        
        TS_ASSERT_THROWS_NOTHING(p_hepa_one_model->ResetForDivision());     
        
        // Set up SimulationTime         
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();                  
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0*CancerParameters::Instance()->GetCriticalHypoxicDuration(), num_steps);
        
        // Create a cell with a simple oxygen-based cell cycle model
        SimpleOxygenBasedCellCycleModel* p_cell_model = new SimpleOxygenBasedCellCycleModel();
        TissueCell necrotic_cell(STEM, HEALTHY, p_cell_model);
        
        // Set up constant oxygen_concentration     
        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);
                          
        // Force the cell to be necrotic
        for (unsigned i = 0 ; i< num_steps ; i++)
        {
            TS_ASSERT(necrotic_cell.GetCellType()!=NECROTIC || 
                      p_simulation_time->GetDimensionalisedTime() >= CancerParameters::Instance()->GetCriticalHypoxicDuration());
            p_simulation_time->IncrementTimeOneStep();
            
            // Note that we need to pass in the updated G1 duration            
            necrotic_cell.ReadyToDivide();
        }
        
        // Test that the cell type is updated to be NECROTIC        
        TS_ASSERT(necrotic_cell.GetCellType()==NECROTIC);          
        TS_ASSERT_EQUALS(p_cell_model->GetCurrentHypoxicDuration(), 2.04);
                  
        CellwiseData<2>::Destroy();
    }    
};

#endif /*TESTCELLCYCLEMODELSSIMPLE_HPP_*/
