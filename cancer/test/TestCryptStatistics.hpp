#ifndef TESTCRYPTSTATISTICS_HPP_
#define TESTCRYPTSTATISTICS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "Tissue.cpp"
#include "CryptStatistics.hpp"
#include "CryptSimulation2d.hpp"
#include "OutputFileHandler.hpp"
#include "TissueCell.hpp"
#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"
#include "SloughingCellKiller.hpp"
#include "CellsGenerator.hpp"


class TestCryptStatistics : public CxxTest::TestSuite
{
public:
    void TestGetSection() throw (Exception)
    {        
        CancerParameters::Instance()->Reset();

        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 0;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);// true = mature cells

        Tissue<2> crypt(*p_mesh, cells);               
        crypt.SetGhostNodes(ghost_node_indices);

        CryptStatistics crypt_statistics(crypt);
        
        std::vector< TissueCell* > test_section=crypt_statistics.GetCryptSection(0.5,1.5,sqrt(3));
        
        //Test the cells are correct
        TS_ASSERT_EQUALS(test_section.size(), 6u);

        unsigned expected_indices[6]={0,1,3,4,7,8};

        for(unsigned i=0; i<test_section.size(); i++)
        {
            TissueCell* cell = test_section[i];
            TS_ASSERT_EQUALS( crypt.GetNodeCorrespondingToCell(*cell)->GetIndex(), 
                              expected_indices[i]);
        }

        
        // Test that we get a valid section when the x-values are the same
        std::vector< TissueCell* > test_section_vertical=crypt_statistics.GetCryptSection(0.5,0.5,sqrt(3));
        
        //Test the cells are correct
        TS_ASSERT_EQUALS(test_section_vertical.size(), 5u);

        unsigned expected_indices_vertical[6]={0,1,3,6,7};

        for(unsigned i=0; i<test_section_vertical.size(); i++)
        {
            TissueCell* cell = test_section_vertical[i];
            TS_ASSERT_EQUALS( crypt.GetNodeCorrespondingToCell(*cell)->GetIndex(), 
                              expected_indices_vertical[i]);
        }
        
        std::vector< TissueCell* > test_section_periodic=crypt_statistics.GetCryptSectionPeriodic(0.5,2.5,sqrt(3));
        
        //Test the cells are correct
        TS_ASSERT_EQUALS(test_section_periodic.size(), 6u);
        
        unsigned expected_indices_periodic[6]={0,1,3,5,6,8};
        
        for(unsigned i=0; i<test_section_periodic.size(); i++)
        {
            TissueCell* cell = test_section_periodic[i];
            TS_ASSERT_EQUALS( crypt.GetNodeCorrespondingToCell(*cell)->GetIndex(), 
                              expected_indices_periodic[i]);
        }
        
        
        std::vector< TissueCell* > test_section_periodic_2=crypt_statistics.GetCryptSectionPeriodic(2.5,0.5,sqrt(3));
        
        //Test the cells are correct
        TS_ASSERT_EQUALS(test_section_periodic_2.size(), 6u);
        
        unsigned expected_indices_periodic_2[6]={0,2,3,5,6,7};
        
        for(unsigned i=0; i<test_section_periodic_2.size(); i++)
        {
            TissueCell* cell = test_section_periodic_2[i];
            TS_ASSERT_EQUALS( crypt.GetNodeCorrespondingToCell(*cell)->GetIndex(), 
                              expected_indices_periodic_2[i]);
        }
        
        // Test an overwritten method
        std::vector< TissueCell* > test_section_periodic_3 = crypt_statistics.GetCryptSectionPeriodic();
        //Test the cells are correct
        TS_ASSERT_EQUALS(test_section_periodic_3.size(), 3u);
        unsigned expected_indices_periodic_3[6]={1,4,8};
        
        for(unsigned i=0; i<test_section_periodic_3.size(); i++)
        {
            TissueCell* cell = test_section_periodic_3[i];
            TS_ASSERT_EQUALS( crypt.GetNodeCorrespondingToCell(*cell)->GetIndex(), 
                              expected_indices_periodic_3[i]);
        }
        
        SimulationTime::Destroy();
    }
    
    
    void TestMakeMeinekeGraphs() throw (Exception)
    {        
        CancerParameters* p_params = CancerParameters::Instance();
        std::string output_directory = "MakeMeinekeGraphs";
        
        //double end_of_simulation = 1.0; // hours
        
        unsigned cells_across = 13;
        unsigned cells_up = 25;
        double crypt_width = 12.1;
        unsigned thickness_of_ghost_layer = 3;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, STOCHASTIC, true);
              
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
                
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory(output_directory);
        double time_of_each_run = simulator.GetDt(); // for each run
        
        // Set simulation to output cell types
        simulator.SetOutputCellTypes(true);
                
        // Set length of simulation here
        simulator.SetEndTime(time_of_each_run);
        AbstractCellKiller<2>* p_cell_killer = new SloughingCellKiller(&simulator.rGetTissue(),0.01);
        simulator.AddCellKiller(p_cell_killer);
        
        // UNUSUAL SET UP HERE /////////////////////////////////////
        
        p_params->SetDampingConstantNormal(1.0);    // normally 1

        // Do not give mutant cells any different movement properties to normal ones
        p_params->SetDampingConstantMutant(p_params->GetDampingConstantNormal());
        
        p_params->SetSpringStiffness(30.0); //normally 15.0;
        // 0.3/30 = 0.01 (i.e. Meineke's values)
        
        simulator.UseJiggledBottomCells();
        
        // END OF UNUSUAL SET UP! //////////////////////////////////
        
        
        // TEST CryptStatistics::GetCryptSectionPeriodic by labelling a column of cells...
        CryptStatistics crypt_statistics(crypt);
        std::vector< TissueCell* > test_section=crypt_statistics.GetCryptSectionPeriodic(8.0,8.0);
        
        for (unsigned i=0; i<test_section.size() ; i++)
        {
            test_section[i]->SetMutationState(LABELLED);
        }
        
        simulator.Solve();
        
        simulator.Save();

        // ... and checking visualization of labelled cells against previous run
        OutputFileHandler handler("MakeMeinekeGraphs",false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/vis_results/results.viznodes";
        TS_ASSERT_EQUALS(system(("cmp " + results_file + " cancer/test/data/MakeMeinekeGraphs/results_from_time_0/vis_results/results.viznodes").c_str()), 0);

        // TEST crypt_statistics::LabelSPhaseCells
        
        // First remove labels.
        for (Tissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            (*cell_iter).SetMutationState(HEALTHY);
        } 
        
        crypt_statistics.LabelSPhaseCells();

        // Iterate over cells checking for correct labels
        unsigned counter =0;
        for (Tissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            bool is_labelled= (*cell_iter).GetMutationState() == LABELLED;
            
            bool in_s_phase = (*cell_iter).GetCellCycleModel()->GetCurrentCellCyclePhase()== S;
            
            TS_ASSERT_EQUALS(is_labelled, in_s_phase);
            
            if(in_s_phase)
            {
                counter++;
            }
                        
        } 
        TS_ASSERT_EQUALS(counter,15u);       
        
        simulator.SetEndTime(2*time_of_each_run);
        simulator.Solve();
        
        // TEST CryptStatistics::GetWhetherCryptSectionCellsAreLabelled
        
        //std::vector<bool> labelled = crypt_statistics.GetWhetherCryptSectionCellsAreLabelled(2.5,0.5,sqrt(3));
        
        // set cells which are not in the crypt section to be in state APC_ONE_HIT, so that we can
        // see the section
        test_section=crypt_statistics.GetCryptSectionPeriodic(8.0,8.0);
                for (Tissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            bool in_section=false;
            for (unsigned vector_index=0; vector_index<test_section.size(); vector_index++)
            {
                if (test_section[vector_index]==&(*cell_iter))
                {
                    in_section=true;
                }
            }
            if (!in_section)
            {
                (*cell_iter).SetMutationState(APC_ONE_HIT);
            }
                        
        } 
        simulator.SetEndTime(3*time_of_each_run);
        simulator.Solve();
        
        std::vector<bool> labelled = crypt_statistics.GetWhetherCryptSectionCellsAreLabelled(8.0,8.0);
        
        for (unsigned vector_index=0; vector_index<labelled.size(); vector_index++)
        {
            if (vector_index==3u)
            {
                TS_ASSERT(labelled[vector_index]);
            }
            else
            {
                TS_ASSERT(!labelled[vector_index]);
            }
        }
        
        delete p_cell_killer;
        delete p_params;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
    }
    
};

#endif /*TESTCRYPTSTATISTICS_HPP_*/
