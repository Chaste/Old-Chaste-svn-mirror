#ifndef TESTCRYPTSIMULATION2DNIGHTLY_HPP_
#define TESTCRYPTSIMULATION2DNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "OutputFileHandler.hpp"
#include "TissueCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "WntGradient.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"
#include "Tissue.cpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "CellsGenerator.hpp"

// Simple cell killer which just kills a single cell.
// The constructor takes in a number, and the killer
// will kill the number-th cell reached using the iterator
// (or the last cell, if n>num_cells)
//
// For testing purposes
class SingleCellCellKiller : public AbstractCellKiller<2>
{
private :
    unsigned mNumber;

public :
    SingleCellCellKiller(Tissue<2>* pTissue, unsigned number)
        : AbstractCellKiller<2>(pTissue),
          mNumber(number)
    {
    }
    
    virtual void TestAndLabelCellsForApoptosisOrDeath()
    {
        if(mpTissue->GetNumRealCells()==0)
        {
            return;
        }
        
        Tissue<2>::Iterator cell_iter = mpTissue->Begin();
       
        for(unsigned i=0; ( (i<mNumber) && (cell_iter!=mpTissue->End()) ); i++)
        {
            ++cell_iter;
        }
        
        cell_iter->Kill();
    }
};


class TestCryptSimulation2dNightly : public CxxTest::TestSuite
{
    void CheckAgainstPreviousRun(std::string resultDirectory, std::string resultSet, unsigned maxCells, unsigned maxElements)
    {
        std::cout << "Comparing " << resultDirectory << std::endl << std::flush;
        
        ColumnDataReader computed_node_results = ColumnDataReader(resultDirectory+"/"+resultSet+"/tab_results",
                                                                  "tabulated_node_results",
                                                                  true);
                                                                  
        ColumnDataReader expected_node_results = ColumnDataReader("cancer/test/data/" + resultDirectory+"Results",
                                                                  "tabulated_node_results",
                                                                  false);
        ColumnDataReader computed_element_results = ColumnDataReader(resultDirectory+"/"+resultSet+"/tab_results",
                                                    "tabulated_element_results",
                                                    true);
                                                    
        ColumnDataReader expected_element_results = ColumnDataReader("cancer/test/data/" + resultDirectory+"Results",
                                                    "tabulated_element_results",
                                                    false);
        
        for (unsigned cell=0; cell<maxCells; cell++)
        {
            std::stringstream cell_type_var_name;
            std::stringstream cell_x_position_var_name;
            std::stringstream cell_y_position_var_name;
            cell_type_var_name << "cell_type_" << cell;
            cell_x_position_var_name << "cell_x_position_" << cell;
            cell_y_position_var_name << "cell_y_position_" << cell;
            
            // Vector of Cell Types
            std::vector<double> expected_cell_types = expected_node_results.GetValues(cell_type_var_name.str());
            std::vector<double> computed_cell_types = computed_node_results.GetValues(cell_type_var_name.str());
            
            //Vector of Cell Positions
            std::vector<double> expected_cell_x_positions = expected_node_results.GetValues(cell_x_position_var_name.str());
            std::vector<double> computed_cell_x_positions = computed_node_results.GetValues(cell_x_position_var_name.str());
            
            std::vector<double> expected_cell_y_positions = expected_node_results.GetValues(cell_y_position_var_name.str());
            std::vector<double> computed_cell_y_positions = computed_node_results.GetValues(cell_y_position_var_name.str());
            
            //Comparing expected and computed vector length
            TS_ASSERT_EQUALS(expected_cell_types.size(), computed_cell_types.size());
            TS_ASSERT_EQUALS(expected_cell_x_positions.size(), computed_cell_x_positions.size());
            TS_ASSERT_EQUALS(expected_cell_y_positions.size(), computed_cell_y_positions.size());
            
            //Walkthrough of the expected and computed vectors
            for (unsigned time_step = 0; time_step < expected_cell_types.size(); time_step++)
            {
                TS_ASSERT_EQUALS(expected_cell_types[time_step], computed_cell_types[time_step]);
                TS_ASSERT_DELTA(expected_cell_x_positions[time_step], computed_cell_x_positions[time_step],1e-6);
                TS_ASSERT_DELTA(expected_cell_y_positions[time_step], computed_cell_y_positions[time_step],1e-6);
            }
        }
        
        for (unsigned element=0; element<maxElements; element++)
        {
            std::stringstream nodeA_var_name;
            std::stringstream nodeB_var_name;
            std::stringstream nodeC_var_name;
            nodeA_var_name << "nodeA_" << element;
            nodeB_var_name << "nodeB_" << element;
            nodeC_var_name << "nodeC_" << element;
            
            // Vector of Node A indexes
            std::vector<double> expected_NodeA_numbers = expected_element_results.GetValues(nodeA_var_name.str());
            std::vector<double> computed_NodeA_numbers = computed_element_results.GetValues(nodeA_var_name.str());
            
            // Vector of Node B indexes
            std::vector<double> expected_NodeB_numbers = expected_element_results.GetValues(nodeB_var_name.str());
            std::vector<double> computed_NodeB_numbers = computed_element_results.GetValues(nodeB_var_name.str());
            
            // Vector of Node C indexes
            std::vector<double> expected_NodeC_numbers = expected_element_results.GetValues(nodeC_var_name.str());
            std::vector<double> computed_NodeC_numbers = computed_element_results.GetValues(nodeC_var_name.str());
            
            TS_ASSERT_EQUALS(expected_NodeA_numbers.size(), computed_NodeA_numbers.size());
            TS_ASSERT_EQUALS(expected_NodeB_numbers.size(), computed_NodeB_numbers.size());
            TS_ASSERT_EQUALS(expected_NodeC_numbers.size(), computed_NodeC_numbers.size());
            
            for (unsigned time_step = 0; time_step < expected_NodeA_numbers.size(); time_step++)
            {
                TS_ASSERT_EQUALS(expected_NodeA_numbers[time_step], computed_NodeA_numbers[time_step]);
                TS_ASSERT_EQUALS(expected_NodeB_numbers[time_step], computed_NodeB_numbers[time_step]);
                TS_ASSERT_EQUALS(expected_NodeC_numbers[time_step], computed_NodeC_numbers[time_step]);
            }
            
        }
    }

    
public:

///////// NON-PERIODIC TESTS - These test the spring system and cell birth etc. ----------

    // Test the spring system. The cells in this test are given an intial
    // age of 2.0 so that their springs are at their natural length
    // i.e. we set birth time=-2.0. 
    // The mesh is initially a set of 10 by 10 squares, each square made
    // up of two triangles. The horizontal and vertical edges (springs) are at rest length, the
    // diagonals are two long, so this means the mesh skews to a (sloughed) parallelogram, each
    // triangle trying to become equilateral.
    //
    // If you want to view the results visually set the end time
    // to 24.0 and it will look like a parallelogram.
    // However we keep the simulation time at 1.0 to make
    // the test short.
    void Test2DSpringSystem() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        
        double crypt_length = 10;
        double crypt_width = 10;
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        RandomNumberGenerator::Instance();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);

        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, mesh, FIXED, false, 0.0, 3.0, 6.5, 8.0);
        
        Tissue<2> crypt(mesh, cells);
        TissueSimulation<2> simulator(crypt);    
                
        // destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        simulator.SetEndTime(1.0);
        TS_ASSERT_THROWS_ANYTHING(simulator.Solve());// fails because output directory not set
                
        // destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        simulator.SetOutputDirectory("Crypt2DSprings");

        simulator.SetEndTime(1.0);
        TS_ASSERT_THROWS_ANYTHING(simulator.SetMaxCells(90));
        simulator.SetMaxCells(400);
        TS_ASSERT_THROWS_ANYTHING(simulator.SetMaxElements(90));
        simulator.SetMaxElements(400);
        
        simulator.SetReMeshRule(false);
        simulator.SetNoBirth(true);

        // destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        simulator.Solve();
        std::vector<double> node_0_location = simulator.GetNodeLocation(0);
        TS_ASSERT_DELTA(node_0_location[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(node_0_location[1], 0.0, 1e-12);
        
        CheckAgainstPreviousRun("Crypt2DSprings","results_from_time_0", 400u, 400u);
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    void Test2DHoneycombMeshNotPeriodic() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
       
        int num_cells_depth = 11;
        int num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
                
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
                
        simulator.SetOutputDirectory("Crypt2DHoneycombMesh");
        simulator.SetEndTime(12.0);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);
       
        simulator.Solve();
        
        CheckAgainstPreviousRun("Crypt2DHoneycombMesh","results_from_time_0", 500u, 1000u);
       
        delete p_sloughing_cell_killer;       
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    void TestMonolayer() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
       
        int num_cells_depth = 11;
        int num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true,-1.0);
                
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        // set the first cell to be logged
        crypt.Begin()->SetLogged();

        TissueSimulation<2> simulator(crypt);

        simulator.SetOutputDirectory("Monolayer");
        simulator.SetEndTime(1);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);

        simulator.SetWriteVoronoiData(true,true);
        
        simulator.Solve();
        
        // check writing of voronoi data
        OutputFileHandler handler("Monolayer",false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "VoronoiAreaAndPerimeter.dat";
        TS_ASSERT_EQUALS(system(("cmp " + results_file + " cancer/test/data/Monolayer/VoronoiAreaAndPerimeter.dat").c_str()), 0);
     
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    //////////////////////////////////////////////////////////////////
    // starting with a small mesh with 1 stem cell and the rest
    // differentiated, check the number of cells at the end of the
    // simulation is as expected.
    //////////////////////////////////////////////////////////////////
    void Test2DCorrectCellNumbers() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        RandomNumberGenerator::Instance();
        
        // check the stem cell cycle time is still 24 hrs, otherwise
        // this test might not pass
        TS_ASSERT_DELTA(p_params->GetStemCellCycleTime(), 24, 1e-12);
        TS_ASSERT_DELTA(p_params->GetTransitCellCycleTime(), 12, 1e-12);
        
        int num_cells_width = 7;
        int num_cells_depth = 5;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        double crypt_width = num_cells_width - 1.0;
        double crypt_length = num_cells_depth - 1.0;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = p_mesh->GetNumAllNodes();
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CellType cell_type;
            unsigned generation;
            double birth_time;
            
            if (i==27) // middle of bottom row of cells
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -1;
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = -1; //hours
            }
            
            TissueCell cell(cell_type, HEALTHY, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DSpringsCorrectCellNumbers");
        simulator.SetEndTime(40); //hours
        simulator.SetMaxCells(800);
        simulator.SetMaxElements(800);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);
        
        simulator.Solve();
        
        // now count the number of each type of cell
        unsigned num_stem = 0;
        unsigned num_transit = 0;
        unsigned num_differentiated = 0;
        
        for (Tissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            CellType type = cell_iter->GetCellType();
                
            if (type==STEM)
            {
                num_stem++;
            }
            else if (type==TRANSIT)
            {
                num_transit++;
            }
            else if (type==DIFFERENTIATED)
            {
                num_differentiated++;
            }
            else
            {
                // shouldn't get here
                TS_ASSERT(false);
            }
        }
        
        TS_ASSERT_EQUALS(num_stem, 1u);
        TS_ASSERT_EQUALS(num_transit, 2u);
        
        TS_ASSERT_LESS_THAN(num_differentiated, 25u);
        TS_ASSERT_LESS_THAN(15u, num_differentiated);
        
        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
////////////////////////////////////////////////////////////////////////////
// PERIODIC TESTS
// 
// These test the system as a whole
// 
////////////////////////////////////////////////////////////////////////////
    
    void Test2DPeriodicNightly() throw (Exception)
    {        
        CancerParameters::Instance()->Reset();

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
               
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DPeriodicNightly");
        simulator.SetEndTime(12.0);
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);
       
        simulator.Solve();
        
        // test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
        
        TS_ASSERT_EQUALS(number_of_cells, 85u);
        TS_ASSERT_EQUALS(number_of_nodes, 133u);
        
        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    
    void TestCrypt2DPeriodicWntNightly() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, WNT, true);
        
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
        
        WntGradient::Instance()->SetType(LINEAR);
        WntGradient::Instance()->SetTissue(crypt);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DPeriodicWntNightly");
        
        simulator.SetEndTime(24.0);
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);
       
        simulator.Solve();
        
        // test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)

        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 96u);
        TS_ASSERT_EQUALS(number_of_nodes, 144u);
        
        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
    }
    
    /*
     * \todo - How is this supposed to test the different viscosities??!!
     * 
     * At least run it twice with different viscosities and check nodes are in different locations.
     */
    void dontTestWithMutantCellsUsingDifferentViscosities() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, WNT, true);
                
        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            CellMutationState mutation_state;

            double x = p_mesh->GetNode(i)->GetPoint().rGetLocation()[0];
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            double dist_from_3_6 = sqrt((x-3)*(x-3)+(y-6)*(y-6));
            
            if(dist_from_3_6<1.1)
            {
                mutation_state = APC_TWO_HIT;
            }
            else
            {
                mutation_state = HEALTHY;
            }
            
            cells[i].SetMutationState(mutation_state);
        }
        
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        WntGradient::Instance()->SetType(LINEAR);
        WntGradient::Instance()->SetTissue(crypt);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DPeriodicMutant");        
        simulator.SetEndTime(12.0);
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);
        
        simulator.Solve();
        
        // test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<bool> ghost_cells = crypt.rGetGhostNodes();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
        
        TS_ASSERT_EQUALS(number_of_nodes,ghost_cells.size());
        
        unsigned number_of_cells = 0;
        unsigned number_of_mutant_cells = 0;
        for (Tissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            number_of_cells++;
            if (cell_iter->GetMutationState()==APC_TWO_HIT)
            {
                number_of_mutant_cells++;
            }
        }
        
        TS_ASSERT_EQUALS(number_of_cells, 93u);
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), number_of_cells);
        TS_ASSERT_EQUALS(number_of_nodes, 135u);
        TS_ASSERT_EQUALS(number_of_mutant_cells, 6u);
        
        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
    }
    

    
    void TestRandomDeathWithPeriodicMesh() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        unsigned cells_across = 7;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
              
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DRandomDeathPeriodic");
        simulator.SetEndTime(4.6);
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);

        AbstractCellKiller<2>* p_random_cell_killer = new RandomCellKiller<2>(&crypt, 0.01);
        simulator.AddCellKiller(p_random_cell_killer);
    
        simulator.Solve();
        
        // there should be no cells left after this amount of time
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 0u);
    
        delete p_random_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
  
    // Sloughing with a sloughing cell killer and not turning into ghost nodes
    // on a non-periodic mesh
    void TestSloughingCellKillerOnNonPeriodicCrypt() throw (Exception)
    {
        CancerParameters::Instance()->Reset();

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
              
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DSloughingDeathNonPeriodic");
        simulator.SetEndTime(4.0);
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);

        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);
        
        simulator.Solve();
    
        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }


    void TestSloughingDeathWithPeriodicMesh() throw (Exception)
    {
        CancerParameters::Instance()->Reset();

        unsigned cells_across = 7;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer,true,crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
              
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DSloughingDeathPeriodic");
        simulator.SetEndTime(4.0);
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);

        AbstractCellKiller<2>* p_cell_killer = new SloughingCellKiller(&crypt);
        simulator.AddCellKiller(p_cell_killer);
               
        simulator.Solve();
        
        std::vector<bool> ghost_node_indices_after = crypt.rGetGhostNodes();
        unsigned num_ghosts=0;
        for (unsigned i=0; i < ghost_node_indices_after.size() ; i++)
        {
            if (ghost_node_indices_after[i])
            {
                num_ghosts++;
            }
        }
        
        // Check no new ghost nodes have been created.
        TS_ASSERT_EQUALS(num_ghosts, ghost_node_indices.size());
        // there should be this number of cells left after this amount of time
        // (we have lost two rows of 7 but had a bit of birth too)
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 85u);
    
        delete p_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }


    void TestWithMultipleCellKillers() throw (Exception)
    {
        CancerParameters::Instance()->Reset();

        unsigned cells_across = 7;
        unsigned cells_up = 11;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
              
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("CryptWithMultipleCellKillers");

        // these killers are defined in this test. They kill the first and second
        // available cell, respectively.
        AbstractCellKiller<2>* p_cell_killer1 = new SingleCellCellKiller(&crypt,0);
        AbstractCellKiller<2>* p_cell_killer2 = new SingleCellCellKiller(&crypt,1);

        simulator.AddCellKiller(p_cell_killer1);
        simulator.AddCellKiller(p_cell_killer2);

        // just enough time to kill off all the cells, as two are killed per timestep
        double dt = 0.01;
        unsigned num_cells = crypt.GetNumRealCells();
        simulator.SetDt(dt);
        simulator.SetEndTime(0.5*dt*num_cells);

        simulator.Solve();
        
        std::vector<bool> ghost_node_indices_after = crypt.rGetGhostNodes();
        unsigned num_ghosts=0;        
        for (unsigned i=0; i < ghost_node_indices_after.size() ; i++)
        {
            if (ghost_node_indices_after[i])
            {
                num_ghosts++;
            }
        }
        
        // Check no new ghost nodes have been created.
        TS_ASSERT_EQUALS(num_ghosts, ghost_node_indices.size());

        // all cells should have been removed in this time
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 0u);
    
        delete p_cell_killer1;
        delete p_cell_killer2;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    
    
    void TestMonolayerWithCutoffPointAndNoGhosts() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
       
        int num_cells_depth = 11;
        int num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true,-1.0);
                
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);

        simulator.SetOutputDirectory("MonolayerCutoffPointNoGhosts");
        simulator.SetEndTime(12.0);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);

        simulator.UseCutoffPoint( sqrt(2) ); // root2 is a sensible choice
        
        simulator.Solve();
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }  
};


#endif /*TESTCRYPTSIMULATION2DNIGHTLY_HPP_*/
