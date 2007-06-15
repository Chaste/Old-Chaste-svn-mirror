#ifndef TESTCRYPTSIMULATION2DPERIODICNIGHTLY_HPP_
#define TESTCRYPTSIMULATION2DPERIODICNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "OutputFileHandler.hpp"
#include "MeinekeCryptCell.hpp"
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
#include "Crypt.cpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"

// Possible types of Cell Cycle Model (just for CreateVectorOfCells method)
typedef enum CellCycleType_
{
    FIXED,
    STOCHASTIC,
    WNT,
    TYSONNOVAK
} CellCycleType;


class TestCryptSimulation2DPeriodicNightly : public CxxTest::TestSuite
{
    void CheckAgainstPreviousRun(std::string resultDirectory, std::string resultSet, unsigned maxCells, unsigned maxElements)
    {
        std::cout << "Comparing " << resultDirectory << std::endl << std::flush;
        
        ColumnDataReader computed_node_results = ColumnDataReader(resultDirectory+"/"+resultSet+"/tab_results",
                                                                  "tabulated_node_results",
                                                                  true);
                                                                  
        ColumnDataReader expected_node_results = ColumnDataReader("models/test/data/" + resultDirectory+"Results",
                                                                  "tabulated_node_results",
                                                                  false);
        ColumnDataReader computed_element_results = ColumnDataReader(resultDirectory+"/"+resultSet+"/tab_results",
                                                    "tabulated_element_results",
                                                    true);
                                                    
        ColumnDataReader expected_element_results = ColumnDataReader("models/test/data/" + resultDirectory+"Results",
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
    
    void CreateVectorOfCells(std::vector<MeinekeCryptCell>& rCells, 
                             ConformingTetrahedralMesh<2,2>& rMesh, 
                             CellCycleType cycleType, 
                             bool randomBirthTimes,
                             double y0 = 0.3,
                             double y1 = 2.0,
                             double y2 = 3.0,
                             double y3 = 4.0)
    {
        RandomNumberGenerator *p_random_num_gen=RandomNumberGenerator::Instance();
        unsigned num_cells = rMesh.GetNumNodes();

        AbstractCellCycleModel* p_cell_cycle_model = NULL;
        double typical_transit_cycle_time;
        double typical_stem_cycle_time;
        
        CancerParameters* p_params = CancerParameters::Instance();
        
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;

            double y = rMesh.GetNode(i)->GetPoint().rGetLocation()[1];
            
            if (cycleType==FIXED)
            {
                p_cell_cycle_model = new FixedCellCycleModel();
                typical_transit_cycle_time = p_params->GetTransitCellCycleTime();
                typical_stem_cycle_time = p_params->GetStemCellCycleTime();
            }
            else if (cycleType==STOCHASTIC)
            {
                p_cell_cycle_model = new StochasticCellCycleModel();
                typical_transit_cycle_time = p_params->GetTransitCellCycleTime();
                typical_stem_cycle_time = p_params->GetStemCellCycleTime();
            }
            else if (cycleType==WNT)
            {
                WntGradient wnt_gradient(LINEAR);
                double wnt = wnt_gradient.GetWntLevel(y);
                p_cell_cycle_model = new WntCellCycleModel(wnt,0);
                typical_transit_cycle_time = 16.0;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else if (cycleType==TYSONNOVAK)
            {
                p_cell_cycle_model = new TysonNovakCellCycleModel();
                typical_transit_cycle_time = 1.25;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else
            {
                EXCEPTION("Cell Cycle Type is not recognised");   
            }
            
            
            double birth_time = -2.0;
            
            if (y <= y0)
            {
                cell_type = STEM;
                generation = 0;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_stem_cycle_time; // hours
                }
            }
            else if (y < y1)
            {
                cell_type = TRANSIT;
                generation = 1;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else if (y < y2)
            {
                cell_type = TRANSIT;
                generation = 2;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else if (y < y3)
            {
                cell_type = TRANSIT;
                generation = 3;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else
            {
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
                if(cycleType==WNT || cycleType==TYSONNOVAK)
                {
                    // There are no fully differentiated cells!
                    cell_type = TRANSIT;
                    
                }
                else
                {
                    cell_type = DIFFERENTIATED;
                }                
                generation = 4;
            }

             MeinekeCryptCell cell(cell_type, HEALTHY, generation, p_cell_cycle_model);
            
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            rCells.push_back(cell);
        }
    }
    
public:

///////// NON-PERIODIC TESTS - These test the spring system and cell birth etc. ----------

    // Test the spring system. The cells in this test are given an intial
    // age of 2.0 so that their springs are at their natural length
    // i.e. we set birth time=-2.0. (There used to be no cells in this test)
    // The mesh is initially a set of 10 by 10 squares, each square made
    // up of two triangles. The horizontal and vertical edges (springs) are at rest length, the
    // diagonals are two long, so this means the mesh skews to a (sloughed) parallelogram, each
    // triangle trying to become equilateral.
    //
    // If you want to view the results visually set the end time
    // to 24.0 and it will look like a parallelogram.
    // However we keep the simulation time at 1.0 to make
    // the test short.
    void Test2DSpringSystemWithSloughing() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        
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
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, mesh, FIXED, false, 0.0, 3.0, 6.5, 8.0);
        
        Crypt<2> crypt(mesh, cells);
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
        simulator.SetNoBirth(true);
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
    
    
    void Test2DSpringsFixedBoundaries() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        
        double crypt_length = 10;
        double crypt_width = 10;
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, mesh, FIXED, true, 0.0, 3.0, 6.5, 8.0);
        Crypt<2> crypt(mesh, cells);
        TissueSimulation<2> simulator(crypt);
            
        simulator.SetOutputDirectory("Crypt2DSpringsFixedBoundaries");
        simulator.SetEndTime(0.2); //hours
        simulator.SetMaxCells(800);
        simulator.SetMaxElements(800);
        simulator.SetFixedBoundaries();

        simulator.Solve();
        CheckAgainstPreviousRun("Crypt2DSpringsFixedBoundaries","results_from_time_0", 400u, 800u);
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    void Test2DHoneycombMeshNotPeriodic() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
       
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
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, FIXED, true);
                
        Crypt<2> crypt(*p_mesh, cells);
        TissueSimulation<2> simulator(crypt);
                
        simulator.SetOutputDirectory("Crypt2DHoneycombMesh");
        simulator.SetEndTime(12.0);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.Solve();
        
        //CheckAgainstPreviousRun("Crypt2DHoneycombMesh","results_from_time_0", 500u, 1000u);
       
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    void TestMonolayer() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
       
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
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, FIXED, true,-1.0);
                
        Crypt<2> crypt(*p_mesh, cells);
        TissueSimulation<2> simulator(crypt);
        
        simulator.SetNoSloughing();
        simulator.SetOutputDirectory("Monolayer");
        simulator.SetEndTime(12.0);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.Solve();
        
        //CheckAgainstPreviousRun("Crypt2DHoneycombMesh","results_from_time_0", 500u, 1000u);
       
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
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
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
            
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        Crypt<2> crypt(*p_mesh, cells);
        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DSpringsCorrectCellNumbers");
        simulator.SetEndTime(40); //hours
        simulator.SetMaxCells(800);
        simulator.SetMaxElements(800);
        
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.Solve();
        
        // now count the number of each type of cell
        std::vector<MeinekeCryptCell> cells_after_simulation = simulator.GetCells();
        std::vector<bool> is_ghost_node = simulator.GetGhostNodes();
        
        int num_stem = 0;
        int num_transit = 0;
        int num_differentiated = 0;
        
        for (unsigned index = 0; index<p_mesh->GetNumAllNodes(); index++)
        {
            if (!is_ghost_node[index])  //!mesh.GetNode(index)->IsDeleted())
            {
                CryptCellType type = cells_after_simulation[index].GetCellType();
                
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
        }
        
        TS_ASSERT_EQUALS(num_stem, 1);
        TS_ASSERT_EQUALS(num_transit, 2);
        
        TS_ASSERT_LESS_THAN(num_differentiated, 23);
        TS_ASSERT_LESS_THAN(15, num_differentiated);
        
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
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, FIXED, true);
               
        Crypt<2> crypt(*p_mesh, cells);
        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DPeriodicNightly");
        
        // Set length of simulation here
        simulator.SetEndTime(12.0);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.Solve();
        
        // test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<MeinekeCryptCell> result_cells = simulator.GetCells();
        std::vector<bool> ghost_cells = simulator.GetGhostNodes();
        unsigned number_of_cells = 0;
        unsigned number_of_nodes = result_cells.size();
        
        TS_ASSERT_EQUALS(result_cells.size(),ghost_cells.size());
        
        for (unsigned i=0 ; i<number_of_nodes ; i++)
        {
            if (!ghost_cells[i])
            {
                number_of_cells++;
            }
        }
        TS_ASSERT_EQUALS(number_of_cells, 85u);
        TS_ASSERT_EQUALS(number_of_nodes, 145u);
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    
    void TestCrypt2DPeriodicWntNightly() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
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
        std::vector<MeinekeCryptCell> cells;
        
        CreateVectorOfCells(cells, *p_mesh, WNT, true);
        
        Crypt<2> crypt(*p_mesh, cells);
        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DPeriodicWntNightly");
        
        // Set length of simulation here
        simulator.SetEndTime(24.0);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        simulator.SetWntGradient(LINEAR);
        
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.Solve();
        
        // test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<MeinekeCryptCell> result_cells = simulator.GetCells();
        std::vector<bool> ghost_cells = simulator.GetGhostNodes();
        unsigned number_of_cells = 0;
        unsigned number_of_nodes = result_cells.size();
        
        TS_ASSERT_EQUALS(result_cells.size(),ghost_cells.size());
        
        for (unsigned i=0 ; i<number_of_nodes ; i++)
        {
            if (!ghost_cells[i])
            {
                number_of_cells++;
            }
        }
        TS_ASSERT_EQUALS(number_of_cells, 96u);
        TS_ASSERT_EQUALS(number_of_nodes, 184u);
        
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
   
    void TestWithMutantCellsUsingDifferentViscosities() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
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
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, WNT, true);
                
        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            CryptCellMutationState mutation_state;

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
        
        Crypt<2> crypt(*p_mesh, cells);
        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DPeriodicMutant");
        
        // Set length of simulation here
        simulator.SetEndTime(12.0);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        simulator.SetWntGradient(LINEAR);
        
        simulator.SetGhostNodes(ghost_node_indices);
         
        simulator.Solve();
        
        // test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<MeinekeCryptCell> result_cells = simulator.GetCells();
        std::vector<bool> ghost_cells = simulator.GetGhostNodes();
        unsigned number_of_cells = 0;
        unsigned number_of_nodes = result_cells.size();
        unsigned number_of_mutant_cells = 0;
        TS_ASSERT_EQUALS(result_cells.size(),ghost_cells.size());
        
        for (unsigned i=0 ; i<number_of_nodes ; i++)
        {
            if (!ghost_cells[i])
            {
                number_of_cells++;
                if (result_cells[i].GetMutationState()==APC_TWO_HIT)
                {
                    number_of_mutant_cells++;
                }
            }
        }
        TS_ASSERT_EQUALS(number_of_cells, 94u);
        TS_ASSERT_EQUALS(number_of_nodes, 149u);
        TS_ASSERT_EQUALS(number_of_mutant_cells, 8u);
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    

    // Death on a non-periodic mesh
    // Massive amount of random death eventually leading to every cell being killed off..
    // Note that birth does occur too.
    void TestRandomDeathOnNonPeriodicCrypt() throw (Exception)
    {
        unsigned cells_across = 7;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, FIXED, true);
              
        Crypt<2> crypt(*p_mesh, cells);
        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DRandomDeathNonPeriodic");
        
        // Set length of simulation here
        simulator.SetEndTime(4.0);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
                
        simulator.SetGhostNodes(ghost_node_indices);

        AbstractCellKiller<2>* p_random_cell_killer = new RandomCellKiller<2>(&crypt, 0.01);
        simulator.AddCellKiller(p_random_cell_killer);

        simulator.Solve();
        
        // there should be no cells left after this amount of time
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 0u);
    
        delete p_random_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    void TestRandomDeathWithPeriodicMesh() throw (Exception)
    {
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
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, FIXED, true);
              
        Crypt<2> crypt(*p_mesh, cells);
        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DRandomDeathPeriodic");
                
        // Set length of simulation here
        simulator.SetEndTime(4.6);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
                
        simulator.SetGhostNodes(ghost_node_indices);

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
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, FIXED, true);
              
        Crypt<2> crypt(*p_mesh, cells);
        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DSloughingDeathNonPeriodic");
        
        // Set length of simulation here
        simulator.SetEndTime(4.0);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
                
        simulator.SetGhostNodes(ghost_node_indices);

        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);

        // switch off normal sloughing
        simulator.SetNoSloughing();
        
        simulator.Solve();
    
        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    void TestSloughingDeathWithPeriodicMesh() throw (Exception)
    {
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
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, FIXED, true);
              
        Crypt<2> crypt(*p_mesh, cells);
        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DSloughingDeathPeriodic");
                
        // Set length of simulation here
        simulator.SetEndTime(4.0);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
                
        simulator.SetGhostNodes(ghost_node_indices);

        AbstractCellKiller<2>* p_cell_killer = new SloughingCellKiller(&crypt);
        simulator.AddCellKiller(p_cell_killer);
        
        simulator.SetNoSloughing();
        
        simulator.Solve();
        
        std::vector<bool> ghost_node_indices_after = simulator.GetGhostNodes();
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

};


#endif /*TESTCRYPTSIMULATION2DPERIODICNIGHTLY_HPP_*/
