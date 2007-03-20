#ifndef TESTCRYPTSIMULATION2DPERIODICNIGHTLY_HPP_
#define TESTCRYPTSIMULATION2DPERIODICNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>
#include "OutputFileHandler.hpp"
#include "CryptSimulation2DPeriodic.hpp"
#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "WntGradient.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"
#include "CryptHoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"

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
        std::cout << "Got to here \n" << std::flush;
        
        
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
        std::cout << "Got to here 2\n" << std::flush;
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
        std::cout << "Got to here 3\n" << std::flush;
    }
    
    
    
public:

    // Test the spring system. There are no cells in this test, therefore no birth, although
    // nodes are sloughed. The mesh is initially a set of 10 by 10 squares, each square made
    // up of two triangles. The horizontal and vertical edges (springs) are at rest length, the
    // diagonals are two long, so this means the mesh skews to a (sloughed) parallelogram, each
    // triangle trying to become equilateral.
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
        
        // throws because start time not set on simulation time
        TS_ASSERT_THROWS_ANYTHING(CryptSimulation2DPeriodic simulator(mesh, std::vector<MeinekeCryptCell>() /*empty*/));
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        CryptSimulation2DPeriodic simulator(mesh, std::vector<MeinekeCryptCell>() /*empty*/);
        
        TS_ASSERT_THROWS_ANYTHING(simulator.Solve());// fails because output directory not set
        
        // destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        simulator.SetOutputDirectory("Crypt2DSprings");
        
        //simulator.SetEndTime(24.0);
        // We need faster tests
        simulator.SetEndTime(1.0);
        TS_ASSERT_THROWS_ANYTHING(simulator.SetMaxCells(90));
        simulator.SetMaxCells(400);
        TS_ASSERT_THROWS_ANYTHING(simulator.SetMaxElements(90));
        simulator.SetMaxElements(400);
        
        simulator.SetReMeshRule(false);
        
        // check an exception is thrown if periodic sim is asked for
        // on a non-periodic mesh
        simulator.SetPeriodicSides(true);
        TS_ASSERT_THROWS_ANYTHING(simulator.Solve());
        
        // destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        simulator.SetPeriodicSides(false);
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
        RandomNumberGenerator *p_random_num_gen=RandomNumberGenerator::Instance();
        
        double crypt_length = 10;
        double crypt_width = 10;
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            
            double y = mesh.GetNode(i)->GetPoint().rGetLocation()[1];
            if (y == 0.0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -p_random_num_gen->ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -p_random_num_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 6.5)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -p_random_num_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 8)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -p_random_num_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
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
        
        CryptSimulation2DPeriodic simulator(mesh,cells);
        simulator.SetOutputDirectory("Crypt2DSpringsFixedBoundaries");
        simulator.SetEndTime(0.2); //hours
        simulator.SetMaxCells(800);
        simulator.SetMaxElements(800);
        simulator.SetFixedBoundaries();
        simulator.SetPeriodicSides(false);
        
        simulator.Solve();
        CheckAgainstPreviousRun("Crypt2DSpringsFixedBoundaries","results_from_time_0", 400u, 800u);
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    void TestWithFixedBirthOnPeriodicMesh() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator *p_random_num_gen=RandomNumberGenerator::Instance();
        
        unsigned cells_across = 7;
        unsigned cells_up = 5;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 3;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
//		double crypt_length = (double)cells_up*(sqrt(3)/2)*crypt_width/(double)cells_across;
        double crypt_length = 4.0;
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = p_mesh->GetNumAllNodes();
        std::cout << "Num Cells = " << num_cells << "\n";
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            
            if (y <= 0.3)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = 0;//-p_random_num_gen->ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -p_random_num_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -p_random_num_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 4)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -p_random_num_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
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
        
        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DPeriodic");
        //simulator.SetEndTime(24.0);
        simulator.SetEndTime(0.2);
        simulator.SetMaxCells(200);
        simulator.SetMaxElements(500);
        simulator.SetNoBirth(false);
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.SetReMeshRule(true);
        
        simulator.Solve();
        CheckAgainstPreviousRun("Crypt2DPeriodic","results_from_time_0", 200u, 500u);
        
        //delete p_mesh;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    void Test2DHoneycombMeshNotPeriodic() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator *p_random_num_gen=RandomNumberGenerator::Instance();
        
        int num_cells_depth = 11;
        int num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;
        
        CryptHoneycombMeshGenerator generator(num_cells_width, num_cells_depth);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = p_mesh->GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            
            if (y == 0.0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -p_random_num_gen->ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -p_random_num_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -p_random_num_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 4)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -p_random_num_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
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
        
        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        //CryptSimulation2D simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DHoneycombMesh");
        simulator.SetEndTime(12.0);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        simulator.SetPeriodicSides(false);
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.Solve();
        
        
        CheckAgainstPreviousRun("Crypt2DHoneycombMesh","results_from_time_0", 500u, 1000u);
        std::cout << "Got to here 4\n" << std::flush;
        
        //delete p_mesh;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        
        std::cout << "Got to here 5\n" << std::flush;
    }
    
    
    //////////////////////////////////////////////////////////////////
    // starting with a small mesh with 1 stem cell and the rest
    // differentiated, check the number of cells at the end of the
    // simulation is as expected.
    //////////////////////////////////////////////////////////////////
    void Test2DCorrectCellNumbers() throw (Exception)
    {
        std::cout << "Got to here 6\n" << std::flush;
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator::Instance();
        
        // check the stem cell cycle time is still 24 hrs, otherwise
        // this test might not pass
        TS_ASSERT_DELTA(p_params->GetStemCellCycleTime(), 24, 1e-12);
        TS_ASSERT_DELTA(p_params->GetTransitCellCycleTime(), 12, 1e-12);
        
        int num_cells_width = 6;
        int num_cells_depth = 5;
        
        CryptHoneycombMeshGenerator generator(num_cells_width, num_cells_depth);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
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
            
            if (i==14) // middle of bottom row of cells
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
        
        CryptSimulation2DPeriodic simulator(*p_mesh,cells);
        simulator.SetOutputDirectory("Crypt2DSpringsCorrectCellNumbers");
        simulator.SetEndTime(40); //hours
        simulator.SetMaxCells(800);
        simulator.SetMaxElements(800);
        
        simulator.SetGhostNodes(ghost_node_indices);
        simulator.SetPeriodicSides(false);
        
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
        TS_ASSERT_LESS_THAN(18, num_differentiated);
        //delete p_mesh;
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    void Test2DPeriodic() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator *p_random_num_gen=RandomNumberGenerator::Instance();
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        double crypt_length = (double)cells_up*(sqrt(3)/2)*crypt_width/(double)cells_across;
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = p_mesh->GetNumAllNodes();
        std::cout << "Num Cells = " << num_cells << "\n";
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            
            if (y <= 0.3)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -p_random_num_gen->ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -p_random_num_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -p_random_num_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 4)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -p_random_num_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
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
        
        
        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DPeriodicNightly");
        
        // Set length of simulation here
        simulator.SetEndTime(12.0);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        // Set to re-mesh and birth
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(false);
        
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.Solve();
        
        std::vector<unsigned> leftBoundary = simulator.GetLeftCryptBoundary();
        std::vector<unsigned> rightBoundary = simulator.GetRightCryptBoundary();
        
        //delete p_mesh;
        
//        std::cout << "Periodic Cell indices at the end of the simulation:\n";
//
//        for(unsigned i=0 ; i<leftBoundary.size(); i++)
//        {
//        	std::cout << "Left " << leftBoundary[i] << ", Right " << rightBoundary[i] << "\n" << std::endl;
//        }

        TS_ASSERT_EQUALS(leftBoundary.size(),12u);
        
        TS_ASSERT_EQUALS(leftBoundary[0], 64u);
        TS_ASSERT_EQUALS(rightBoundary[0], 70u);
        TS_ASSERT_EQUALS(leftBoundary[1], 95u);
        TS_ASSERT_EQUALS(rightBoundary[1], 101u);
        TS_ASSERT_EQUALS(leftBoundary[2], 109u);
        TS_ASSERT_EQUALS(rightBoundary[2], 115u);
        TS_ASSERT_EQUALS(leftBoundary[3], 124u);
        TS_ASSERT_EQUALS(rightBoundary[3], 130u);
        TS_ASSERT_EQUALS(leftBoundary[4], 138u);
        TS_ASSERT_EQUALS(rightBoundary[4], 300u);
        TS_ASSERT_EQUALS(leftBoundary[5], 139u);
        TS_ASSERT_EQUALS(rightBoundary[5], 145u);
        TS_ASSERT_EQUALS(leftBoundary[6], 154u);
        TS_ASSERT_EQUALS(rightBoundary[6], 160u);
        TS_ASSERT_EQUALS(leftBoundary[7], 169u);
        TS_ASSERT_EQUALS(rightBoundary[7], 175u);
        TS_ASSERT_EQUALS(leftBoundary[8], 184u);
        TS_ASSERT_EQUALS(rightBoundary[8], 190u);
        TS_ASSERT_EQUALS(leftBoundary[9], 199u);
        TS_ASSERT_EQUALS(rightBoundary[9], 205u);
        TS_ASSERT_EQUALS(leftBoundary[10], 319u);
        TS_ASSERT_EQUALS(rightBoundary[10], 71u);
        TS_ASSERT_EQUALS(leftBoundary[11], 322u);
        TS_ASSERT_EQUALS(rightBoundary[11], 102u);
        
        //CheckAgainstPreviousRun("Crypt2DPeriodicNightly","results_from_time_0", 500u, 1000u);
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    
    void TestWithWntDependentCells() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        RandomNumberGenerator *p_random_num_gen=RandomNumberGenerator::Instance();
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        double crypt_length = (double)cells_up*(sqrt(3)/2)*crypt_width/(double)cells_across;
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = p_mesh->GetNumAllNodes();
        std::cout << "Num Cells = " << num_cells << std::endl;
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            double typical_wnt_cycle = 16.0;
            if (y <= 0.3)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -p_random_num_gen->ranf()*typical_wnt_cycle; //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -p_random_num_gen->ranf()*typical_wnt_cycle; //hours
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -p_random_num_gen->ranf()*typical_wnt_cycle; //hours
            }
            else if (y < 4)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -p_random_num_gen->ranf()*typical_wnt_cycle; //hours
            }
            else
            {	// There are no fully differentiated cells in a Wnt simulation!
                cell_type = TRANSIT;
                generation = 4;
                birth_time = -p_random_num_gen->ranf()*typical_wnt_cycle; //hours
            }
            WntGradient wnt_gradient(LINEAR);
            double wnt = wnt_gradient.GetWntLevel(y);
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new WntCellCycleModel(wnt,0));
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        std::cout << "Cells Ready." << std::endl;
        
        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DPeriodicWntNightly");
        
        // Set length of simulation here
        simulator.SetEndTime(12.0);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        // Set to re-mesh and birth
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(false);
        simulator.SetWntGradient(LINEAR);
        
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.Solve();
        
        std::vector<unsigned> leftBoundary = simulator.GetLeftCryptBoundary();
        std::vector<unsigned> rightBoundary = simulator.GetRightCryptBoundary();
        std::vector<unsigned> cryptBoundary = simulator.GetCryptBoundary();
        
        //delete p_mesh;
//        std::cout << "Periodic Cell indices at the end of the simulation:\n";
//
//        for(unsigned i=0 ; i<leftBoundary.size(); i++)
//        {
//        	std::cout << "Left " << leftBoundary[i] << ", Right " << rightBoundary[i] << "\n" << std::endl;
//        }
        TS_ASSERT_EQUALS(cryptBoundary.size(),39u);
        TS_ASSERT_EQUALS(leftBoundary.size(),13u);
        
        TS_ASSERT_EQUALS(leftBoundary[0], 64u);
        TS_ASSERT_EQUALS(rightBoundary[0], 70u);
        TS_ASSERT_EQUALS(leftBoundary[1], 78u);
        TS_ASSERT_EQUALS(rightBoundary[1], 84u);
        TS_ASSERT_EQUALS(leftBoundary[2], 95u);
        TS_ASSERT_EQUALS(rightBoundary[2], 101u);
        TS_ASSERT_EQUALS(leftBoundary[3], 109u);
        TS_ASSERT_EQUALS(rightBoundary[3], 115u);
        TS_ASSERT_EQUALS(leftBoundary[4], 124u);
        TS_ASSERT_EQUALS(rightBoundary[4], 130u);
        TS_ASSERT_EQUALS(leftBoundary[5], 139u);
        TS_ASSERT_EQUALS(rightBoundary[5], 145u);
        TS_ASSERT_EQUALS(leftBoundary[6], 154u);
        TS_ASSERT_EQUALS(rightBoundary[6], 160u);
        TS_ASSERT_EQUALS(leftBoundary[7], 169u);
        TS_ASSERT_EQUALS(rightBoundary[7], 175u);
        TS_ASSERT_EQUALS(leftBoundary[8], 184u);
        TS_ASSERT_EQUALS(rightBoundary[8], 190u);
        TS_ASSERT_EQUALS(leftBoundary[9], 199u);
        TS_ASSERT_EQUALS(rightBoundary[9], 205u);
        TS_ASSERT_EQUALS(leftBoundary[10], 214u);
        TS_ASSERT_EQUALS(rightBoundary[10], 220u);
        TS_ASSERT_EQUALS(leftBoundary[11], 300u);
        TS_ASSERT_EQUALS(rightBoundary[11], 116u);
        TS_ASSERT_EQUALS(leftBoundary[12], 325u);
        TS_ASSERT_EQUALS(rightBoundary[12], 85u);
        //CheckAgainstPreviousRun("Crypt2DPeriodicWntNightly","results_from_time_0", 500u, 1000u);
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    // This is strange test -- all cells divide within a quick time, it gives
    // good testing of the periodic boundaries though...
    void TestCrypt2DTysonNovakNightly() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        RandomNumberGenerator *p_random_num_gen=RandomNumberGenerator::Instance();
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        double crypt_length = (double)cells_up*(sqrt(3)/2)*crypt_width/(double)cells_across;
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = p_mesh->GetNumAllNodes();
        std::cout << "Num Cells = " << num_cells << std::endl;
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            double typical_Tyson_Novak_cycle = 1.25;
            if (y <= 0.3)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -p_random_num_gen->ranf()*typical_Tyson_Novak_cycle; //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -p_random_num_gen->ranf()*typical_Tyson_Novak_cycle; //hours
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -p_random_num_gen->ranf()*typical_Tyson_Novak_cycle; //hours
            }
            else if (y < 4)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -p_random_num_gen->ranf()*typical_Tyson_Novak_cycle; //hours
            }
            else
            {  // There are no fully differentiated cells in a Wnt simulation!
                cell_type = TRANSIT;
                generation = 4;
                birth_time = -p_random_num_gen->ranf()*typical_Tyson_Novak_cycle; //hours
            }
            //double wnt = 1.0 - y/p_params->GetCryptLength();
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new TysonNovakCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        std::cout << "Cells Ready." << std::endl;
        
        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DTysonNovakNightly");
        
        // Set length of simulation here
        simulator.SetEndTime(0.5);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        // Set to re-mesh and birth
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(false);
        
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.SetDt(0.001);
        
        simulator.Solve();
        
        //CheckAgainstPreviousRun("Crypt2DPeriodicTysonNovak", 500u, 1000u);
        std::vector<unsigned> leftBoundary = simulator.GetLeftCryptBoundary();
        std::vector<unsigned> rightBoundary = simulator.GetRightCryptBoundary();
        
        std::cout << "Periodic Cell indices at the end of the simulation:\n";
        for (unsigned i=0 ; i<leftBoundary.size(); i++)
        {
            std::cout << "Left " << leftBoundary[i] << ", Right " << rightBoundary[i] << "\n" << std::endl;
        }
        
        TS_ASSERT_EQUALS(leftBoundary.size(),18u);
        
        TS_ASSERT_EQUALS(leftBoundary[0], 64u);
        TS_ASSERT_EQUALS(rightBoundary[0], 70u);
        TS_ASSERT_EQUALS(leftBoundary[1], 79u);
        TS_ASSERT_EQUALS(rightBoundary[1], 85u);
        TS_ASSERT_EQUALS(leftBoundary[2], 94u);
        TS_ASSERT_EQUALS(rightBoundary[2], 100u);
        TS_ASSERT_EQUALS(leftBoundary[3], 108u);
        TS_ASSERT_EQUALS(rightBoundary[3], 114u);
        TS_ASSERT_EQUALS(leftBoundary[4], 123u);
        TS_ASSERT_EQUALS(rightBoundary[4], 129u);
        TS_ASSERT_EQUALS(leftBoundary[5], 153u);
        TS_ASSERT_EQUALS(rightBoundary[5], 130u);
        TS_ASSERT_EQUALS(leftBoundary[6], 154u);
        TS_ASSERT_EQUALS(rightBoundary[6], 160u);
        TS_ASSERT_EQUALS(leftBoundary[7], 168u);
        TS_ASSERT_EQUALS(rightBoundary[7], 337u);
        TS_ASSERT_EQUALS(leftBoundary[8], 169u);
        TS_ASSERT_EQUALS(rightBoundary[8], 175u);
        TS_ASSERT_EQUALS(leftBoundary[9], 184u);
        TS_ASSERT_EQUALS(rightBoundary[9], 190u);
        TS_ASSERT_EQUALS(leftBoundary[10], 199u);
        TS_ASSERT_EQUALS(rightBoundary[10], 205u);
        TS_ASSERT_EQUALS(leftBoundary[11], 229u);
        TS_ASSERT_EQUALS(rightBoundary[11], 235u);
        TS_ASSERT_EQUALS(leftBoundary[12], 339u);
        TS_ASSERT_EQUALS(rightBoundary[12], 221u);
        TS_ASSERT_EQUALS(leftBoundary[13], 349u);
        TS_ASSERT_EQUALS(rightBoundary[13], 144u);
        TS_ASSERT_EQUALS(leftBoundary[14], 357u);
        TS_ASSERT_EQUALS(rightBoundary[14], 101u);
        TS_ASSERT_EQUALS(leftBoundary[15], 378u);
        TS_ASSERT_EQUALS(rightBoundary[15], 71u);
        TS_ASSERT_EQUALS(leftBoundary[16], 382u);
        TS_ASSERT_EQUALS(rightBoundary[16], 206u);
        TS_ASSERT_EQUALS(leftBoundary[17], 393u);
        TS_ASSERT_EQUALS(rightBoundary[17], 114u);
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        
    }
};


#endif /*TESTCRYPTSIMULATION2DPERIODICNIGHTLY_HPP_*/
