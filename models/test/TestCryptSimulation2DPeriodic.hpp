#ifndef TESTCRYPTSIMULATION2DPERIODIC_HPP_
#define TESTCRYPTSIMULATION2DPERIODIC_HPP_

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

// Possible types of Cell Cycle Model (just for CreateVectorOfCells method)
typedef enum CellCycleType_
{
    FIXED,
    STOCHASTIC,
    WNT,
    TYSONNOVAK
} CellCycleType;


class TestCryptSimulation2DPeriodic : public CxxTest::TestSuite
{
    void CheckAgainstPreviousRun(std::string resultDirectory,std::string resultSet, unsigned maxCells, unsigned maxElements)
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
            
            
            double birth_time = 0.0;
            
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

    void Test2DCylindrical() throw (Exception)
    {        
        std::string output_directory = "Crypt2DCylindrical";
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 2;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer,true);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, FIXED, true);
               
        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory(output_directory);
        
        // Set length of simulation here
        simulator.SetEndTime(0.1);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        // turn off the old periodic handling - the mesh should deal with it now.
        simulator.SetCylindrical(); // this will disappear once all the methods are overwritten in the mesh class
        simulator.SetNoBirth(true);
        
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.Solve();

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    // This is a rubbish test - all cells start at birthTime = 0.
    // So bizarrely the crypt shrinks as the rest lengths are shortened! 
    // But at least it uses Wnt cell cycle and runs reasonably quickly...
    // For a better test with more randomly distributed cell ages see the Nightly test pack.
    void TestWithWntDependentCells() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells 
        std::vector<MeinekeCryptCell> cells;        
        CreateVectorOfCells(cells, *p_mesh, WNT, false);
                
        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DPeriodicWnt");
        
        // Set length of simulation here
        simulator.SetEndTime(0.3);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        simulator.SetWntGradient(LINEAR);
        
        simulator.SetGhostNodes(ghost_node_indices);
        
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        CheckAgainstPreviousRun("Crypt2DPeriodicWnt","results_from_time_0", 500u, 1000u);
        
        // The following commented out lines provide test results for TestLoad() to check against.
        
//		std::vector<unsigned> left_boundary = simulator.GetLeftCryptBoundary();
//      std::vector<unsigned> right_boundary = simulator.GetRightCryptBoundary();
//
//      std::cout << "Periodic Cell indices at the end of the simulation:\n";
//      for(unsigned i=0 ; i<left_boundary.size(); i++)
//      {
//          std::cout << "Left " << left_boundary[i] << ", Right " << right_boundary[i] << "\n" << std::endl;
//      }
//
//      // A node on the top edge - used for testing load function
//      std::vector<double> node_248_location = simulator.GetNodeLocation(248);
//      std::vector<double> node_219_location = simulator.GetNodeLocation(219);
//
//      std::cout << "Node 248 location: x = " << node_248_location[0] << ",y = "
//                << node_248_location[1] << std::endl;
//      std::cout << "Node 219 location: x = " << node_219_location[0] << ",y = "
//                << node_219_location[1] << std::endl;
        
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    

    // Testing Save (based on previous test)
    void TestSave() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000); 
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, WNT, false);
        
        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DPeriodicWntSaveAndLoad");
        
        // Our full end time is 0.2, here we run for half the time
        simulator.SetEndTime(0.1);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);

        simulator.SetWntGradient(LINEAR);
        
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.Solve();
        
        // save the results..
        simulator.Save();
        
        CheckAgainstPreviousRun("Crypt2DPeriodicWnt","results_from_time_0", 500u, 1000u);
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    

    // Testing Load (based on previous test)
    void TestLoad() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        RandomNumberGenerator::Instance();
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = p_mesh->GetNumAllNodes();
        std::cout << "Num Cells = " << num_cells << std::endl;
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            unsigned generation = 0;
            double wnt_level = 0;
            unsigned mutation_state = 0;
            MeinekeCryptCell cell(STEM, HEALTHY, generation, new WntCellCycleModel(wnt_level,mutation_state));
            cell.SetNodeIndex(i);
            cells.push_back(cell);
        }
        CryptSimulation2DPeriodic simulator(*p_mesh, cells);

        // Load the simulation from the TestSave method above and
        // run it from 0.1 to 0.2
        simulator.Load("Crypt2DPeriodicWntSaveAndLoad",0.1);
        
        simulator.SetEndTime(0.2);
        
        simulator.Solve();
        
        // save that then reload
        // and run from 0.2 to 0.3.
        
        simulator.Save();
        
        simulator.Load("Crypt2DPeriodicWntSaveAndLoad",0.2);
        
        simulator.SetEndTime(0.3);
        
        simulator.Solve();
        
        // compare the results with what they should be after complete of 0.3 hours
        // still worth comparing visually with the output from Crypt2DPeriodicWnt
        std::vector<unsigned> left_boundary = simulator.GetLeftCryptBoundary();
        std::vector<unsigned> right_boundary = simulator.GetRightCryptBoundary();
        
        TS_ASSERT_EQUALS(left_boundary.size(),12u);
        TS_ASSERT_EQUALS(left_boundary[10], 229u);
        TS_ASSERT_EQUALS(right_boundary[10], 221u);
        
        std::vector<double> node_248_location = simulator.GetNodeLocation(248);
        std::vector<double> node_219_location = simulator.GetNodeLocation(219);
        
        TS_ASSERT_DELTA(node_248_location[0], 4.00000 , 1e-5);
        TS_ASSERT_DELTA(node_248_location[1], 8.09225 , 1e-5);
        TS_ASSERT_DELTA(node_219_location[0], 5.00000 , 1e-5);
        TS_ASSERT_DELTA(node_219_location[1], 7.69802 , 1e-5);
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        
        // When the mesh is archived we need a good test here
        // to ensure these results are the same as the ones
        // from TestWithWntDependentCells().
    }
    
    
    
    // This is a rubbish test - all cells start at birthTime = 0.
    // So bizarrely the crypt shrinks as the rest lengths are shortened! But at least it uses Wnt
    // cell cycle and runs reasonably quickly...
    // For a better test with more randomly distributed cell ages see the Nightly test pack.
    void TestWithWntDependentCellsAndAMutation() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, WNT, false);

        // Set a stem cell to be an evil cancer cell and see what happens
        cells[67].SetMutationState(APC_TWO_HIT);
        
        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DMutation");
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        simulator.SetWntGradient(LINEAR);
        
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.SetEndTime(0.05);
        
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    
    // This is strange test -- all cells divide within a quick time, it gives
    // good testing of the periodic boundaries though...
    void TestWithTysonNovakCells() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in T&N
        p_params->SetMaxTransitGenerations(1000);
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, TYSONNOVAK, true);
        
        
        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DPeriodicTysonNovak");
        
        // Set length of simulation here
        simulator.SetEndTime(0.1);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.SetDt(0.001);
        
        simulator.Solve();
        
        std::vector<unsigned> left_boundary = simulator.GetLeftCryptBoundary();
        std::vector<unsigned> right_boundary = simulator.GetRightCryptBoundary();
        
        std::cout << "Periodic Cell indices at the end of the simulation:\n";
        for (unsigned i=0 ; i<left_boundary.size(); i++)
        {
            std::cout << "Left " << left_boundary[i] << ", Right " << right_boundary[i] << "\n" << std::endl;
        }
        
        TS_ASSERT_EQUALS(left_boundary.size(),14u);
        
        TS_ASSERT_EQUALS(left_boundary[0], 64u);
        TS_ASSERT_EQUALS(right_boundary[0], 70u);
        TS_ASSERT_EQUALS(left_boundary[1], 79u);
        TS_ASSERT_EQUALS(right_boundary[1], 85u);
        TS_ASSERT_EQUALS(left_boundary[2], 94u);
        TS_ASSERT_EQUALS(right_boundary[2], 100u);
        TS_ASSERT_EQUALS(left_boundary[3], 108u);
        TS_ASSERT_EQUALS(right_boundary[3], 114u);
        TS_ASSERT_EQUALS(left_boundary[4], 124u);
        TS_ASSERT_EQUALS(right_boundary[4], 341u);
        TS_ASSERT_EQUALS(left_boundary[5], 137u);
        TS_ASSERT_EQUALS(right_boundary[5], 329u);
        TS_ASSERT_EQUALS(left_boundary[6], 138u);
        TS_ASSERT_EQUALS(right_boundary[6], 130u);
        TS_ASSERT_EQUALS(left_boundary[7], 153u);
        TS_ASSERT_EQUALS(right_boundary[7], 144u);
        TS_ASSERT_EQUALS(left_boundary[8], 168u);
        TS_ASSERT_EQUALS(right_boundary[8], 337u);
        TS_ASSERT_EQUALS(left_boundary[9], 169u);
        TS_ASSERT_EQUALS(right_boundary[9], 175u);
        TS_ASSERT_EQUALS(left_boundary[10], 184u);
        TS_ASSERT_EQUALS(right_boundary[10], 190u);
        TS_ASSERT_EQUALS(left_boundary[11], 199u);
        TS_ASSERT_EQUALS(right_boundary[11], 205u);
        TS_ASSERT_EQUALS(left_boundary[12], 229u);
        TS_ASSERT_EQUALS(right_boundary[12], 235u);
        TS_ASSERT_EQUALS(left_boundary[13], 339u);
        TS_ASSERT_EQUALS(right_boundary[13], 221u);
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    

    void TestCalculateCryptBoundaries()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(0.0,-2.0,0.0) ;
        //Create Vector of ghost nodes
        std::vector<unsigned> ghost_node_indices;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint().rGetLocation()[0];
            double y = mesh.GetNode(i)->GetPoint().rGetLocation()[1];
            if ((x<2.0)||(x>8.0)||(y>6.0)||(y<0.0))
            {
                ghost_node_indices.push_back(i);
            }
        }
        
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->SetCryptLength(6.0);
        p_params->SetCryptWidth(6.0);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        CryptSimulation2DPeriodic simulator(mesh);
        simulator.SetGhostNodes(ghost_node_indices);
        
        simulator.CalculateCryptBoundary();
        
        //simulator.DetectNaughtyCellsJoiningPeriodicEdges();
        
        std::vector<unsigned> calculated_boundary_nodes  = simulator.GetCryptBoundary();
        std::vector<unsigned> actual_boundary_nodes(24);
        
        actual_boundary_nodes[0] = 24;
        actual_boundary_nodes[1] = 25;
        actual_boundary_nodes[2] = 26;
        actual_boundary_nodes[3] = 27;
        actual_boundary_nodes[4] = 28;
        actual_boundary_nodes[5] = 29;
        actual_boundary_nodes[6] = 30;
        actual_boundary_nodes[7] = 35;
        actual_boundary_nodes[8] = 41;
        actual_boundary_nodes[9] = 46;
        actual_boundary_nodes[10] = 52;
        actual_boundary_nodes[11] = 57;
        actual_boundary_nodes[12] = 63;
        actual_boundary_nodes[13] = 68;
        actual_boundary_nodes[14] = 74;
        actual_boundary_nodes[15] = 79;
        actual_boundary_nodes[16] = 85;
        actual_boundary_nodes[17] = 90;
        actual_boundary_nodes[18] = 91;
        actual_boundary_nodes[19] = 92;
        actual_boundary_nodes[20] = 93;
        actual_boundary_nodes[21] = 94;
        actual_boundary_nodes[22] = 95;
        actual_boundary_nodes[23] = 96;
        
        TS_ASSERT_EQUALS(actual_boundary_nodes.size(),calculated_boundary_nodes.size());
        
        
        for (unsigned i=0; i<calculated_boundary_nodes.size(); i++)
        {
            TS_ASSERT_EQUALS(actual_boundary_nodes[i],calculated_boundary_nodes[i]);
        }
        
        std::vector<unsigned> calculated_left_boundary_nodes = simulator.GetLeftCryptBoundary();
        std::vector<unsigned> calculated_right_boundary_nodes = simulator.GetRightCryptBoundary();
        
        std::vector<unsigned> actual_left_boundary_nodes(7);
        
        actual_left_boundary_nodes[0] = 24;
        actual_left_boundary_nodes[1] = 35;
        actual_left_boundary_nodes[2] = 46;
        actual_left_boundary_nodes[3] = 57;
        actual_left_boundary_nodes[4] = 68;
        actual_left_boundary_nodes[5] = 79;
        actual_left_boundary_nodes[6] = 90;
        
        std::vector<unsigned> actual_right_boundary_nodes(7);
        
        actual_right_boundary_nodes[0] = 24+6;
        actual_right_boundary_nodes[1] = 35+6;
        actual_right_boundary_nodes[2] = 46+6;
        actual_right_boundary_nodes[3] = 57+6;
        actual_right_boundary_nodes[4] = 68+6;
        actual_right_boundary_nodes[5] = 79+6;
        actual_right_boundary_nodes[6] = 90+6;
        
        
        TS_ASSERT_EQUALS(actual_left_boundary_nodes.size(),calculated_left_boundary_nodes.size());
        TS_ASSERT_EQUALS(actual_right_boundary_nodes.size(),calculated_right_boundary_nodes.size());
        
        
        for (unsigned i=0; i<calculated_left_boundary_nodes.size(); i++)
        {
            TS_ASSERT_EQUALS(actual_left_boundary_nodes[i],calculated_left_boundary_nodes[i]);
            TS_ASSERT_EQUALS(actual_right_boundary_nodes[i],calculated_right_boundary_nodes[i]);
        }
        SimulationTime::Destroy();
    }
    
    
    void TestPrivateFunctionsOf2DCryptSimulation() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator::Instance();
        
        double crypt_length = 9.3;
        double crypt_width = 10.0;
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);        
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, mesh, FIXED, false, 0.0, 3.0, 6.5, 8.0);
        
        cells[60].SetBirthTime(-50.0);
        
        CryptSimulation2DPeriodic simulator(mesh,cells);
        
        simulator.SetFixedBoundaries();
        simulator.SetPeriodicSides(false);
        
        unsigned num_deaths = simulator.DoCellRemoval();
        unsigned num_births = simulator.DoCellBirth();
        Node<2> *p_node = mesh.GetNode(60);
        Element<2,2>* p_element = simulator.FindElementForBirth(p_node, 60u);
                                                                
        TS_ASSERT_EQUALS(num_births, 1u);
        TS_ASSERT_EQUALS(num_deaths,11u);
        TS_ASSERT_EQUALS(p_element->GetIndex(),110u);
        
        p_params->SetCryptLength(10.1);
        CryptSimulation2DPeriodic simulator2(mesh,cells);
        
        simulator2.SetFixedBoundaries();
        simulator2.SetPeriodicSides(false);
        
        num_deaths = simulator2.DoCellRemoval();
        TS_ASSERT_EQUALS(num_deaths,0u);
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    
    void TestPrivateFunctionsOf2DCryptSimulationOnHoneycombMesh() throw (Exception)
    {
        /*
         ************************************************************************
         ************************************************************************ 
         *     Set up a simulation class to run the individual tests on.
         ************************************************************************
         ************************************************************************ 
         */
        unsigned cells_across2 = 6;
        unsigned cells_up2 = 5;
        double crypt_width2 = 6.0;
        unsigned thickness_of_ghost_layer2 = 3;
        
        CryptHoneycombMeshGenerator generator(cells_across2, cells_up2, crypt_width2,thickness_of_ghost_layer2);
        ConformingTetrahedralMesh<2,2>* p_mesh2=generator.GetMesh();
        std::vector<unsigned> ghost_node_indices2 = generator.GetGhostNodeIndices();
        unsigned num_cells2 = p_mesh2->GetNumAllNodes();
        
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator::Instance();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells by iterating through the mesh nodes
        std::vector<MeinekeCryptCell> cells2;
        for (unsigned i=0; i<num_cells2; i++)
        {
            double birth_time;
            CryptCellType cell_type;
            unsigned generation;
            double y = p_mesh2->GetNode(i)->GetPoint().rGetLocation()[1];
            if (y == 0.0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -2.0; //hours - doesn't matter for stem cell;
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -2.0; //hours
            }
            else if (y < 6.5)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -2.0;  //hours
            }
            else if (y < 8)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -2.0;  //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = -2.0;  //hours
            }
            
            CryptCellMutationState mutation_state;
            if(i!=60)
            {
                mutation_state = HEALTHY;
            }
            else
            {
                mutation_state = APC_TWO_HIT;
            }
            
            WntGradient wnt_gradient(LINEAR);
            double wnt = wnt_gradient.GetWntLevel(y);
            
            
            
            MeinekeCryptCell cell(cell_type, mutation_state, generation, new WntCellCycleModel(wnt,0));
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells2.push_back(cell);
        }
        
        CryptSimulation2DPeriodic simulator3(*p_mesh2,cells2);
        simulator3.SetGhostNodes(ghost_node_indices2);
        simulator3.SetPeriodicSides(false);
        
        simulator3.SetMaxCells(400);
        simulator3.SetMaxElements(400);
        simulator3.SetOutputDirectory("TestPrivateMemberDirectory");
        std::string output_directory = "TestPrivateMemberDirectory";
        ColumnDataWriter tabulated_node_writer(output_directory+"/tab_results", "tabulated_node_results");
        ColumnDataWriter tabulated_element_writer(output_directory+"/tab_results", "tabulated_element_results");
        
        /*
         ************************************************************************
         ************************************************************************ 
         *  Test results file writers
         ************************************************************************
         ************************************************************************ 
         */
        node_writer_ids_t node_writer_ids;
        TS_ASSERT_THROWS_NOTHING(simulator3.SetupNodeWriter(tabulated_node_writer, node_writer_ids));
        
        element_writer_ids_t element_writer_ids;
        TS_ASSERT_THROWS_NOTHING(simulator3.SetupElementWriter(tabulated_element_writer, element_writer_ids));
        
        OutputFileHandler output_file_handler(output_directory);
        out_stream p_node_file = output_file_handler.OpenOutputFile("results.viznodes");
        out_stream p_element_file = output_file_handler.OpenOutputFile("results.vizelements");
        unsigned tabulated_output_counter = 0;
        
        simulator3.WriteResultsToFiles(tabulated_node_writer, node_writer_ids,
                                       tabulated_element_writer, element_writer_ids,
                                       *p_node_file, *p_element_file,
                                       tabulated_output_counter==0,
                                       true);
        /*
         ************************************************************************
         ************************************************************************ 
         *  Test Calculate Velocities on each node
         ************************************************************************
         ************************************************************************ 
         */
                
        std::vector<c_vector<double, 2> > velocities_on_each_node(p_mesh2->GetNumAllNodes());
        
        velocities_on_each_node = simulator3.CalculateVelocitiesOfEachNode();
        bool is_a_ghost_node;

        
        for (unsigned i=0; i<p_mesh2->GetNumAllNodes(); i++)
        {
            is_a_ghost_node = false;
            for (unsigned j=0; j<ghost_node_indices2.size(); j++)
            {
                if (ghost_node_indices2[j]==i)
                {
                    is_a_ghost_node = true;
                }
            }
            if (!is_a_ghost_node)
            {
                TS_ASSERT_DELTA(velocities_on_each_node[i][0], 0.0, 1e-4);
                TS_ASSERT_DELTA(velocities_on_each_node[i][1], 0.0, 1e-4);
            }
        }
        
        // Move a node along the x-axis and calculate the force exerted on a neighbour
        c_vector<double,2> old_point = p_mesh2->GetNode(59)->rGetLocation();
        Point<2> new_point;
        new_point.rGetLocation()[0] = old_point[0]+0.5;
        new_point.rGetLocation()[1] = old_point[1];
        p_mesh2->SetNode(59, new_point, false);
        velocities_on_each_node = simulator3.CalculateVelocitiesOfEachNode();
        TS_ASSERT_DELTA(velocities_on_each_node[60][0], 0.5*p_params->GetSpringStiffness()/p_params->GetDampingConstantMutant(), 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[60][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(velocities_on_each_node[59][0], (-3+4.0/sqrt(7))*p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal(), 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[59][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(velocities_on_each_node[58][0], 0.5*p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal(), 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[58][1], 0.0, 1e-4);
        /*
         ************************************************************************
         ************************************************************************ 
         *  Test Calculate force on a spring
         ************************************************************************
         ************************************************************************ 
         */
        
        c_vector<double,2> force_on_spring ; // between nodes 59 and 60
        
        // Find one of the elements that nodes 59 and 60 live on
        Point<2> new_point2;
        new_point2.rGetLocation()[0] = new_point[0] + 0.01;
        new_point2.rGetLocation()[1] = new_point[1] + 0.01 ;
        
        unsigned elem_index = p_mesh2->GetContainingElementIndex(new_point2,false);
        Element<2,2>* p_element = p_mesh2->GetElement(elem_index);
        
        force_on_spring = simulator3.CalculateForceInThisSpring(p_element,1,0);
        
        TS_ASSERT_DELTA(force_on_spring[0], 0.5*p_params->GetSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);
        
        /*
         ************************************************************************
         ************************************************************************ 
         *  Test Calculate force on a boundary spring
         ************************************************************************
         ************************************************************************ 
         */
        
        c_vector<double,2> force_on_boundary_spring ; 
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator it;
        
        it = p_mesh2->GetBoundaryElementIteratorBegin();
               
        // find a boundary element which is not deleted
        BoundaryElement<1,2>* p_edge;
        bool found=false;
        while (it != p_mesh2->GetBoundaryElementIteratorEnd() && !found )
        {   
            p_edge=(*it);
            found = !p_edge->IsDeleted();    
            it++;
        }
        
        TS_ASSERT(it != p_mesh2->GetBoundaryElementIteratorEnd());
        
        force_on_boundary_spring = simulator3.CalculateForceInThisBoundarySpring(p_edge);
        TS_ASSERT_DELTA(force_on_boundary_spring[0], 0.0, 1e-10);
        TS_ASSERT_DELTA(force_on_boundary_spring[1], 0.0, 1e-10);
        
        // move one of the nodes of the boundary element
        unsigned node_index = p_edge->GetNode(0)->GetIndex();
        c_vector<double,2> node_location = p_edge->GetNode(0)->rGetLocation();
        node_location[0] += 0.5;
        p_mesh2->SetNode(node_index, Point<2>(node_location), false);
        
        // check force
        force_on_boundary_spring = simulator3.CalculateForceInThisBoundarySpring(p_edge);
        TS_ASSERT_DELTA(force_on_boundary_spring[0], -0.5*p_params->GetSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(force_on_boundary_spring[1], 0.0, 1e-4);

        
        /*
         ************************************************************************
         ************************************************************************ 
         *  Test UpdateNodePositions
         ************************************************************************
         ************************************************************************ 
         */
        
        Point<2> point_of_node60 = p_mesh2->GetNode(60)->rGetLocation();
        
        simulator3.SetDt(0.01);
        simulator3.UpdateNodePositions(velocities_on_each_node);
        
        TS_ASSERT_DELTA(p_mesh2->GetNode(60)->rGetLocation()[0],point_of_node60.rGetLocation()[0]+force_on_spring[0]/p_params->GetDampingConstantMutant() *0.01, 1e-4);
        TS_ASSERT_DELTA(p_mesh2->GetNode(60)->rGetLocation()[1],point_of_node60.rGetLocation()[1], 1e-4);
        
        /*
         ************************************************************************
         ************************************************************************ 
         * Test UpdateCellTypes 
         ************************************************************************
         ************************************************************************ 
         */
        
        std::vector<MeinekeCryptCell> cells3;
        simulator3.SetWntGradient(LINEAR);
        simulator3.UpdateCellTypes();
        cells3 = simulator3.GetCells();
        
        std::vector<bool> is_node_a_ghost = simulator3.GetGhostNodes();
        
        for (unsigned i=0; i<num_cells2; i++)
        {
            if (!is_node_a_ghost[i])
            {
                CryptCellType cell_type;
                cell_type = cells3[i].GetCellType();
                if (!cell_type==STEM)
                {
                    //std::cout << "Cell type = " << cell_type << std::endl;
                    WntCellCycleModel *p_this_model = static_cast<WntCellCycleModel*>(cells3[i].GetCellCycleModel());
                    double beta_cat_level = p_this_model->GetProteinConcentrations()[6]+ p_this_model->GetProteinConcentrations()[7];
                    //std::cout << "Cell " << i << ", beta-cat = " << beta_cat_level << std::endl;
                    if (beta_cat_level > 0.4127)
                    {
                        TS_ASSERT_EQUALS(cell_type,TRANSIT);
                    }
                    else
                    {
                        TS_ASSERT_EQUALS(cell_type,DIFFERENTIATED);
                    }
                }
            }
        }
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    

    
};

#endif /*TESTCRYPTSIMULATION2DPERIODIC_HPP_*/
