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

class TestCryptSimulation2DPeriodic : public CxxTest::TestSuite
{
	void CheckAgainstPreviousRun(std::string resultDirectory, unsigned maxCells, unsigned maxElements)
    {
        std::cout << "Comparing " << resultDirectory << std::endl << std::flush;
        
        ColumnDataReader computed_node_results = ColumnDataReader(resultDirectory+"Results",
                                                             "tabulated_node_results",
                                                             true);
                                                             
        ColumnDataReader expected_node_results = ColumnDataReader("models/test/data/" + resultDirectory+"Results",
                                                             "tabulated_node_results",
                                                             false);
        ColumnDataReader computed_element_results = ColumnDataReader(resultDirectory+"Results",
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

        RandomNumberGenerator random_num_gen; // passed into crypt sim for coverage
        
        // throws because start time not set on simulation time
        TS_ASSERT_THROWS_ANYTHING(CryptSimulation2DPeriodic simulator(mesh, std::vector<MeinekeCryptCell>() /*empty*/, &random_num_gen));
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        CryptSimulation2DPeriodic simulator(mesh, std::vector<MeinekeCryptCell>() /*empty*/, &random_num_gen);

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
              
        CheckAgainstPreviousRun("Crypt2DSprings", 400u, 400u);
    }
    
    
    void Test2DSpringsFixedBoundaries() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator random_num_gen;
        
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
                birth_time = -random_num_gen.ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 6.5)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 8)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
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
        CheckAgainstPreviousRun("Crypt2DSpringsFixedBoundaries", 400u, 800u);
    }
	
	void TestWithFixedBirthOnPeriodicMesh() throw (Exception)
    {
    	CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator random_num_gen;
        
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
                birth_time = 0;//-random_num_gen.ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 4)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
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
		
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        CheckAgainstPreviousRun("Crypt2DPeriodic", 200u, 500u);
    }
    
	// This is a rubbish test - all cells start at birthTime = 0.
	// So bizarrely the crypt shrinks as the rest lengths are shortened! But at least it uses Wnt
	// cell cycle and runs reasonably quickly...
	// For a better test with more randomly distributed cell ages see the Nightly test pack.
    void TestWithWntDependentCells() throw (Exception)
    {
    	CancerParameters *p_params = CancerParameters::Instance();
		// There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        RandomNumberGenerator random_num_gen;
        
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
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else if (y < 3)
            {
            	cell_type = TRANSIT;
                generation = 2;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else if (y < 4)
            {
            	cell_type = TRANSIT;
                generation = 3;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else
            {	
                // There are no fully differentiated cells in a Wnt simulation!
                cell_type = TRANSIT;
                generation = 4;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            WntGradient wnt_gradient(LINEAR);
            double wnt = wnt_gradient.GetWntLevel(y);
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new WntCellCycleModel(wnt,0));
            cell.SetNodeIndex(i);
            cell.SetBirthTime(0.0);
            cells.push_back(cell);
        }
        std::cout << "Cells Ready." << std::endl;

        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DPeriodicWnt");

		// Set length of simulation here
        simulator.SetEndTime(0.2);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        // Set to re-mesh and birth
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(false);
        simulator.SetWntGradient(LINEAR);
        
        simulator.SetGhostNodes(ghost_node_indices);
        
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        CheckAgainstPreviousRun("Crypt2DPeriodicWnt", 500u, 1000u);
    }
    
    // Testing Save (based on previous test)
    void TestSave() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        RandomNumberGenerator random_num_gen;
        
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
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else if (y < 4)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else
            {   
                // There are no fully differentiated cells in a Wnt simulation!
                cell_type = TRANSIT;
                generation = 4;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            WntGradient wnt_gradient(LINEAR);
            double wnt = wnt_gradient.GetWntLevel(y);
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new WntCellCycleModel(wnt,0));
            cell.SetNodeIndex(i);
            cell.SetBirthTime(0.0);
            cells.push_back(cell);
        }
        std::cout << "Cells Ready." << std::endl;

        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DPeriodicWntSaveAndLoad");

        // Our full end time is 0.2, here we run for half the time
        simulator.SetEndTime(0.1);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        // Set to re-mesh and birth
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(false);
        simulator.SetWntGradient(LINEAR);
        
        simulator.SetGhostNodes(ghost_node_indices);
        

        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // save the results..        
        simulator.Save();

        CheckAgainstPreviousRun("Crypt2DPeriodicWnt", 500u, 1000u);
    }
    
    // Testing Load (based on previous test)
    void xTestLoad() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        RandomNumberGenerator random_num_gen;
        
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
        
        std::cout << "crypt width = " << p_params->GetCryptWidth() << "\n" << std::flush;
        std::cout << "crypt length = " << p_params->GetCryptLength() << "\n" << std::flush;
        
        
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
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else if (y < 4)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else
            {   
                // There are no fully differentiated cells in a Wnt simulation!
                cell_type = TRANSIT;
                generation = 4;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            WntGradient wnt_gradient(LINEAR);
            double wnt = wnt_gradient.GetWntLevel(y);
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new WntCellCycleModel(wnt,0));
            cell.SetNodeIndex(i);
            cell.SetBirthTime(0.0);
            cells.push_back(cell);
        }
        std::cout << "Cells Ready." << std::endl;

        CryptSimulation2DPeriodic simulator(*p_mesh, cells);

        simulator.Load();

        simulator.SetEndTime(0.2);


        simulator.Solve();

        CheckAgainstPreviousRun("Crypt2DPeriodicWnt", 500u, 1000u);
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
        RandomNumberGenerator random_num_gen;
        
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
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else if (y < 3)
            {
            	cell_type = TRANSIT;
                generation = 2;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else if (y < 4)
            {
            	cell_type = TRANSIT;
                generation = 3;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else
            {	
                // There are no fully differentiated cells in a Wnt simulation!
                cell_type = TRANSIT;
                generation = 4;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            WntGradientType this_type=LINEAR;
            WntGradient wnt_gradient(this_type);
            double wnt = wnt_gradient.GetWntLevel(y);
            CryptCellMutationState this_state=HEALTHY;
            MeinekeCryptCell cell(cell_type, this_state, generation, new WntCellCycleModel(wnt,0));
            cell.SetNodeIndex(i);
            cell.SetBirthTime(0.0);
            cells.push_back(cell);
        }
        // Set a stem cell to be an evil cancer cell and see what happens
        cells[67].SetMutationState(APC_TWO_HIT);
        
        std::cout << "Cells Ready." << std::endl;

        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DMutation");
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        // Set to re-mesh and birth
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(false);
        simulator.SetWntGradient(LINEAR);
        
        simulator.SetGhostNodes(ghost_node_indices);
                
        simulator.SetEndTime(0.05);

        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }
    
    // This is strange test -- all cells divide within a quick time, it gives
    // good testing of the periodic boundaries though...
    void TestWithTysonNovakCells() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        RandomNumberGenerator random_num_gen;
        
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
                birth_time = -random_num_gen.ranf()*typical_Tyson_Novak_cycle; //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -random_num_gen.ranf()*typical_Tyson_Novak_cycle; //hours
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -random_num_gen.ranf()*typical_Tyson_Novak_cycle; //hours
            }
            else if (y < 4)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -random_num_gen.ranf()*typical_Tyson_Novak_cycle; //hours
            }
            else
             {  // There are no fully differentiated cells in a Wnt simulation!
                cell_type = TRANSIT;
                generation = 4;
                birth_time = -random_num_gen.ranf()*typical_Tyson_Novak_cycle; //hours
            }
            //double wnt = 1.0 - y/p_params->GetCryptLength();
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new TysonNovakCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        std::cout << "Cells Ready." << std::endl;

        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DPeriodicTysonNovak");

        // Set length of simulation here
        simulator.SetEndTime(0.1);
        
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
        for(unsigned i=0 ; i<leftBoundary.size(); i++)
        {
        	std::cout << "Left " << leftBoundary[i] << ", Right " << rightBoundary[i] << "\n" << std::endl;
        }
        
        TS_ASSERT_EQUALS(leftBoundary.size(),14u);
      		
       	TS_ASSERT_EQUALS(leftBoundary[0], 64u);
		TS_ASSERT_EQUALS(rightBoundary[0], 70u);
		TS_ASSERT_EQUALS(leftBoundary[1], 79u);
		TS_ASSERT_EQUALS(rightBoundary[1], 85u);
		TS_ASSERT_EQUALS(leftBoundary[2], 94u);
		TS_ASSERT_EQUALS(rightBoundary[2], 100u);
		TS_ASSERT_EQUALS(leftBoundary[3], 108u);
		TS_ASSERT_EQUALS(rightBoundary[3], 114u);
		TS_ASSERT_EQUALS(leftBoundary[4], 124u);
		TS_ASSERT_EQUALS(rightBoundary[4], 341u);
		TS_ASSERT_EQUALS(leftBoundary[5], 137u);
		TS_ASSERT_EQUALS(rightBoundary[5], 329u);
		TS_ASSERT_EQUALS(leftBoundary[6], 138u);
		TS_ASSERT_EQUALS(rightBoundary[6], 130u);
		TS_ASSERT_EQUALS(leftBoundary[7], 153u);
		TS_ASSERT_EQUALS(rightBoundary[7], 144u);
		TS_ASSERT_EQUALS(leftBoundary[8], 168u);
		TS_ASSERT_EQUALS(rightBoundary[8], 337u);
		TS_ASSERT_EQUALS(leftBoundary[9], 169u);
		TS_ASSERT_EQUALS(rightBoundary[9], 175u);
		TS_ASSERT_EQUALS(leftBoundary[10], 184u);
		TS_ASSERT_EQUALS(rightBoundary[10], 190u);
		TS_ASSERT_EQUALS(leftBoundary[11], 199u);
		TS_ASSERT_EQUALS(rightBoundary[11], 205u);
		TS_ASSERT_EQUALS(leftBoundary[12], 229u);
		TS_ASSERT_EQUALS(rightBoundary[12], 235u);
		TS_ASSERT_EQUALS(leftBoundary[13], 339u);
		TS_ASSERT_EQUALS(rightBoundary[13], 221u);
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
        
        
        for(unsigned i=0; i<calculated_boundary_nodes.size(); i++)
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
        
        
        for(unsigned i=0; i<calculated_left_boundary_nodes.size(); i++)
        {
            TS_ASSERT_EQUALS(actual_left_boundary_nodes[i],calculated_left_boundary_nodes[i]);
            TS_ASSERT_EQUALS(actual_right_boundary_nodes[i],calculated_right_boundary_nodes[i]);
        }
    }
    
    
    
    void TestPrivateFunctionsOf2DCryptSimulation() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator random_num_gen;
        
        double crypt_length = 9.3;
        double crypt_width = 10.0;
        
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
            //double birth_time;
            
            double y = mesh.GetNode(i)->GetPoint().rGetLocation()[1];
            if (y == 0.0)
            {
                cell_type = STEM;
                generation = 0;
                //birth_time = -random_num_gen.ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 1;
                //birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 6.5)
            {
                cell_type = TRANSIT;
                generation = 2;
               // birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 8)
            {
                cell_type = TRANSIT;
                generation = 3;
                //birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                //birth_time = -1; //hours
            }
            
            
            
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            if ( i == 60u)
            {
               cell.SetBirthTime(-50.0 );  
            }
            
            cells.push_back(cell);
        }
        
        
        CryptSimulation2DPeriodic simulator(mesh,cells,&random_num_gen);
        //simulator.SetOutputDirectory("ghg");
        //simulator.SetMaxCells(200);
        //simulator.SetMaxElements(400);
        simulator.SetFixedBoundaries();
        simulator.SetPeriodicSides(false);
        // sim time destroy needs to be below where we move cells around 
        // (because it makes new ones from copies of old ones??)
        
        //std::cout << "d test \n " << std::endl;
//        for (unsigned i=0; i<mesh.GetNumAllNodes(); i++)
//        {
//            std::cout << " i " << i << "forces x " << forces_on_each_node[i][0] << ", y " << forces_on_each_node[i][1] << "\n" << std::endl;
//            //TS_ASSERT_DELTA(forces_on_each_node[i][0], 0.0, 1e-8);
//            //TS_ASSERT_DELTA(forces_on_each_node[i][1], 0.0, 1e-8);
//        }
        
        unsigned num_deaths = simulator.DoCellRemoval();
        unsigned num_births = simulator.DoCellBirth();
        Node<2> *p_node = mesh.GetNode(60);
        Element<2,2>* p_element = simulator.FindElementForBirth(p_node, 60u,
                                      false, 0u);
        
        TS_ASSERT_EQUALS(num_births, 1u);
        TS_ASSERT_EQUALS(num_deaths,11u);
        TS_ASSERT_EQUALS(p_element->GetIndex(),110u);
        
        p_params->SetCryptLength(10.1);
        CryptSimulation2DPeriodic simulator2(mesh,cells,&random_num_gen);
        //simulator.SetOutputDirectory("ghg");
        
        simulator2.SetFixedBoundaries();
        simulator2.SetPeriodicSides(false);
        // sim time destroy needs to be below where we move cells around 
        // (because it makes new ones from copies of old ones??)
        
        
        num_deaths = simulator2.DoCellRemoval();
        TS_ASSERT_EQUALS(num_deaths,0u);
        
        SimulationTime::Destroy();
    }   
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void TestPrivateFunctionsOf2DCryptSimulationOnHoneycombMesh() throw (Exception)
  {      
        unsigned cells_across2 = 6;
        unsigned cells_up2 = 5;
        double crypt_width2 = 6.0;
        unsigned thickness_of_ghost_layer2 = 3;
        
        CryptHoneycombMeshGenerator generator(cells_across2, cells_up2, crypt_width2,thickness_of_ghost_layer2);  
        ConformingTetrahedralMesh<2,2>* p_mesh2=generator.GetMesh(); 
        std::vector<unsigned> ghost_node_indices2 = generator.GetGhostNodeIndices(); 
        unsigned num_cells2 = p_mesh2->GetNumAllNodes();
        
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator random_num_gen;
        
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
            
            
            
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new WntCellCycleModel(0.0));
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells2.push_back(cell);
        }
        
        CryptSimulation2DPeriodic simulator3(*p_mesh2,cells2,&random_num_gen);
        simulator3.SetGhostNodes(ghost_node_indices2);
        simulator3.SetPeriodicSides(false);
        
        simulator3.SetMaxCells(400);
        simulator3.SetMaxElements(400);
        simulator3.SetOutputDirectory("TestPrivateMemberDirectory");
        std::string output_directory = "TestPrivateMemberDirectory";
        ColumnDataWriter tabulated_node_writer(output_directory+"Results", "tabulated_node_results");
        ColumnDataWriter tabulated_element_writer(output_directory+"Results", "tabulated_element_results");
        
        
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
                                
        
        
        std::vector<std::vector<double> > forces_on_each_node(p_mesh2->GetNumAllNodes());
        
        forces_on_each_node = simulator3.CalculateForcesOnEachNode();
        //std::cout << "d test \n " << std::endl;
        bool is_a_ghost_node;
        
        
        
        
        for (unsigned i=0; i<p_mesh2->GetNumAllNodes(); i++)
        {
            //std::cout << " i " << i << "forces x " << forces_on_each_node[i][0] << ", y " << forces_on_each_node[i][1] << "\n" << std::endl;
            is_a_ghost_node = false;
            for (unsigned j=0; j<ghost_node_indices2.size(); j++)
            {
                if(ghost_node_indices2[j]==i)
                {
                    is_a_ghost_node = true;
                }
            }
            if(!is_a_ghost_node)
            {
                TS_ASSERT_DELTA(forces_on_each_node[i][0], 0.0, 1e-4);
                TS_ASSERT_DELTA(forces_on_each_node[i][1], 0.0, 1e-4);
            }
        }
        
        // Move a node along the x-axis and calculate the force exerted on a neighbour
        c_vector<double,2> old_point = p_mesh2->GetNode(59)->rGetLocation();
        Point<2> new_point;
        new_point.rGetLocation()[0] = old_point[0]+0.5;
        new_point.rGetLocation()[1] = old_point[1] ;
        p_mesh2->SetNode(59, new_point, false);
        forces_on_each_node = simulator3.CalculateForcesOnEachNode();
        TS_ASSERT_DELTA(forces_on_each_node[60][0], 0.5*p_params->GetMeinekeLambda(), 1e-4);
        TS_ASSERT_DELTA(forces_on_each_node[60][1], 0.0, 1e-4);
  
        c_vector<double,2> force_on_spring ; // between nodes 59 and 60
        
        // Find one of the elements that nodes 59 and 60 live on
        Point<2> new_point2;
        new_point2.rGetLocation()[0] = new_point[0]+0.01;
        new_point2.rGetLocation()[1] = new_point[1] + 0.01 ;
        
        unsigned elem_index = p_mesh2->GetContainingElementIndex(new_point2,false);
        Element<2,2>* p_element = p_mesh2->GetElement(elem_index);
        
        force_on_spring = simulator3.CalculateForceInThisSpring(p_element,1,0);
  
        TS_ASSERT_DELTA(force_on_spring[0], 0.5*p_params->GetMeinekeLambda(), 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);   
        
        Point<2> point_of_node60 = p_mesh2->GetNode(60)->rGetLocation();
        
        simulator3.SetDt(0.01);
        simulator3.UpdateNodePositions(forces_on_each_node);
                
        TS_ASSERT_DELTA(p_mesh2->GetNode(60)->rGetLocation()[0],point_of_node60.rGetLocation()[0]+force_on_spring[0]*0.01, 1e-4);
        TS_ASSERT_DELTA(p_mesh2->GetNode(60)->rGetLocation()[1],point_of_node60.rGetLocation()[1], 1e-4);
        
        
        
        std::vector<MeinekeCryptCell> cells3;
        simulator3.UpdateCellTypes();
        cells3 = simulator3.GetCells();
        
        for (unsigned i=0; i<num_cells2; i++)
        {
            CryptCellType cell_type;
            cell_type = cells3[i].GetCellType();
            TS_ASSERT_EQUALS(cell_type,TRANSIT);
            
        }
//      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//      THESE SHOULD BE DIFFERENTIATED NOT TRANSIT
//      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



//        CalculateForceInThisBoundarySpring



        SimulationTime::Destroy();
    }
        
};

#endif /*TESTCRYPTSIMULATION2DPERIODIC_HPP_*/
