#ifndef TESTCRYPTSIMULATION2D_HPP_
#define TESTCRYPTSIMULATION2D_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>
#include "OutputFileHandler.hpp"
#include "CryptSimulation2D.hpp"
#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"


class TestCryptSimulation2D : public CxxTest::TestSuite
{
    void Make2dCryptMesh(std::string meshFilename, unsigned numNodesAlongWidth, unsigned numNodesAlongLength, double width, double length, double x0=0.0, double y0=0.0)
    {
        OutputFileHandler output_file_handler("CryptMesh");
        out_stream p_node_file = output_file_handler.OpenOutputFile(meshFilename+".node");
        (*p_node_file) << std::scientific;
        
        out_stream p_elem_file = output_file_handler.OpenOutputFile(meshFilename+".ele");
        (*p_elem_file) << std::scientific;
        
        unsigned num_nodes            = numNodesAlongWidth*numNodesAlongLength;
        unsigned num_elem_along_width = numNodesAlongWidth-1;
        unsigned num_elem_along_length = numNodesAlongLength-1;
        unsigned num_elem             = 2*num_elem_along_width*num_elem_along_length;
        unsigned num_edges            = 3*num_elem_along_width*num_elem_along_length + num_elem_along_width + num_elem_along_length;
        
        (*p_node_file) << num_nodes << "\t2\t0\t1" << std::endl;
        unsigned node = 0;
        for (unsigned i = 0; i < numNodesAlongLength; i++)
        {
            for (unsigned j = 0; j < numNodesAlongWidth; j++)
            {
                int b = 0;
                if ((i==0) || (i==numNodesAlongLength-1) || (j==0) || (j==numNodesAlongWidth-1))
                {
                    b = 1;
                }
                double x = x0 + width*((double)j + 0.25*(1+ pow(-1,i+1)))/(num_elem_along_width) ;
                
                double y = y0 + length*(sqrt(3)/2)*(double)i/(num_elem_along_length);
                
                (*p_node_file) << node++ << "\t" << x << "\t" << y << "\t" << b << std::endl;
            }
        }
        p_node_file->close();
        
        out_stream p_edge_file = output_file_handler.OpenOutputFile(meshFilename+".edge");
        (*p_node_file) << std::scientific;
        
        (*p_elem_file) << num_elem << "\t3\t0" << std::endl;
        (*p_edge_file) << num_edges << "\t3\t0\t1" << std::endl;
        
        unsigned elem = 0;
        unsigned edge = 0;
        for (unsigned i = 0; i < num_elem_along_length; i++)
        {
            for (unsigned j = 0; j < num_elem_along_width; j++)
            {
                int node0 =     i*numNodesAlongWidth + j;
                int node1 =     i*numNodesAlongWidth + j+1;
                int node2 = (i+1)*numNodesAlongWidth + j;
                if (i%2 != 0)
                {
                    node2 = node2 + 1;
                }
                
                (*p_elem_file) << elem++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;
                
                int horizontal_edge_is_boundary_edge = 0;
                int vertical_edge_is_boundary_edge = 0;
                if (i==0)
                {
                    horizontal_edge_is_boundary_edge = 1;
                }
                if (j==0)
                {
                    vertical_edge_is_boundary_edge = 1;
                }
                
                (*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 <<  "\t" << horizontal_edge_is_boundary_edge << std::endl;
                (*p_edge_file) << edge++ << "\t" << node1 << "\t" << node2 <<  "\t" << 0 << std::endl;
                (*p_edge_file) << edge++ << "\t" << node2 << "\t" << node0 <<  "\t" << vertical_edge_is_boundary_edge << std::endl;
                
                node0 = i*numNodesAlongWidth + j + 1;
                
                if (i%2 != 0)
                {
                    node0 = node0 - 1;
                }
                node1 = (i+1)*numNodesAlongWidth + j+1;
                node2 = (i+1)*numNodesAlongWidth + j;
                
                (*p_elem_file) << elem++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;
            }
        }
        
        for (unsigned i = 0; i < num_elem_along_length; i++)
        {
            int node0 = (i+1)*numNodesAlongWidth-1;
            int node1 = (i+2)*numNodesAlongWidth-1;
            (*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << 1 << std::endl;
        }
        
        for (unsigned j = 0; j < num_elem_along_width; j++)
        {
            int node0 =  numNodesAlongWidth*(numNodesAlongLength-1) + j;
            int node1 =  numNodesAlongWidth*(numNodesAlongLength-1) + j+1;
            (*p_edge_file) << edge++ << "\t" << node1 << "\t" << node0 << "\t" << 1 << std::endl;
        }
        
        p_elem_file->close();
        p_edge_file->close();
    }
    
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
        srandom(0);
        
        double crypt_length = 10;
        double crypt_width = 10;
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
       
        CryptSimulation2D simulator(mesh);
        simulator.SetOutputDirectory("Crypt2DSprings");

        simulator.SetEndTime(24.0);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(400);

        simulator.SetReMeshRule(false);
        simulator.Solve();
              
        CheckAgainstPreviousRun("Crypt2DSprings", 400u, 400u);
    }
    
    // note - there is no remeshing here so it will crash if run for too long
    void Test2DSpringsWithCells() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        srandom(0);
        double crypt_length = 10;
        double crypt_width = 10;
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            
            double y = mesh.GetNodeAt(i)->GetPoint().rGetLocation()[1];
            if (y == 0.0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 5.0)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetTransitCellCycleTime(); //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0; //hours
            }
            MeinekeCryptCell cell(cell_type, 0.0, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        

        CryptSimulation2D simulator(mesh, cells);
        simulator.SetOutputDirectory("Crypt2DSpringsWithCells");

        simulator.SetEndTime(0.45*24.0);

        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);

        simulator.SetReMeshRule(false);
                
        // throws anything because not working at the moment
        simulator.Solve();
        CheckAgainstPreviousRun("Crypt2DSpringsWithCells", 400u, 800u);
    }
    
    
    // not being run because takes a few minutes to run
    void DO_NOT________TestWithBirthOnHoneycombMesh() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        srandom(0);  // this is BAD, mkay, no way?
        
        double crypt_length = 10.0;
        double crypt_width = 5.0;
        
        Make2dCryptMesh("2D_crypt_mesh", 8, 15, crypt_width+2, crypt_length+4, -1.0, -sqrt(3)/2.0);
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<2,2> mesh_reader(testoutput_dir+"/CryptMesh/2D_crypt_mesh");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        std::vector<int> ghost_node_indices;
        for (int i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNodeAt(i)->GetPoint().rGetLocation()[0];
            double y = mesh.GetNodeAt(i)->GetPoint().rGetLocation()[1];
            if ((x<0)||(x>crypt_width)||(y>crypt_length)||(y<0))
            {
               ghost_node_indices.push_back(i);
            }
        }
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            double y = mesh.GetNodeAt(i)->GetPoint().rGetLocation()[1];
            if (y == 0.0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 5.0)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetTransitCellCycleTime(); //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0; //hours
            }
            
            MeinekeCryptCell cell(cell_type, 0.0, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        CryptSimulation2D simulator(mesh, cells);
        simulator.SetOutputDirectory("Crypt2DHoneycombMesh");
        simulator.SetEndTime(24.0);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        
        simulator.SetGhostNodes(ghost_node_indices);

        simulator.Solve();
        CheckAgainstPreviousRun("Crypt2DHoneycombMesh", 400u, 800u);
    }
    
    
    void Test2DSpringsFixedBoundaries() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        srandom(0);
        double crypt_length = 10;
        double crypt_width = 10;
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
          
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            
            double y = mesh.GetNodeAt(i)->GetPoint().rGetLocation()[1];
            if (y == 0.0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 6.5)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 8)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetTransitCellCycleTime(); //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0; //hours
            }
            
            MeinekeCryptCell cell(cell_type, 0.0, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        CryptSimulation2D simulator(mesh,cells);
        simulator.SetOutputDirectory("Crypt2DSpringsFixedBoundaries");
        simulator.SetEndTime(0.5); //hours
        simulator.SetMaxCells(800);
        simulator.SetMaxElements(800);
        simulator.SetFixedBoundaries();
        
        simulator.Solve();
        CheckAgainstPreviousRun("Crypt2DSpringsFixedBoundaries", 400u, 800u);
    }
    
    void TestCalculateCryptBoundaries()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(0.0,-2.0,0.0) ;
        //Create Vector of ghost nodes
        std::vector<int> ghost_node_indices;
        for (int i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNodeAt(i)->GetPoint().rGetLocation()[0];
            double y = mesh.GetNodeAt(i)->GetPoint().rGetLocation()[1];
            if ((x<2.0)||(x>8.0)||(y>6.0)||(y<0.0))
            {
               ghost_node_indices.push_back(i);
            }
        }
        
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->SetCryptLength(6.0);
        
        CryptSimulation2D simulator(mesh);
        simulator.SetGhostNodes(ghost_node_indices);
        simulator.CalculateCryptBoundary();
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
        
        //* Must uncomment and include this test
        
        //TS_ASSERT_EQUALS(actual_left_boundary_nodes.size(),calculated_left_boundary_nodes.size());
        //TS_ASSERT_EQUALS(actual_right_boundary_nodes.size(),calculated_right_boundary_nodes.size());
        
        
        //for(unsigned i=0; i<calculated_left_boundary_nodes.size(); i++)
        //{
        //    std::cout<< "calculated_left_boundary_nodes "<< calculated_left_boundary_nodes[i] <<"\n" << std::flush;
            
            //TS_ASSERT_EQUALS(actual_left_boundary_nodes[i],calculated_left_boundary_nodes[i]);
            //TS_ASSERT_EQUALS(actual_right_boundary_nodes[i],calculated_right_boundary_nodes[i]);
       // }
        
    }
};

#endif /*TESTCRYPTSIMULATION2D_HPP_*/

