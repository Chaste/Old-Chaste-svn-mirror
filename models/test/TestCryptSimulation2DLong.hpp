#ifndef TESTCRYPTSIMULATION2DLONG_HPP_
#define TESTCRYPTSIMULATION2DLONG_HPP_

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


class TestCryptSimulation2DLong : public CxxTest::TestSuite
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

    void TestWithBirthOnHoneycombMesh() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator random_num_gen;
        
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
                birth_time = -random_num_gen.ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
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
    
    
    //////////////////////////////////////////////////////////////////
    // starting with a small mesh with 1 stem cell and the rest
    // differentiated, check the number of cells at the end of the 
    // simulation is as expected.
    //////////////////////////////////////////////////////////////////
    void Test2DCorrectCellNumbers() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator random_num_gen;
        
        // check the stem cell cycle time is still 24 hrs, otherwise
        // this test might not pass
        TS_ASSERT_DELTA(p_params->GetStemCellCycleTime(), 24, 1e-12);  
        TS_ASSERT_DELTA(p_params->GetTransitCellCycleTime(), 12, 1e-12);  
                
        double crypt_length = 4.0;
        double crypt_width = 5.0;
        
        Make2dCryptMesh("2D_crypt_mesh", 8, 9, crypt_width+2, crypt_length+4, -1.0, -sqrt(3)/2.0);
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
            
            if (i==11) // middle of bottom row of cells
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
            
            MeinekeCryptCell cell(cell_type, 0.0, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        CryptSimulation2D simulator(mesh,cells);
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
        
        for (int index = 0; index<mesh.GetNumAllNodes(); index++)
        {
            if(!is_ghost_node[index])   //!mesh.GetNodeAt(index)->IsDeleted())
            {
                CryptCellType type = cells_after_simulation[index].GetCellType();
                
                if(type==STEM)
                {
                    num_stem++;
                }
                else if(type==TRANSIT)
                {
                    num_transit++;
                }
                else if(type==DIFFERENTIATED)
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
    } 
    
    
};


#endif /*TESTCRYPTSIMULATION2DLONG_HPP_*/
