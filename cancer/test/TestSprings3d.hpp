#ifndef TESTSPRINGS3D_HPP_
#define TESTSPRINGS3D_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "TrianglesMeshWriter.cpp"
#include <cmath>
#include <ctime>
#include <vector>
#include "OutputFileHandler.hpp"
#include "CancerParameters.hpp"
#include "NodeMap.hpp"
#include "TissueCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "ColumnDataReader.hpp"
#include "SimulationTime.hpp"

class TestSprings3d : public CxxTest::TestSuite
{
   ConformingTetrahedralMesh<3,3> Make3dMesh(unsigned width=3, unsigned height=3, unsigned depth=3)
   {
       ConformingTetrahedralMesh<3,3> mesh;
       mesh.ConstructCuboid(width,height,depth,true);
       TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMesh");
       mesh_writer.WriteFilesUsingMesh(mesh);

       return mesh;
   }
   
   void CheckAgainstPreviousRun3D(std::string resultDirectory,std::string resultSet, unsigned maxCells, unsigned maxElements)
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
            std::stringstream cell_z_position_var_name;
            cell_type_var_name << "cell_type_" << cell;
            cell_x_position_var_name << "cell_x_position_" << cell;
            cell_y_position_var_name << "cell_y_position_" << cell;
            cell_z_position_var_name << "cell_z_position_" << cell;
            
            // Vector of Cell Types
            std::vector<double> expected_cell_types = expected_node_results.GetValues(cell_type_var_name.str());
            std::vector<double> computed_cell_types = computed_node_results.GetValues(cell_type_var_name.str());
            
            //Vector of Cell Positions
            std::vector<double> expected_cell_x_positions = expected_node_results.GetValues(cell_x_position_var_name.str());
            std::vector<double> computed_cell_x_positions = computed_node_results.GetValues(cell_x_position_var_name.str());
            
            
            std::vector<double> expected_cell_y_positions = expected_node_results.GetValues(cell_y_position_var_name.str());
            std::vector<double> computed_cell_y_positions = computed_node_results.GetValues(cell_y_position_var_name.str());
            
            std::vector<double> expected_cell_z_positions = expected_node_results.GetValues(cell_z_position_var_name.str());
            std::vector<double> computed_cell_z_positions = computed_node_results.GetValues(cell_z_position_var_name.str());
            
            //Comparing expected and computed vector length
            TS_ASSERT_EQUALS(expected_cell_types.size(), computed_cell_types.size());
            TS_ASSERT_EQUALS(expected_cell_x_positions.size(), computed_cell_x_positions.size());
            TS_ASSERT_EQUALS(expected_cell_y_positions.size(), computed_cell_y_positions.size());
            TS_ASSERT_EQUALS(expected_cell_z_positions.size(), computed_cell_z_positions.size());
            
            //Walkthrough of the expected and computed vectors
            for (unsigned time_step = 0; time_step < expected_cell_types.size(); time_step++)
            {
                TS_ASSERT_EQUALS(expected_cell_types[time_step], computed_cell_types[time_step]);
                TS_ASSERT_DELTA(expected_cell_x_positions[time_step], computed_cell_x_positions[time_step],1e-6);
                TS_ASSERT_DELTA(expected_cell_y_positions[time_step], computed_cell_y_positions[time_step],1e-6);
                TS_ASSERT_DELTA(expected_cell_z_positions[time_step], computed_cell_z_positions[time_step],1e-6);
            }
        }
        
        for (unsigned element=0; element<maxElements; element++)
        {
            std::stringstream nodeA_var_name;
            std::stringstream nodeB_var_name;
            std::stringstream nodeC_var_name;
            std::stringstream nodeD_var_name;
            
            nodeA_var_name << "nodeA_" << element;
            nodeB_var_name << "nodeB_" << element;
            nodeC_var_name << "nodeC_" << element;
            nodeD_var_name << "nodeD_" << element;
            
            // Vector of Node A indexes
            std::vector<double> expected_NodeA_numbers = expected_element_results.GetValues(nodeA_var_name.str());
            std::vector<double> computed_NodeA_numbers = computed_element_results.GetValues(nodeA_var_name.str());
            
            // Vector of Node B indexes
            std::vector<double> expected_NodeB_numbers = expected_element_results.GetValues(nodeB_var_name.str());
            std::vector<double> computed_NodeB_numbers = computed_element_results.GetValues(nodeB_var_name.str());
            
            // Vector of Node C indexes
            std::vector<double> expected_NodeC_numbers = expected_element_results.GetValues(nodeC_var_name.str());
            std::vector<double> computed_NodeC_numbers = computed_element_results.GetValues(nodeC_var_name.str());
            
            // Vector of Node D indexes
            std::vector<double> expected_NodeD_numbers = expected_element_results.GetValues(nodeD_var_name.str());
            std::vector<double> computed_NodeD_numbers = computed_element_results.GetValues(nodeD_var_name.str());
            
            TS_ASSERT_EQUALS(expected_NodeA_numbers.size(), computed_NodeA_numbers.size());
            TS_ASSERT_EQUALS(expected_NodeB_numbers.size(), computed_NodeB_numbers.size());
            TS_ASSERT_EQUALS(expected_NodeC_numbers.size(), computed_NodeC_numbers.size());
            TS_ASSERT_EQUALS(expected_NodeD_numbers.size(), computed_NodeD_numbers.size());
            
            for (unsigned time_step = 0; time_step < expected_NodeA_numbers.size(); time_step++)
            {
                TS_ASSERT_EQUALS(expected_NodeA_numbers[time_step], computed_NodeA_numbers[time_step]);
                TS_ASSERT_EQUALS(expected_NodeB_numbers[time_step], computed_NodeB_numbers[time_step]);
                TS_ASSERT_EQUALS(expected_NodeC_numbers[time_step], computed_NodeC_numbers[time_step]);
                TS_ASSERT_EQUALS(expected_NodeD_numbers[time_step], computed_NodeD_numbers[time_step]);
            }
        }
    }
    

   

public:

//   void Test3dSpringMesh() throw (Exception)
//   {
//
//        ConformingTetrahedralMesh<3,3> mesh = Make3dMesh();
//        
////        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
////        ConformingTetrahedralMesh<3,3> mesh;
////        mesh.ConstructFromMeshReader(mesh_reader);
////        
//        
//        {
//            TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMeshStart");
//            mesh_writer.WriteFilesUsingMesh(mesh);
//
//        }
//        
//        
//        
//        for (unsigned i=0 ; i<2 ; i++)
//        {
////            std::cout << "step = " << i << std::endl;
//            std::vector<std::vector<double> > forces = CalculateForcesOnEachNode(mesh);
//            UpdateNodePositions(forces,mesh);
//            
////            NodeMap map(mesh.GetNumAllNodes());
////            mesh.ReMesh(map);
//
////            std::ostringstream time_stamp;
////            time_stamp << i;
////            TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMesh"+time_stamp.str());
////            mesh_writer.WriteFilesUsingMesh(mesh);
//        }
//        
//        
//        
//        TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMeshEnd");
//        mesh_writer.WriteFilesUsingMesh(mesh);
//        
//        /*
//         * Use the command:
//         * 
//         * tetview 3dSpringMesh
//         * 
//         * to visualize final output.
//         */
//        
//   }

    void TestOne3dElement() throw (Exception)
    {
        CancerParameters::Instance()->Reset();

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element");
        
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TrianglesMeshWriter<3,3> mesh_writer1("","3dSpringTetrahedronMeshStart");
        mesh_writer1.WriteFilesUsingMesh(mesh);
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        std::vector<TissueCell> cells;
        CryptCellType cell_type;
        cell_type = STEM;
        
        TissueCell cell(cell_type, HEALTHY, 0u, new FixedCellCycleModel());
        for(unsigned i=0; i< mesh.GetNumNodes(); i++)
        {
            cell.SetNodeIndex(i);
            cell.SetBirthTime(-50.0);
            cells.push_back(cell);
        }
        
        Crypt<3> crypt(mesh,cells);
        TissueSimulation<3> simulator(crypt);
        simulator.SetMaxCells(40);
        simulator.SetMaxElements(40);

        // Test forces on springs
        unsigned nodeA = 0, nodeB = 1 ;
        Element<3,3>* p_element = mesh.GetElement(0);
        c_vector<double, 3> force = simulator.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(nodeA),p_element->GetNodeGlobalIndex(nodeB));
        for(unsigned i=0; i < 3;i++)
        {
            TS_ASSERT_DELTA(force[i],0.0,1e-6);
        }
        
        //  Test forces on nodes        
        for (unsigned i=0 ; i<1 ; i++)
        {
            std::vector<c_vector<double,3> > velocities = simulator.CalculateVelocitiesOfEachNode();
            simulator.UpdateNodePositions(velocities);
            
            for (unsigned j=0; j<4; j++)
            {
                for(unsigned k=0;k<3;k++)
                {
                    TS_ASSERT_DELTA(velocities[j](k),0.0,1e-6);
                }
            }
            
        }
        TrianglesMeshWriter<3,3> mesh_writer("","3dSpringTetrahedronMeshEnd");
        mesh_writer.WriteFilesUsingMesh(mesh);    
        
        /*
         ************************************************************************
         ************************************************************************ 
         *  Scale entire mesh and check that forces are correctly calculated
         ************************************************************************
         ************************************************************************ 
         */
        CancerParameters *p_params = CancerParameters::Instance();
        double scale_factor = 1.5;
        
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            c_vector<double,3> old_point = mesh.GetNode(i)->rGetLocation();
            ChastePoint<3> new_point;
            new_point.rGetLocation()[0] = scale_factor*old_point[0];
            new_point.rGetLocation()[1] = scale_factor*old_point[1];
            new_point.rGetLocation()[2] = scale_factor*old_point[2];
            mesh.SetNode(i, new_point, false);            
        }
               
        std::vector<c_vector<double,3> > new_velocities = simulator.CalculateVelocitiesOfEachNode();
        simulator.UpdateNodePositions(new_velocities);
        for (unsigned j=0; j<4; j++)
        {
            for(unsigned k=0;k<3;k++)
            {
                TS_ASSERT_DELTA(fabs(new_velocities[j](k)),p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(scale_factor-1)*sqrt(2),1e-6);
            }
        }
        
        /*
         ************************************************************************
         ************************************************************************ 
         *  Move one node and check that forces are correctly calculated
         ************************************************************************
         ************************************************************************ 
         */
        ConformingTetrahedralMesh<3,3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);
        Crypt<3> crypt2(mesh2,cells);
        TissueSimulation<3> simulator2(crypt2);
        
        simulator2.SetMaxCells(40);
        simulator2.SetMaxElements(40);
        
        c_vector<double,3> old_point = mesh2.GetNode(0)->rGetLocation();
        ChastePoint<3> new_point;
        new_point.rGetLocation()[0] = 0.0;
        new_point.rGetLocation()[1] = 0.0;
        new_point.rGetLocation()[2] = 0.0;
        mesh2.SetNode(0, new_point, false);   
          
        /*
         ************************************************************************
         ************************************************************************ 
         *  Test forces on springs
         ************************************************************************
         ************************************************************************ 
         */
        unsigned nodeA2 = 0, nodeB2 = 1 ;
        Element<3,3>* p_element2 = mesh2.GetElement(0);
        c_vector<double,3> force2 = simulator2.CalculateForceBetweenNodes(p_element2->GetNodeGlobalIndex(nodeA2),p_element2->GetNodeGlobalIndex(nodeB2));
        
        for(unsigned i=0; i < 3;i++)
        {
            TS_ASSERT_DELTA(fabs(force2[i]),p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }
        
        std::vector<c_vector<double,3> > new_velocities2 = simulator2.CalculateVelocitiesOfEachNode();
        
        for(unsigned k=0;k<3;k++)
        {
            TS_ASSERT_DELTA(new_velocities2[0](k),p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }
        simulator2.UpdateNodePositions(new_velocities2);
        TrianglesMeshWriter<3,3> mesh_writer2("","3dSpringTetrahedronMeshEnd2");
        mesh_writer2.WriteFilesUsingMesh(mesh2);  
        
        
        /*
         ************************************************************************
         ************************************************************************ 
         *  Test Cell Birth
         ************************************************************************
         ************************************************************************ 
         */                
        ConformingTetrahedralMesh<3,3> mesh3;
        mesh3.ConstructFromMeshReader(mesh_reader);
        
        std::vector<TissueCell> cells2;
        
        for(unsigned i=0; i< mesh.GetNumNodes()-1; i++)
        {
            cell.SetNodeIndex(i);
            cell.SetBirthTime(0.0);
            cells2.push_back(cell);
        }
        // Setting last cell to undergo cell birth.
        cell.SetNodeIndex(mesh.GetNumNodes()-1);
        cell.SetBirthTime(-50.0);
        cells2.push_back(cell);
        
        Crypt<3> crypt3(mesh3,cells2);
        TissueSimulation<3> simulator3(crypt3);
       
        TrianglesMeshWriter<3,3> mesh_writer3("Test3DCellBirth","StartMesh");
        mesh_writer3.WriteFilesUsingMesh(mesh3);
        
        simulator3.SetMaxCells(10);
        simulator3.SetMaxElements(25);
        simulator3.SetOutputDirectory("Test3DCellBirth");
        
        // Set to re-mesh
        simulator3.SetReMeshRule(true);
        simulator3.SetEndTime(1.0);
        
        simulator3.Solve();
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        
        TrianglesMeshWriter<3,3> mesh_writer4("Test3DCellBirth","EndMesh",false); 
        mesh_writer4.WriteFilesUsingMesh(mesh3);
        
        TS_ASSERT_EQUALS(mesh3.GetNumNodes(),5u);
        TS_ASSERT_EQUALS(mesh3.GetNumElements(),3u);
    } 
    

    void TestPrivateFunctionsOfSpheroidSimulation3D() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        RandomNumberGenerator::Instance();
        
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            cell_type = STEM;
            generation = 0;
            TissueCell cell(cell_type, HEALTHY, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            if ( i == 50u)
            {
                cell.SetBirthTime(-50.0 );
            }
            
            cells.push_back(cell);
        }
        
        Crypt<3> crypt(mesh,cells);
        TissueSimulation<3> simulator(crypt);
        
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(2400);
        
        unsigned num_births = simulator.DoCellBirth();
                                                                
        TS_ASSERT_EQUALS(num_births, 1u);
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
       
 
    void TestSolveMethodSpheroidSimulation3D() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        RandomNumberGenerator *p_random_num_gen=RandomNumberGenerator::Instance();
                
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<3,3> mesh_writer("TestSolveMethodSpheroidSimulation3DMesh","StartMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);

        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            cell_type = STEM;
            generation = 0;
            TissueCell cell(cell_type, HEALTHY, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);    
            cell.SetBirthTime(p_random_num_gen->ranf()*p_params->GetStemCellCycleTime());            
            cells.push_back(cell);
        } 
        
        Crypt<3> crypt(mesh,cells);
        TissueSimulation<3> simulator(crypt);
        
        simulator.SetMaxCells(1000);
        simulator.SetMaxElements(2500);
        simulator.SetOutputDirectory("TestSolveMethodSpheroidSimulation3D");
        
        // Set to re-mesh
        simulator.SetReMeshRule(true);
        simulator.SetEndTime(0.1);
        
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        
        TrianglesMeshWriter<3,3> mesh_writer2("TestSolveMethodSpheroidSimulation3DMesh","EndMesh",false); 
        mesh_writer2.WriteFilesUsingMesh(mesh);      
    }
    
 
    void TestGhostNodesSpheroidSimulation3D() throw (Exception)
    {
        double start_time = std::clock();
        
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        RandomNumberGenerator *p_random_num_gen=RandomNumberGenerator::Instance();
                       
        unsigned width = 3;
        unsigned height = 3;               
        unsigned depth = 3;      
                 
        ConformingTetrahedralMesh<3,3> mesh = Make3dMesh(width,height,depth);
        TrianglesMeshWriter<3,3> mesh_writer("TestGhostNodesSpheroidSimulation3D","StartMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);
                       
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
      
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<TissueCell> cells;
                  
        c_vector<double, 3> spheroid_centre;
        spheroid_centre[0] = 0.5*((double) width);
        spheroid_centre[1] = 0.5*((double) height);
        spheroid_centre[2] = 0.5*((double) depth);
        
        std::set<unsigned> ghost_node_indices;

        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            
            c_vector<double, 3> node_location = mesh.GetNode(i)->rGetLocation();
            
            unsigned min_spatial_dimension;
            if (width <= height && width <= depth)
            {
                min_spatial_dimension = width;
            }
            else
            {
                if (height <= depth)
                {
                    min_spatial_dimension = height;
                }
                else
                {
                    min_spatial_dimension = depth;
                }
            }
            if ( norm_2(node_location - spheroid_centre) > 0.5*sqrt(3)*1.01*((double) min_spatial_dimension)/3.0 )
            {
                ghost_node_indices.insert(i);
            }
                
            cell_type = STEM;
            generation = 0;
            TissueCell cell(cell_type, HEALTHY, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);    
            cell.SetBirthTime(-p_random_num_gen->ranf()*p_params->GetStemCellCycleTime());      
            cells.push_back(cell);
        } 
                    
        TS_ASSERT(ghost_node_indices.size() < num_cells);
        TS_ASSERT(ghost_node_indices.size() > 0)
        
        Crypt<3> crypt(mesh,cells);
        crypt.SetGhostNodes(ghost_node_indices);        

        TissueSimulation<3> simulator(crypt);

        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        simulator.SetOutputDirectory("TestGhostNodesSpheroidSimulation3D");
        
        // Set to re-mesh
        simulator.SetReMeshRule(true);
        simulator.SetEndTime(0.1);

        // TODO: make sure th right number of new cells are born                        
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        
        TrianglesMeshWriter<3,3> mesh_writer2("TestGhostNodesSpheroidSimulation3D","EndMesh",false); 
        mesh_writer2.WriteFilesUsingMesh(mesh);
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        
        double end_time = std::clock();
        double elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "Time of simulation " << elapsed_time << "\n" << std::flush; 
    }
};

#endif /*TESTSPRINGS3D_HPP_*/

