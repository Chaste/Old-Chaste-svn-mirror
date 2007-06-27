#ifndef TESTCRYPT_HPP_
#define TESTCRYPT_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <vector>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "OutputFileHandler.hpp"
#include "MeinekeCryptCell.hpp"
#include "CancerParameters.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"
#include "Crypt.cpp"

class TestCrypt : public CxxTest::TestSuite
{
    
private: 
    // test construction (without ghost nodes), accessors and iterator
    template<unsigned DIM>
    void TestSimpleCrypt(std::string meshFilename)
    {
        // set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // create a simple mesh
        TrianglesMeshReader<DIM,DIM> mesh_reader(meshFilename);
        ConformingTetrahedralMesh<DIM,DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            double birth_time = 0.0-i;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        // create the crypt
        Crypt<DIM> crypt(mesh, cells);
        
        TS_ASSERT_EQUALS(crypt.rGetMesh().GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(crypt.rGetCells().size(), cells.size());
        
        unsigned counter = 0;
        for (typename Crypt<DIM>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            // test operator* and that cells are in sync
            TS_ASSERT_EQUALS((*cell_iter).GetNodeIndex(), counter);
 
            // test operator-> and that cells are in sync
            TS_ASSERT_DELTA(cell_iter->GetAge(), (double)counter, 1e-12);
 
            // test GetNode on the iterator
            TS_ASSERT_EQUALS(cell_iter.GetNode()->GetIndex(), mesh.GetNode(counter)->GetIndex());
 
            // test iter.GetNode()->GetIndex() is consistent with cell.GetNodeIndex()
            TS_ASSERT_EQUALS((*cell_iter).GetNodeIndex(), cell_iter.GetNode()->GetIndex());
 
            // test rGetLocation on the iterator
            for(unsigned space_index=0; space_index<DIM; space_index++)
            {
                TS_ASSERT_EQUALS(cell_iter.rGetLocation()[space_index], 
                                 mesh.GetNode(counter)->rGetLocation()[space_index]);
            }
  
            counter++;
        }
        
        TS_ASSERT_EQUALS(counter, crypt.GetNumRealCells());
        
        SimulationTime::Destroy();
    }


public:

    // test construction, accessors and Iterator
    void TestSimpleCrypt1d2d3d() throw(Exception)
    {
        TestSimpleCrypt<1>("mesh/test/data/1D_0_to_1_10_elements");
        TestSimpleCrypt<2>("mesh/test/data/square_4_elements");
        TestSimpleCrypt<3>("mesh/test/data/cube_136_elements");
    }
    
    void TestValidate()
    {
    	// set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            double birth_time = 0.0-i;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        cells[0].SetNodeIndex(1);

		// fails as no cell or ghost correponding to node 0        
        TS_ASSERT_THROWS_ANYTHING(Crypt<2> crypt2(mesh, cells));

        SimulationTime::Destroy();
    }
    
    
    // test with ghost nodes, incl the Iterator doesn't loop over ghost nodes
    void TestCryptWithGhostNodes() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);        
        
        unsigned num_cells_depth = 11;
        unsigned num_cells_width = 6;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);

        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(0);
            cells.push_back(cell);
        }
        
        // create a crypt, with no ghost nodes at the moment
        Crypt<2> crypt(*p_mesh,cells);

        // iterator should loop over all nodes
        unsigned counter = 0;
        for (Crypt<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            counter++;
        }
        TS_ASSERT_EQUALS(counter, p_mesh->GetNumNodes());
        
        // get ghost node indices from the generator and set up the ghost nodes
        // vector of bools
        std::vector<bool> is_ghost_node(p_mesh->GetNumNodes(),false);
        for(std::set<unsigned>::iterator iter = ghost_node_indices.begin();
            iter!=ghost_node_indices.end();
            iter++)
        {
            is_ghost_node[*iter] = true;
        }
        
        // set ghost nodes
        crypt.SetGhostNodes(is_ghost_node);
        
        TS_ASSERT_EQUALS(crypt.rGetGhostNodes(), is_ghost_node);
        
        // check the iterator doesn't loop over ghost nodes
        counter = 0;
        for (Crypt<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            unsigned node_index = cell_iter->GetNodeIndex();
            TS_ASSERT_EQUALS(is_ghost_node[node_index], false);
            counter++;
        }
        
        TS_ASSERT_EQUALS(counter, crypt.GetNumRealCells());
        
        // check counter = num_nodes - num_ghost_nodes
        TS_ASSERT_EQUALS(counter + ghost_node_indices.size(), p_mesh->GetNumNodes());
        
        SimulationTime::Destroy();
    }
    
    void TestMoveCellAndAddCell()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);        
        
        // create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(0);
            cells.push_back(cell);
        }
        // create a crypt, with no ghost nodes at the moment
        Crypt<2> crypt(mesh, cells);

        //////////////////
        // test move cell
        //////////////////
        
        // move node 0 by a small amount
        Crypt<2>::Iterator cell_iter = crypt.Begin();
        c_vector<double,2> new_location = cell_iter.rGetLocation();
        new_location[0] += 1e-2;
        new_location[1] += 1e-2;
        Point<2> new_location_point(new_location);
        crypt.MoveCell(cell_iter, new_location_point);
        
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[0], new_location[0], 1e-12);
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[1], new_location[1], 1e-12);

        //////////////////
        // test add cell
        //////////////////
        unsigned old_num_nodes = mesh.GetNumNodes();
        unsigned old_num_cells = crypt.rGetCells().size();

        // create a new cell, DON'T set the node index, set birth time=-1
        MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
        cell.SetBirthTime(-1);
        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 2;
        new_cell_location[1] = 2;
        
        crypt.AddCell(cell,new_cell_location); 

        // crypt should have updated mesh and cells
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), old_num_nodes+1);
        TS_ASSERT_EQUALS(crypt.rGetCells().size(), old_num_cells+1);
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), old_num_nodes+1);

        // same test via crypt class
        TS_ASSERT_EQUALS(crypt.rGetMesh().GetNumNodes(), old_num_nodes+1);
        TS_ASSERT_EQUALS(crypt.rGetCells().size(), old_num_cells+1);
        
        // check the location of the new node
        TS_ASSERT_DELTA(mesh.GetNode(old_num_nodes)->rGetLocation()[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(mesh.GetNode(old_num_nodes)->rGetLocation()[1], 2.0, 1e-12);

        // check the index of the new cell
        MeinekeCryptCell& new_cell = crypt.rGetCells().back();
        TS_ASSERT_EQUALS(new_cell.GetNodeIndex(), old_num_nodes);
        
        SimulationTime::Destroy();
    }
    
    void TestRemoveDeadCellsAndReMesh()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);        
        
        // create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(0);
            cells.push_back(cell);
        }
        
        cells[27].StartApoptosis();
        
        // create a crypt, with some random ghost nodes
        Crypt<2> crypt(mesh,cells);

        std::vector<bool> is_ghost_node(mesh.GetNumNodes(), false);
        for(unsigned i=0; i<10; i++)
        {
            is_ghost_node[i] = true;
        }
        is_ghost_node[80]=true;
        
        crypt.SetGhostNodes(is_ghost_node);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(crypt.rGetCells().size(), 81u);

        // num real cells should be num_nodes (81) - num_ghosts (11) = 70
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 70u);

        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);
        p_simulation_time->IncrementTimeOneStep();
        unsigned num_removed = crypt.RemoveDeadCells();
        
        TS_ASSERT_EQUALS(num_removed, 1u);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 80u);
        TS_ASSERT_EQUALS(crypt.rGetCells().size(), 80u);
        TS_ASSERT_DIFFERS(crypt.rGetCells().size(), cells.size()); // Crypt now copies cells
        
        // num real cells should be num_nodes (81) - num_ghosts (11) - One deleted node = 69
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 69u);
        
        TS_ASSERT_EQUALS(crypt.rGetGhostNodes().size(), mesh.GetNumAllNodes()); 

        crypt.ReMesh();

        // num real cells should be num_new_nodes (80) - num_ghosts (11)
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 69u);

        // test size of ghost nodes vector is correct - mesh.GetNumNodes() ?
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh.GetNumAllNodes()); 
        TS_ASSERT_EQUALS(crypt.rGetGhostNodes().size(), mesh.GetNumNodes()); 
        
        // nodes 0-9 should not been renumbered so are still ghost nodes.
        // the ghost node at node 80 is now at 79 as node 27 was deleted..
        for(unsigned i=0; i<crypt.rGetGhostNodes().size(); i++)
        {
            // true (ie should be a ghost) if i<10 or i==79, else false
            TS_ASSERT_EQUALS(crypt.rGetGhostNodes()[i], ((i<10)||(i==79))); 
        }
        
        // finally, check the cells node indices have updated..
        
        SimulationTime::Destroy(); 
    }
    
    
    void TestAddAndRemoveAndAddWithOutRemesh()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);        
        
        // create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(0);
            cells.push_back(cell);
        }
        
        cells[27].StartApoptosis();
        
        // create a crypt, with some random ghost nodes
        Crypt<2> crypt(mesh,cells);

        std::vector<bool> is_ghost_node(mesh.GetNumNodes(), false);
        for(unsigned i=0; i<10; i++)
        {
            is_ghost_node[i] = true;
        }
        is_ghost_node[80]=true;
        
        crypt.SetGhostNodes(is_ghost_node);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(crypt.rGetCells().size(), 81u);
        
        MeinekeCryptCell new_cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
        new_cell.SetBirthTime(0);
        
        c_vector<double,2> new_location;
        new_location[0] = 0.3433453454443;
        new_location[0] = 0.3435346344234;
        crypt.AddCell(new_cell, new_location);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(crypt.rGetCells().size(), 82u);
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 71u);

        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);
        p_simulation_time->IncrementTimeOneStep();

        unsigned num_removed = crypt.RemoveDeadCells();
        TS_ASSERT_EQUALS(num_removed, 1u);


        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(crypt.rGetCells().size(), 81u);
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 70u);

        MeinekeCryptCell new_cell2(STEM, HEALTHY, 0, new FixedCellCycleModel());
        new_cell2.SetBirthTime(0);
        
        c_vector<double,2> new_location2;
        new_location2[0] = 0.6433453454443;
        new_location2[0] = 0.6435346344234;
        crypt.AddCell(new_cell2, new_location2);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(crypt.rGetCells().size(), 82u);
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 71u);
        
        SimulationTime::Destroy();
    }

    // test update ghost node positions
    
    void TestOutputWriters()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);        
        
        // create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(0);
            cells.push_back(cell);
        }
     
        Crypt<2> crypt(mesh,cells);
        
        std::string output_directory = "TestCryptWriters";
        OutputFileHandler output_file_handler(output_directory, false);
        ColumnDataWriter tabulated_node_writer(output_directory, "tab_node_results");
        ColumnDataWriter tabulated_element_writer(output_directory, "tab_elem_results", false);

        TS_ASSERT_THROWS_NOTHING(crypt.SetupTabulatedWriters(tabulated_node_writer, tabulated_element_writer));
                
        out_stream p_node_file = output_file_handler.OpenOutputFile("results.viznodes");
        out_stream p_element_file = output_file_handler.OpenOutputFile("results.vizelements");
      
        crypt.WriteResultsToFiles(tabulated_node_writer,
                                  tabulated_element_writer,
                                  *p_node_file,
                                  *p_element_file,
                                  true,
                                  true);
        p_node_file->close();                          
        p_element_file->close();                          

        // compare output with saved files of what they should look like                           
        std::string results_dir = output_file_handler.GetTestOutputDirectory();

        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "results.vizelements  models/test/data/TestCryptWriters/results.vizelements").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "results.viznodes     models/test/data/TestCryptWriters/results.viznodes").c_str()), 0);

        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "tab_node_results.dat models/test/data/TestCryptWriters/tab_node_results.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "tab_elem_results.dat models/test/data/TestCryptWriters/tab_elem_results.dat").c_str()), 0);
        SimulationTime::Destroy();        
    }
    
    void TestSpringIterator2d() throw(Exception)
    {
        // set up expected results for the honeycombmesh created below
        // the following are the edges which do not contain a ghost node
        std::set < std::set < unsigned > > expected_node_pairs;        
        unsigned expected_node_pairs_array[] = {5,6,
                                                5,9,
                                                5,10,
                                                6,10,
                                                9,10 };

        for (unsigned i=0; i<10; i=i+2)
        {
            std::set < unsigned > node_pair;
            node_pair.insert(expected_node_pairs_array[i]);
            node_pair.insert(expected_node_pairs_array[i+1]);
            expected_node_pairs.insert(node_pair);
        }
        
        // set up simple crypt with honeycomb mesh
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);        
        
        unsigned num_cells_depth = 2;
        unsigned num_cells_width = 2;
        unsigned thickness_of_ghosts = 1;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghosts, false);

        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(0);
            cells.push_back(cell);
        }
        
        // create a crypt, with no ghost nodes at the moment
        Crypt<2> crypt(*p_mesh,cells);
        
        // get ghost node indices from the generator and set up the ghost nodes
        // vector of bools
        std::vector<bool> is_ghost_node(p_mesh->GetNumNodes(),false);
        for(std::set<unsigned>::iterator iter = ghost_node_indices.begin();
            iter!=ghost_node_indices.end();
            iter++)
        {
            is_ghost_node[*iter] = true;
        }
        
        // set ghost nodes
        crypt.SetGhostNodes(is_ghost_node);
                
        // check that we can iterate over the set of springs
        std::set< std::set< unsigned > > springs_visited;
        
        for (Crypt<2>::SpringIterator spring_iterator=crypt.SpringsBegin();
             spring_iterator!=crypt.SpringsEnd();
             ++spring_iterator)
        {
            std::set<unsigned> node_pair;
            node_pair.insert(spring_iterator.GetNodeA()->GetIndex());
            node_pair.insert(spring_iterator.GetNodeB()->GetIndex());
            
            TS_ASSERT_EQUALS(springs_visited.find(node_pair), springs_visited.end());
            springs_visited.insert(node_pair);
            
            TS_ASSERT_EQUALS(spring_iterator.rGetCellA().GetNodeIndex(), spring_iterator.GetNodeA()->GetIndex());
            TS_ASSERT_EQUALS(spring_iterator.rGetCellB().GetNodeIndex(), spring_iterator.GetNodeB()->GetIndex());
        }
        
         TS_ASSERT_EQUALS(springs_visited, expected_node_pairs);

        SimulationTime::Destroy();        
    }

    // 3d test with some ghost nodes
    void TestSpringIterator3d() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);        
        
        // create a simple mesh
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(0);
            cells.push_back(cell);
        }
        
        // create a crypt, with no ghost nodes at the moment
        Crypt<3> crypt(mesh,cells);
        
        // make nodes 0-10 ghost nodes
        std::vector<bool> is_ghost_node(mesh.GetNumNodes(),false);
        for(unsigned i=0; i<10; i++)
        {
            is_ghost_node[i] = true;
        }
        // set ghost nodes
        crypt.SetGhostNodes(is_ghost_node);

                
        // check that we can iterate over the set of springs
        std::set< std::set< unsigned > > springs_visited;
        
        for (Crypt<3>::SpringIterator spring_iterator=crypt.SpringsBegin();
             spring_iterator!=crypt.SpringsEnd();
             ++spring_iterator)
        {
            std::set<unsigned> node_pair;
            node_pair.insert(spring_iterator.GetNodeA()->GetIndex());
            node_pair.insert(spring_iterator.GetNodeB()->GetIndex());
            
            TS_ASSERT_EQUALS(springs_visited.find(node_pair), springs_visited.end());
            springs_visited.insert(node_pair);
            
            TS_ASSERT_EQUALS(spring_iterator.rGetCellA().GetNodeIndex(), spring_iterator.GetNodeA()->GetIndex());
            TS_ASSERT_EQUALS(spring_iterator.rGetCellB().GetNodeIndex(), spring_iterator.GetNodeB()->GetIndex());
        }
        
        // set up expected node pairs
        std::set< std::set<unsigned> > expected_node_pairs;
        for(unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            Element<3,3>* p_element = mesh.GetElement(i);
            for(unsigned j=0; j<4; j++)
            {
                for(unsigned k=0; k<4; k++)
                {
                    unsigned node_A = p_element->GetNodeGlobalIndex(j);
                    unsigned node_B = p_element->GetNodeGlobalIndex(k);
                    
                    // if nodeA or node_B are <10 they will have been labelled a ghost node
                    // above
                    if(node_A != node_B && node_A>=10 && node_B>=10)
                    {
                        std::set<unsigned> node_pair;
                        node_pair.insert(node_A);
                        node_pair.insert(node_B);
                        
                        expected_node_pairs.insert(node_pair);
                    }
                }
            }
        }
        
        TS_ASSERT_EQUALS(springs_visited, expected_node_pairs);

        SimulationTime::Destroy();        
    }
};



#endif /*TESTCRYPT_HPP_*/
