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
#include "CryptHoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"
#include "Crypt.cpp"
#include "CryptHoneycombMeshGenerator.hpp"

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
    
    // test with ghost nodes, incl the Iterator doesn't loop over ghost nodes
    void TestCryptWithGhostNodes() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);        
        
        unsigned num_cells_depth = 11;
        unsigned num_cells_width = 6;
        
        CryptHoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 6.0, 2u, false);

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
        std::cout << "counter = " << counter << "\n";
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
        Crypt<2> crypt(mesh,cells);

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
        unsigned old_num_cells = cells.size();

        // create a new cell, DON'T set the node index, set birth time=-1
        MeinekeCryptCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
        cell.SetBirthTime(-1);
        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 2;
        new_cell_location[1] = 2;
        
        crypt.AddCell(cell,new_cell_location); 

        // crypt should have updated mesh and cells
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), old_num_nodes+1);
        TS_ASSERT_EQUALS(cells.size(), old_num_cells+1);

        // same test via crypt class
        TS_ASSERT_EQUALS(crypt.rGetMesh().GetNumNodes(), old_num_nodes+1);
        TS_ASSERT_EQUALS(crypt.rGetCells().size(), old_num_cells+1);
        
        // check the location of the new node
        TS_ASSERT_DELTA(mesh.GetNode(old_num_nodes)->rGetLocation()[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(mesh.GetNode(old_num_nodes)->rGetLocation()[1], 2.0, 1e-12);

        // check the index of the new cell
        TS_ASSERT_EQUALS(cells[cells.size()-1].GetNodeIndex(), old_num_nodes);
        
        SimulationTime::Destroy();
    }
    
    // test remove dead cells
    void TestRemoveDeadCells()
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
        // create a crypt, with no ghost nodes at the moment
        Crypt<2> crypt(mesh,cells);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(crypt.rGetCells().size(), 81u);

        // check the iterator still works 
        unsigned counter = 0;
        for (Crypt<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            counter++;
        }
        TS_ASSERT_EQUALS(counter, 81u);

        
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);
        p_simulation_time->IncrementTimeOneStep();
        crypt.RemoveDeadCells();
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 80u);
        TS_ASSERT_EQUALS(crypt.rGetCells().size(), 80u);
        TS_ASSERT_EQUALS(crypt.rGetCells().size(), cells.size());
        
     //   crypt.Remesh();
        // check the iterator still works 
        counter = 0;
        for (Crypt<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            counter++;
        }
        TS_ASSERT_EQUALS(counter, 80u);
        
        // test size of ghost nodes vector is correct - mesh.GetNumAllNodes() ?
        
        // crypt.Remesh();
        
        // test size of ghost nodes vector is correct - mesh.GetNumNodes() ?
        
        SimulationTime::Destroy();       
    }
    // test add and remove, remove and add
    
    
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
    }
};



#endif /*TESTCRYPT_HPP_*/
