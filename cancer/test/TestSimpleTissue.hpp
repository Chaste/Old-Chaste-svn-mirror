#ifndef TESTSIMPLETISSUE_HPP_
#define TESTSIMPLETISSUE_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cmath>
#include <vector>
#include "SimpleTissue.cpp"
#include "TrianglesMeshReader.cpp"
#include "OutputFileHandler.hpp"
#include "CellsGenerator.hpp"
#include "FixedCellCycleModel.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestSimpleTissue : public AbstractCancerTestSuite
{    
private: 

    template<unsigned DIM>
    std::vector<Node<DIM> > SetUpNodes(ConformingTetrahedralMesh<DIM,DIM>* pMesh)
    {
        std::vector<Node<DIM> > nodes;
        
        for (unsigned i=0; i<pMesh->GetNumNodes(); i++)
        {
            nodes.push_back(*(pMesh->GetNode(i)));
        }
        return nodes;
    }
    
    template<unsigned DIM>
    std::vector<TissueCell> SetUpCells(ConformingTetrahedralMesh<DIM,DIM>* pMesh)
    {
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<pMesh->GetNumNodes(); i++)
        {            
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            double birth_time = 0.0-i;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        return cells;
    }
    
    template<unsigned DIM>
    void TestSimpleSimpleTissue(std::string meshFilename)
    {        
        // Create a simple mesh
        TrianglesMeshReader<DIM,DIM> mesh_reader(meshFilename);
        ConformingTetrahedralMesh<DIM,DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Get node vector from mesh
        std::vector<Node<DIM> > nodes = SetUpNodes(&mesh);
                
        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);

        // Create the tissue
        SimpleTissue<DIM> simple_tissue(nodes, cells);
        
        TS_ASSERT_EQUALS(simple_tissue.rGetNodes().size(), mesh.GetNumNodes());        
        TS_ASSERT_EQUALS(simple_tissue.rGetCells().size(), cells.size());
        
        unsigned counter = 0;
        for (typename SimpleTissue<DIM>::Iterator cell_iter = simple_tissue.Begin();
             cell_iter != simple_tissue.End();
             ++cell_iter)
        {
            // Test operator* and that cells are in sync
            TS_ASSERT_EQUALS((*cell_iter).GetNodeIndex(), counter);
 
            // Test operator-> and that cells are in sync
            TS_ASSERT_DELTA(cell_iter->GetAge(), (double)counter, 1e-12);
 
            // Test GetNode() on the iterator
            TS_ASSERT_EQUALS(cell_iter.GetNode()->GetIndex(), mesh.GetNode(counter)->GetIndex());
 
            // Test iter.GetNode()->GetIndex() is consistent with cell.GetNodeIndex()
            TS_ASSERT_EQUALS((*cell_iter).GetNodeIndex(), cell_iter.GetNode()->GetIndex());
 
            // Test rGetLocation() on the iterator
            for(unsigned space_index=0; space_index<DIM; space_index++)
            {
                TS_ASSERT_EQUALS(cell_iter.rGetLocation()[space_index], 
                                 mesh.GetNode(counter)->rGetLocation()[space_index]);
            }  
            counter++;
        }
        
        TS_ASSERT_EQUALS(counter, simple_tissue.GetNumRealCells());
    }
    

public:

    // Test construction, accessors and Iterator
    void TestSimpleTissue1d2d3d() throw(Exception)
    {
        TestSimpleSimpleTissue<1>("mesh/test/data/1D_0_to_1_10_elements");
        TestSimpleSimpleTissue<2>("mesh/test/data/square_4_elements");
        TestSimpleSimpleTissue<3>("mesh/test/data/cube_136_elements");
    }
        
    void TestValidate()
    {        
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Get node vector from mesh
        std::vector<Node<2> > nodes = SetUpNodes(&mesh);
                
        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);
        cells[0].SetNodeIndex(1);

        // Fails as no cell or ghost correponding to node 0        
        TS_ASSERT_THROWS_ANYTHING(SimpleTissue<2> simple_tissue(nodes, cells));
    }
    
    
    void TestMoveCellAndAddCell()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Get node vector from mesh
        std::vector<Node<2> > nodes = SetUpNodes(&mesh);
                
        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);

        // Create a tissue
        SimpleTissue<2> simple_tissue(nodes, cells);

        // Test move cell by moving node 0 by a small amount
        
        SimpleTissue<2>::Iterator cell_iter = simple_tissue.Begin();
        c_vector<double,2> new_location = cell_iter.rGetLocation();
        new_location[0] += 1e-2;
        new_location[1] += 1e-2;
        ChastePoint<2> new_location_point(new_location);
        simple_tissue.MoveCell(cell_iter, new_location_point);
        
        TS_ASSERT_DELTA(simple_tissue.GetNode(0)->rGetLocation()[0], new_location[0], 1e-12);
        TS_ASSERT_DELTA(simple_tissue.GetNode(0)->rGetLocation()[1], new_location[1], 1e-12);

        // Test AddNode

        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = 1.71;
        new_point.rGetLocation()[1] = 1.72;
        
        unsigned num_nodes = simple_tissue.GetNumNodes();
        
        Node<2>* p_node = new Node<2>(num_nodes, new_point, false);
        unsigned new_node_index = simple_tissue.AddNode(p_node);
        
        TS_ASSERT_EQUALS(new_node_index, num_nodes);
        TS_ASSERT_DELTA(simple_tissue.GetNode(num_nodes)->rGetLocation()[0], 1.71, 1e-12);
        TS_ASSERT_DELTA(simple_tissue.GetNode(num_nodes)->rGetLocation()[1], 1.72, 1e-12);
        
        simple_tissue.mNodes.pop_back();        
        
//        \todo: testRemoveNode
 
        // Test AddCell

        unsigned old_num_nodes = simple_tissue.rGetNodes().size();
        unsigned old_num_cells = simple_tissue.rGetCells().size();

        // Create a new cell, DON'T set the node index, set birth time=-1
        TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
        cell.SetBirthTime(-1);
        
        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 2.0;
        new_cell_location[1] = 2.0;
        
        simple_tissue.AddCell(cell,new_cell_location); 

        // Tissue should have updated nodes and cells
        TS_ASSERT_EQUALS(simple_tissue.GetNumNodes(), old_num_nodes+1);
        TS_ASSERT_EQUALS(simple_tissue.rGetCells().size(), old_num_cells+1);
        TS_ASSERT_EQUALS(simple_tissue.GetNumRealCells(), old_num_nodes+1);
     
        // Check the location of the new node
        TS_ASSERT_DELTA(simple_tissue.GetNode(old_num_nodes)->rGetLocation()[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(simple_tissue.GetNode(old_num_nodes)->rGetLocation()[1], 2.0, 1e-12);

        // Check the index of the new cell
        TissueCell& new_cell = simple_tissue.rGetCells().back();
        TS_ASSERT_EQUALS(new_cell.GetNodeIndex(), old_num_nodes);
        
        delete p_node;
    }
    
    
    void TestRemoveDeadCellsAndUpdateNodeCellMap()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);
        
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Get node vector from mesh
        std::vector<Node<2> > nodes = SetUpNodes(&mesh);
                
        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);
        
        cells[27].StartApoptosis();

        // Create a tissue
        SimpleTissue<2> simple_tissue(nodes, cells);

        // Test we have the right numbers of nodes and cells
        TS_ASSERT_EQUALS(simple_tissue.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(simple_tissue.GetNumRealCells(), 81u);

        p_simulation_time->IncrementTimeOneStep();
        
        unsigned num_removed = simple_tissue.RemoveDeadCells();
        
        // Test that one cell has been removed
        TS_ASSERT_EQUALS(num_removed, 1u);
        TS_ASSERT_EQUALS(simple_tissue.GetNumRealCells(), 80u);
        
        // Test that one node has been removed
        TS_ASSERT_EQUALS(simple_tissue.GetNumNodes(), 80u);
        
        // Test that we have removed the correct node
        unsigned index = 0;
        for (SimpleTissue<2>::Iterator cell_iter = simple_tissue.Begin();
             cell_iter != simple_tissue.End();
             ++cell_iter)
        {
            if (index < 27)
            {
                TS_ASSERT_EQUALS(cell_iter->GetNodeIndex(), index);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->GetNodeIndex(), index+1);                
            }
            index++;
        }
        
        // Test that each cell's node index has been correctly updated
        // The cell at node 28 should now be at 27 as node 27 was deleted
        simple_tissue.UpdateNodeCellMap();
        
        // Test that we have removed the correct node
        index = 0;
        for (SimpleTissue<2>::Iterator cell_iter = simple_tissue.Begin();
             cell_iter != simple_tissue.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_iter->GetNodeIndex(), index);
            index++;
        }
    }
    
    void TestAddAndRemoveAndAddWithOutUpdatingNodeCellMap()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);
        
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Get node vector from mesh
        std::vector<Node<2> > nodes = SetUpNodes(&mesh);
                
        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);
        
        // Make one cell start apoptosis
        cells[27].StartApoptosis();
        
        TissueCell new_cell(STEM, HEALTHY, new FixedCellCycleModel());
        new_cell.SetBirthTime(0);
        
        // Create a tissue
        SimpleTissue<2> simple_tissue(nodes, cells);
        
        // Test we have the right numbers of nodes and cells
        TS_ASSERT_EQUALS(simple_tissue.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(simple_tissue.GetNumRealCells(), 81u);
        
        // Add a cell to the tissue
        c_vector<double,2> new_location;
        new_location[0] = 0.3433453454443;
        new_location[0] = 0.3435346344234;
        simple_tissue.AddCell(new_cell, new_location);
        
        // Test that the numbers of nodes and cells has been updated
        TS_ASSERT_EQUALS(simple_tissue.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(simple_tissue.GetNumRealCells(), 82u);
        
        p_simulation_time->IncrementTimeOneStep();

        // Test that the apoptotic cell has been removed
        unsigned num_removed = simple_tissue.RemoveDeadCells();
        TS_ASSERT_EQUALS(num_removed, 1u);
        TS_ASSERT_EQUALS(simple_tissue.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(simple_tissue.GetNumRealCells(), 81u);

        // Add another cell to the tissue
        TissueCell new_cell2(STEM, HEALTHY, new FixedCellCycleModel());
        new_cell2.SetBirthTime(0);
        
        c_vector<double,2> new_location2;
        new_location2[0] = 0.6433453454443;
        new_location2[0] = 0.6435346344234;
        simple_tissue.AddCell(new_cell2, new_location2);

        // Test that the numbers of nodes and cells has been updated
        TS_ASSERT_EQUALS(simple_tissue.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(simple_tissue.GetNumRealCells(), 82u);
    }
    
    void TestSettingCellAncestors() throw (Exception)
    {        
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Get node vector from mesh
        std::vector<Node<2> > nodes = SetUpNodes(&mesh);
                
        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);
                
        // Create a tissue
        SimpleTissue<2> simple_tissue(nodes, cells);
        
        // Test that the tissue makes all cells fix the node index as ancestor
        simple_tissue.SetCellAncestorsToNodeIndices();
        
        unsigned counter = 0;
        
        for (SimpleTissue<2>::Iterator cell_iter = simple_tissue.Begin();
             cell_iter != simple_tissue.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_iter->GetAncestor(), cell_iter->GetNodeIndex());
            counter ++;
        }
        TS_ASSERT_EQUALS(counter, 5u);
        
        // Test that we can recover remaining number of ancestors
        std::set<unsigned> remaining_ancestors = simple_tissue.GetCellAncestors();
        TS_ASSERT_EQUALS(remaining_ancestors.size(), 5u);
        
        // Test that the set correctly represents a monoclonal population
        for(SimpleTissue<2>::Iterator cell_iter = simple_tissue.Begin();
             cell_iter != simple_tissue.End();
             ++cell_iter)
        {  
            // Set all cells to have the same ancestor...
            cell_iter->SetAncestor(1u);
        }
        remaining_ancestors = simple_tissue.GetCellAncestors();
        TS_ASSERT_EQUALS(remaining_ancestors.size(), 1u);
    }
    
    void TestGetLocationOfCell() throw (Exception)
    {        
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Get node vector from mesh
        std::vector<Node<2> > nodes = SetUpNodes(&mesh);
                
        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);
                
        // Create a tissue
        SimpleTissue<2> simple_tissue(nodes, cells);
                
        // Loop over nodes
        for (SimpleTissue<2>::Iterator cell_iter = simple_tissue.Begin();
             cell_iter != simple_tissue.End();
             ++cell_iter)
        {
            // Record node location
            c_vector<double , 2> node_location = cell_iter.GetNode()->rGetLocation();
            
            // Get cell at each node
            TissueCell& r_cell = simple_tissue.rGetCellAtNodeIndex(cell_iter.GetNode()->GetIndex());      
            
            // Test GetLocationOfCell()
            TS_ASSERT_DELTA(node_location[0] , simple_tissue.GetLocationOfCell(r_cell)[0] , 1e-9);
            TS_ASSERT_DELTA(node_location[1] , simple_tissue.GetLocationOfCell(r_cell)[1] , 1e-9);
        }
    }
      
    void TestSimpleTissueOutputWriters()
    {        
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Get node vector from mesh
        std::vector<Node<2> > nodes = SetUpNodes(&mesh);
                
        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);

        // Create a tissue
        SimpleTissue<2> simple_tissue(nodes, cells);
        
        // For coverage of WriteResultsToFiles()
        simple_tissue.rGetCellAtNodeIndex(0).SetCellType(TRANSIT);
        simple_tissue.rGetCellAtNodeIndex(0).SetMutationState(LABELLED);
        simple_tissue.rGetCellAtNodeIndex(1).SetCellType(DIFFERENTIATED);
        simple_tissue.rGetCellAtNodeIndex(1).SetMutationState(APC_ONE_HIT);
        simple_tissue.rGetCellAtNodeIndex(2).SetMutationState(APC_TWO_HIT);
        simple_tissue.rGetCellAtNodeIndex(3).SetMutationState(BETA_CATENIN_ONE_HIT);
        simple_tissue.rGetCellAtNodeIndex(4).SetCellType(NECROTIC);
        simple_tissue.rGetCellAtNodeIndex(4).StartApoptosis();
        
        std::string output_directory = "TestSimpleTissueWriters";
        OutputFileHandler output_file_handler(output_directory, false);
                
        TS_ASSERT_THROWS_NOTHING(simple_tissue.CreateOutputFiles(output_directory, false, true, true, false, true));

        simple_tissue.WriteResultsToFiles(true, true, false, true);

        TS_ASSERT_THROWS_NOTHING(simple_tissue.CloseOutputFiles(true, true, false, true));

        // Compare output with saved files of what they should look like                           
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes     cancer/test/data/TestSimpleTissueWriters/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellmutationstates.dat     cancer/test/data/TestSimpleTissueWriters/cellmutationstates.dat").c_str()), 0);
        
        // Test the GetCellMutationStateCount function
        c_vector<unsigned,5> cell_mutation_states = simple_tissue.GetCellMutationStateCount();
        for (unsigned i=1; i<NUM_CELL_MUTATION_STATES; i++)
        {
            TS_ASSERT_EQUALS(cell_mutation_states[i], 1u);
        }
        
        // Test the GetCellTypeCount function
        c_vector<unsigned,5> cell_types = simple_tissue.GetCellTypeCount();
        TS_ASSERT_EQUALS(cell_types[0], 2u);
        TS_ASSERT_EQUALS(cell_types[1], 1u);
        TS_ASSERT_EQUALS(cell_types[2], 1u);
        TS_ASSERT_EQUALS(cell_types[3], 1u);
        
        // For coverage
        simple_tissue.SetCellAncestorsToNodeIndices();
        TS_ASSERT_THROWS_NOTHING(simple_tissue.WriteResultsToFiles(true, false, false, true));
    }
    
    void TestWritingCellCyclePhases()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Get node vector from mesh
        std::vector<Node<2> > nodes = SetUpNodes(&mesh);
        
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {            
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            double birth_time;
            if (i==1)
            {
                birth_time = -0.5;                
            }
            else if (i==2)
            {
                birth_time = -1.5;
            }
            else if (i==3)
            {
                birth_time = -15.5;
            }
            else
            {
                birth_time = -23.5;
            }
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        cells[0].SetCellType(DIFFERENTIATED);
        
        // Create a tissue
        SimpleTissue<2> simple_tissue(nodes, cells);
        
        // Cells have been given birth times of 0, -1, -2, -3, -4.
        // loop over them to run to time 0.0;
        for (MeshBasedTissue<2>::Iterator cell_iter = simple_tissue.Begin();
             cell_iter != simple_tissue.End();
             ++cell_iter)
        {                
            cell_iter->ReadyToDivide();
        }
        
        std::string output_directory = "TestWritingCellCyclePhases";
        OutputFileHandler output_file_handler(output_directory, false); 
               
        simple_tissue.CreateOutputFiles(output_directory, false, false, false, false, true);        
        simple_tissue.WriteResultsToFiles(false, false, false, true);
        simple_tissue.CloseOutputFiles(false, false, false, true);
        
        // Test the GetCellCyclePhaseCount function
        c_vector<unsigned,5> cell_cycle_phases = simple_tissue.GetCellCyclePhaseCount();
        TS_ASSERT_EQUALS(cell_cycle_phases[0], 1u);
        TS_ASSERT_EQUALS(cell_cycle_phases[1], 3u);
        TS_ASSERT_EQUALS(cell_cycle_phases[2], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phases[3], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phases[4], 1u);
   }
    
    void TestArchivingTissue() throw (Exception)
    {    
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "simpletissue.arch";
        
        // Archive a simple tissue 
        {
            // Need to set up time 
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            
            // Create a simple mesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
            ConformingTetrahedralMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            
            // Get node vector from mesh
            std::vector<Node<2> > nodes = SetUpNodes(&mesh);
            
            // Set up cells, one for each node. Get each a birth time of -node_index,
            // so the age = node_index
            std::vector<TissueCell> cells = SetUpCells(&mesh);
                    
            // Create a tissue
            SimpleTissue<2>* const p_tissue = new SimpleTissue<2>(nodes, cells);

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // loop over them to run to time 0.0;
            for (SimpleTissue<2>::Iterator cell_iter = p_tissue->Begin();
                cell_iter != p_tissue->End();
                ++cell_iter)
            {                
                cell_iter->ReadyToDivide();
            }
                        
            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            // Write the tissue to the archive
            output_arch << static_cast<const SimulationTime&> (*p_simulation_time);
            output_arch << p_tissue;
            SimulationTime::Destroy();
            delete p_tissue;
        }
        
        // Restore simple tissue
        {
            // Need to set up time 
            unsigned num_steps = 10;
            
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();
            
            SimpleTissue<2>* p_tissue;
                                                   
            // Restore the tissue
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            input_arch >> *p_simulation_time;
            input_arch >> p_tissue;
            
            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // this checks that individual cells and their models are archived.
            unsigned counter = 0u;
            
            for (SimpleTissue<2>::Iterator cell_iter = p_tissue->Begin();
                 cell_iter != p_tissue->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(), (double)(counter), 1e-7);
                counter++;
            }

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetDimensionalisedTime(), 0.0);
            
            // Check the tissue has been restored
            TS_ASSERT_EQUALS(p_tissue->rGetCells().size(), 5u);

            // Check number of nodes
            std::vector<Node<2> > nodes = p_tissue->rGetNodes();
            TS_ASSERT_EQUALS(nodes.size(), 5u);
            
            // Check some node positions
            TS_ASSERT_EQUALS(nodes[3].GetIndex(), 3u);
            TS_ASSERT_EQUALS(nodes[4].GetIndex(), 4u);
            
            TS_ASSERT_DELTA(nodes[3].rGetLocation()[0], 0.0, 1e-9);
            TS_ASSERT_DELTA(nodes[3].rGetLocation()[1], 1.0, 1e-9);
            TS_ASSERT_DELTA(nodes[4].rGetLocation()[0], 0.5, 1e-9);
            TS_ASSERT_DELTA(nodes[4].rGetLocation()[1], 0.5, 1e-9);
            
            TS_ASSERT_EQUALS(nodes[3].IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(nodes[4].IsBoundaryNode(), false);
            
            delete p_tissue;
        }
    }
    
};

#endif /*TESTSIMPLETISSUE_HPP_*/
