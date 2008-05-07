/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef TESTTISSUESIMULATIONWITHSIMPLETISSUE_HPP_
#define TESTTISSUESIMULATIONWITHSIMPLETISSUE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "TissueSimulation.hpp"
#include "RandomCellKiller.hpp" 
#include "CellsGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestTissueSimulationWithSimpleTissue : public AbstractCancerTestSuite
{
private:

    template<unsigned DIM>
    std::vector<Node<DIM> > SetUpNodes(ConformingTetrahedralMesh<DIM,DIM>* pMesh)
    {
        std::vector<Node<DIM> > nodes;
        
        for(unsigned i=0; i<pMesh->GetNumNodes(); i++)
        {
            nodes.push_back(*(pMesh->GetNode(i)));
        }        
        return nodes;
    }
    
    template<unsigned DIM>
    std::vector<TissueCell> SetUpCells(ConformingTetrahedralMesh<DIM,DIM>* pMesh)
    {   
        std::vector<TissueCell> cells;
        for(unsigned i=0; i<pMesh->GetNumNodes(); i++)
        {
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                (CancerParameters::Instance()->GetStemCellG1Duration()
                                    + CancerParameters::Instance()->GetSG2MDuration() );        
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }        
        return cells;
    }

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCancerTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCancerTestSuite::tearDown();
    }

public:

    /** 
     * Create a simulation of a SimpleTissue with a SimpleTissueMechanicsSystem. 
     * Test that no exceptions are thrown, and write the results to file.
     */
    void TestSimpleMonolayer() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        
        // Get node vector from mesh
        std::vector<Node<2> > nodes = SetUpNodes(p_mesh);
                
        // Set up cells, one for each node. Get each a random birth time.
        std::vector<TissueCell> cells = SetUpCells(p_mesh);

        // Create a simple tissue
        SimpleTissue<2> simple_tissue(nodes, cells);
        
        // For coverage, construct tissue simulation without passing in a mechanics system
        TissueSimulation<2> simulator(simple_tissue);
        simulator.SetOutputDirectory("TestTissueSimulationWithSimpleTissue");
        simulator.SetEndTime(10.0);
        
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        
        // Check that nothing's gone badly wrong by testing that nodes 
        // aren't too close together
        double min_distance_between_cells = 1.0; 
        
        for (unsigned i=0; i<simulator.rGetTissue().GetNumNodes(); i++)
        {
            for (unsigned j=i+1; j<simulator.rGetTissue().GetNumNodes(); j++)
            {
                double distance = norm_2(simulator.rGetTissue().GetNode(i)->rGetLocation()-simulator.rGetTissue().GetNode(j)->rGetLocation());
                if (distance < min_distance_between_cells)
                {
                    min_distance_between_cells = distance;
                }
            }
        }
        
        TS_ASSERT(min_distance_between_cells > CancerParameters::Instance()->GetDivisionSeparation() );
    }

    // results: with a few cells and small end times, Simple was twice as fast as meineke
    //          with 10000 cells, and t_end=0.05, (fixed cell cycle) takes 6.5 mins
    //          => 2 hours real time to do 1hr simulation time
    //   run commented test before to see how meineke does with 10000 cells     

//    void TestSimpleMonolayer2() throw (Exception)
//    {
//        // Create a simple mesh
//        int num_cells_depth = 100;
//        int num_cells_width = 100;        
//        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
//        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
//        
//        // Get node vector from mesh
//        std::vector<Node<2> > nodes = SetUpNodes(p_mesh);
//                
//        // Set up cells, one for each node. Get each a random birth time.
//        std::vector<TissueCell> cells = SetUpCells(p_mesh);
//
//        // Create a simple tissue
//        MeshBasedTissue<2> tissue(*p_mesh, cells);
//        
//        Meineke2001SpringSystem<2> mechanics_system(tissue);
//
//        // Create a tissue simulation
//        TissueSimulation<2> simulator(tissue, &mechanics_system);
//        simulator.SetOutputDirectory("TestSimpleMonolayer2");
//        simulator.SetEndTime(0.05);
//        
//        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
//    }
    
    /** 
     * Create a simulation of a SimpleTissue with a SimpleTissueMechanicsSystem
     * and a CellKiller. Test that no exceptions are thrown, and write the results to file.
     * 
     *  ** NOT RUN as currently fails with positions of some cells becoming NaN ** 
     */
    void failingTestCellDeath() throw (Exception)
    {
        std::string test_folder = "TestTissueSimulationWithSimpleTissueCellDeath";
        LogFile::Instance()->Set(2, test_folder, "log.dat");
        
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        
        // Get node vector from mesh
        std::vector<Node<2> > nodes = SetUpNodes(p_mesh);

        // Set up cells, one for each node. Get each a random birth time.
        std::vector<TissueCell> cells = SetUpCells(p_mesh);

        // Create a simple tissue
        SimpleTissue<2> simple_tissue(nodes, cells);
        
        // Create simple tissue mechanics system (with default cutoff=1.5)
        SimpleTissueMechanicsSystem<2> mechanics_system(simple_tissue);

        // Create a tissue simulation
        TissueSimulation<2> simulator(simple_tissue, &mechanics_system);
        simulator.SetOutputDirectory(test_folder);
        simulator.SetEndTime(0.5);    
      
        // Add cell killer
        RandomCellKiller<2> random_cell_killer(&simple_tissue, 0.05);
        simulator.AddCellKiller(&random_cell_killer);
        
        simulator.Solve();

        // The first change in cell numbers is a death at t=0.258333
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);// this needs changing to suit actual answer when test works!
        LogFile::Close();
    }

    /**
     * Test archiving of a TissueSimulation that uses a SimpleTissue.
     */  
    void TestArchiving() throw (Exception)
    {        
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        
        // Get node vector from mesh
        std::vector<Node<2> > nodes = SetUpNodes(p_mesh);
                
        // Set up cells, one for each node. Get each a random birth time.
        std::vector<TissueCell> cells = SetUpCells(p_mesh);

        // Create a simple tissue
        SimpleTissue<2> simple_tissue(nodes, cells);
        
        // Create simple tissue mechanics system (with default cutoff=1.5)
        SimpleTissueMechanicsSystem<2> mechanics_system(simple_tissue);

        // Create a tissue simulation
        TissueSimulation<2> simulator(simple_tissue, &mechanics_system);
        simulator.SetOutputDirectory("TestTissueSimulationWithSimpleTissueSaveAndLoad");
        simulator.SetEndTime(0.5);
        
        simulator.Solve();
        
        TS_ASSERT_THROWS_NOTHING(simulator.Save());

        TissueSimulation<2>* p_simulator = TissueSimulation<2>::Load("TestTissueSimulationWithSimpleTissueSaveAndLoad", 0.5);
        
        p_simulator->SetEndTime(1.0);
        
        TS_ASSERT_THROWS_NOTHING(p_simulator->Solve());
        
        ///\todo: test results against previous test, once cell death 
        ///       and a stable force law are implemented (see #642 and #678)
        delete p_simulator;
    }
    
};

#endif /*TESTTISSUESIMULATIONWITHSIMPLETISSUE_HPP_*/

