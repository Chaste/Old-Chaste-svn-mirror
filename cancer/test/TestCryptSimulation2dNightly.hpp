#ifndef TESTCRYPTSIMULATION2DNIGHTLY_HPP_
#define TESTCRYPTSIMULATION2DNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "CryptSimulation2d.hpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "OutputFileHandler.hpp"
#include "TissueCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "WntConcentration.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "MeshBasedTissue.cpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "CellsGenerator.hpp"
#include "AbstractCancerTestSuite.hpp"

// Simple cell killer which just kills a single cell.
// The constructor takes in a number, and the killer
// will kill the number-th cell reached using the iterator
// (or the last cell, if n>num_cells)
//
// For testing purposes
class SingleCellCellKiller : public AbstractCellKiller<2>
{
private :
    unsigned mNumber;

public :
    SingleCellCellKiller(AbstractTissue<2>* pTissue, unsigned number)
        : AbstractCellKiller<2>(pTissue),
          mNumber(number)
    {
    }
    
    virtual void TestAndLabelCellsForApoptosisOrDeath()
    {
        if(mpTissue->GetNumRealCells()==0)
        {
            return;
        }
        
        MeshBasedTissue<2>::Iterator cell_iter = mpTissue->Begin();
       
        for(unsigned i=0; ( (i<mNumber) && (cell_iter!=mpTissue->End()) ); i++)
        {
            ++cell_iter;
        }
        
        cell_iter->Kill();
    }
};


class TestCryptSimulation2dNightly : public AbstractCancerTestSuite
{
private:
  
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

///////// NON-PERIODIC TESTS - These test the spring system and cell birth etc. ----------

    // Test the spring system. The cells in this test are given an intial
    // age of 2.0 so that their springs are at their natural length
    // i.e. we set birth time=-2.0. 
    // The mesh is initially a set of 10 by 10 squares, each square made
    // up of two triangles. The horizontal and vertical edges (springs) are at rest length, the
    // diagonals are two long, so this means the mesh skews to a (sloughed) parallelogram, each
    // triangle trying to become equilateral.
    //
    // If you want to view the results visually set the end time
    // to 24.0 and it will look like a parallelogram.
    // However we keep the simulation time at 1.0 to make
    // the test short.
    void Test2DSpringSystem() throw (Exception)
    {
        double crypt_length = 10;
        double crypt_width = 10;
        CancerParameters::Instance()->SetCryptLength(crypt_length);
        CancerParameters::Instance()->SetCryptWidth(crypt_width);
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, mesh, FIXED, false, 0.0, 3.0, 6.5, 8.0);
        
        MeshBasedTissueWithGhostNodes<2> crypt(mesh, cells);
        CryptSimulation2d simulator(crypt);    
                
        // Destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        simulator.SetEndTime(1.0);
        TS_ASSERT_THROWS_ANYTHING(simulator.Solve());// fails because output directory not set
                
        // Destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        
        simulator.SetOutputDirectory("Crypt2DSprings");

        simulator.SetEndTime(1.0);
        
        simulator.SetReMeshRule(false);
        simulator.SetNoBirth(true);

        // Destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        
        simulator.Solve();
        std::vector<double> node_0_location = simulator.GetNodeLocation(0);
        TS_ASSERT_DELTA(node_0_location[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(node_0_location[1], 0.0, 1e-12);
        
        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt2DSprings",false);
        std::string node_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        TS_ASSERT_EQUALS(system(("diff " + node_results_file + " cancer/test/data/Crypt2DSpringsResults/results.viznodes").c_str()), 0);
         
        std::string elem_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizelements";
        TS_ASSERT_EQUALS(system(("diff " + elem_results_file + " cancer/test/data/Crypt2DSpringsResults/results.vizelements").c_str()), 0);          
    }
    
    void Test2DHoneycombMeshNotPeriodic() throw (Exception)
    {
        int num_cells_depth = 11;
        int num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        CancerParameters::Instance()->SetCryptLength(crypt_length);
        CancerParameters::Instance()->SetCryptWidth(crypt_width);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
                
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(crypt);
                
        simulator.SetOutputDirectory("Crypt2DHoneycombMesh");
        simulator.SetEndTime(12.0);
        
        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);
       
        simulator.Solve();
        
        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt2DHoneycombMesh",false);
        std::string node_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        TS_ASSERT_EQUALS(system(("diff " + node_results_file + " cancer/test/data/Crypt2DHoneycombMeshResults/results.viznodes").c_str()), 0);
         
        std::string elem_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizelements";
        TS_ASSERT_EQUALS(system(("diff " + elem_results_file + " cancer/test/data/Crypt2DHoneycombMeshResults/results.vizelements").c_str()), 0);   
    }
    
    void TestMonolayer() throw (Exception)
    {
        int num_cells_depth = 11;
        int num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        CancerParameters::Instance()->SetCryptLength(crypt_length);
        CancerParameters::Instance()->SetCryptWidth(crypt_width);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true,-1.0);
                
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);
        crypt.SetWriteVoronoiData(true,true);
        
        // Set the first cell to be logged
        crypt.Begin()->SetLogged();

        CryptSimulation2d simulator(crypt);

        simulator.SetOutputDirectory("Monolayer");
        simulator.SetEndTime(1);
        
        simulator.Solve();
        
        //Check writing of voronoi data
        OutputFileHandler handler("Monolayer",false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizvoronoi";
        TS_ASSERT_EQUALS(system(("diff " + results_file + " cancer/test/data/Monolayer/VoronoiAreaAndPerimeter.dat").c_str()), 0);
    }
    
    //////////////////////////////////////////////////////////////////
    // starting with a small mesh with 1 stem cell and the rest
    // differentiated, check the number of cells at the end of the
    // simulation is as expected.
    //////////////////////////////////////////////////////////////////
    void Test2DCorrectCellNumbers() throw (Exception)
    {
        CancerParameters* p_params = CancerParameters::Instance();
        
        // Check the stem cell cycle time is still 24 hrs, otherwise
        // this test might not pass
        TS_ASSERT_DELTA(p_params->GetStemCellG1Duration(), 14, 1e-12);
        TS_ASSERT_DELTA(p_params->GetTransitCellG1Duration(), 2, 1e-12);
        TS_ASSERT_DELTA(p_params->GetSG2MDuration(), 10, 1e-12);
        
        int num_cells_width = 7;
        int num_cells_depth = 5;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        double crypt_width = num_cells_width - 1.0;
        double crypt_length = num_cells_depth - 1.0;
        
        CancerParameters::Instance()->SetCryptLength(crypt_length);
        CancerParameters::Instance()->SetCryptWidth(crypt_width);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = p_mesh->GetNumAllNodes();
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CellType cell_type;
            unsigned generation;
            double birth_time;
            
            if (i==27) // middle of bottom row of cells
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
            
            TissueCell cell(cell_type, HEALTHY, new FixedCellCycleModel());
            cell.GetCellCycleModel()->SetGeneration(generation);
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DSpringsCorrectCellNumbers");
        simulator.SetEndTime(40); // hours
        
        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);
        
        simulator.Solve();
        
        // Now count the number of each type of cell
        unsigned num_stem = 0;
        unsigned num_transit = 0;
        unsigned num_differentiated = 0;
        
        for (MeshBasedTissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            CellType type = cell_iter->GetCellType();
                
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
        
        TS_ASSERT_EQUALS(num_stem, 1u);
        TS_ASSERT_EQUALS(num_transit, 2u);
        
        TS_ASSERT_LESS_THAN(num_differentiated, 25u);
        TS_ASSERT_LESS_THAN(15u, num_differentiated);
    }
    
    
    
////////////////////////////////////////////////////////////////////////////
// PERIODIC TESTS
// 
// These test the system as a whole
// 
////////////////////////////////////////////////////////////////////////////
    
    void Test2DPeriodicNightly() throw (Exception)
    {        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
       
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
               
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DPeriodicNightly");
        simulator.SetEndTime(12.0);
        
        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);
       
        simulator.Solve();
        
        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
        
        TS_ASSERT_EQUALS(number_of_cells, 85u);
        TS_ASSERT_EQUALS(number_of_nodes, 133u);
    }
    
    
    void TestCrypt2DPeriodicWntNightly() throw (Exception)
    {        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, WNT, true);
        
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);
        
        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(crypt);

        CryptSimulation2d simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DPeriodicWntNightly");
        
        simulator.SetEndTime(24.0);
        
        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);
       
        simulator.Solve();
        
        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)

        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 94u);
        TS_ASSERT_EQUALS(number_of_nodes, 142u);

        WntConcentration::Destroy();
        
        CancerEventHandler::Headings();
        CancerEventHandler::Report();
    }
    

    // This test is dontTest-ed out and not run every night as it doesn't really test
    // anything. It does show how to set up a mutant simulation. Mutant viscosities 
    // are tested elsewhere directly. 
    void dontRunTestWithMutantCellsUsingDifferentViscosities() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, WNT, true);
                
        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            CellMutationState mutation_state;

            double x = p_mesh->GetNode(i)->GetPoint().rGetLocation()[0];
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            double dist_from_3_6 = sqrt((x-3)*(x-3)+(y-6)*(y-6));
            
            if(dist_from_3_6<1.1)
            {
                mutation_state = APC_TWO_HIT;
            }
            else
            {
                mutation_state = HEALTHY;
            }
            
            cells[i].SetMutationState(mutation_state);
        }
        
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(crypt);

        CryptSimulation2d simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DPeriodicMutant");        
        simulator.SetEndTime(12.0);
        
        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);
        
        simulator.Solve();
        
        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<bool> ghost_cells = crypt.rGetGhostNodes();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
        
        TS_ASSERT_EQUALS(number_of_nodes,ghost_cells.size());
        
        unsigned number_of_cells = 0;
        unsigned number_of_mutant_cells = 0;
        for (MeshBasedTissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            number_of_cells++;
            if (cell_iter->GetMutationState()==APC_TWO_HIT)
            {
                number_of_mutant_cells++;
            }
        }
 
        WntConcentration::Destroy();
    }
    
        
    void TestRandomDeathWithPeriodicMesh() throw (Exception)
    {
        unsigned cells_across = 7;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
              
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DRandomDeathPeriodic");
        simulator.SetEndTime(4.6);

        RandomCellKiller<2> random_cell_killer(&crypt, 0.01);
        simulator.AddCellKiller(&random_cell_killer);
    
        simulator.Solve();
        
        // There should be no cells left after this amount of time
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 0u);
    }
    
  
    // Sloughing with a sloughing cell killer and not turning into ghost nodes
    // on a non-periodic mesh
    void TestSloughingCellKillerOnNonPeriodicCrypt() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
              
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DSloughingDeathNonPeriodic");
        simulator.SetEndTime(4.0);

        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);
        
        simulator.Solve();
    }


    void TestSloughingDeathWithPeriodicMesh() throw (Exception)
    {
        unsigned cells_across = 7;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer,true,crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
              
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DSloughingDeathPeriodic");
        simulator.SetEndTime(4.0);

        SloughingCellKiller cell_killer(&crypt);
        simulator.AddCellKiller(&cell_killer);
               
        simulator.Solve();
        
        std::vector<bool> ghost_node_indices_after = crypt.rGetGhostNodes();
        unsigned num_ghosts=0;
        for (unsigned i=0; i < ghost_node_indices_after.size() ; i++)
        {
            if (ghost_node_indices_after[i])
            {
                num_ghosts++;
            }
        }
        
        // Check no new ghost nodes have been created.
        TS_ASSERT_EQUALS(num_ghosts, ghost_node_indices.size());
        
        // There should be this number of cells left after this amount of time
        // (we have lost two rows of 7 but had a bit of birth too)
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 85u);
    }


    void TestWithMultipleCellKillers() throw (Exception)
    {
        CancerParameters::Instance()->Reset();

        unsigned cells_across = 7;
        unsigned cells_up = 11;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
              
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(crypt);
        
        simulator.SetOutputDirectory("CryptWithMultipleCellKillers");

        // These killers are defined in this test. They kill the first and second
        // available cell, respectively.
        SingleCellCellKiller cell_killer1(&crypt,0);
        SingleCellCellKiller cell_killer2(&crypt,1);

        simulator.AddCellKiller(&cell_killer1);
        simulator.AddCellKiller(&cell_killer2);

        // Just enough time to kill off all the cells, as two are killed per timestep
        double dt = 0.01;
        unsigned num_cells = crypt.GetNumRealCells();
        simulator.SetDt(dt);
        simulator.SetEndTime(0.5*dt*num_cells);

        simulator.Solve();
        
        std::vector<bool> ghost_node_indices_after = crypt.rGetGhostNodes();
        unsigned num_ghosts=0;        
        for (unsigned i=0; i < ghost_node_indices_after.size() ; i++)
        {
            if (ghost_node_indices_after[i])
            {
                num_ghosts++;
            }
        }
        
        // Check no new ghost nodes have been created.
        TS_ASSERT_EQUALS(num_ghosts, ghost_node_indices.size());

        // All cells should have been removed in this time
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 0u);
    }
        
        
    void TestMonolayerWithCutoffPointAndNoGhosts() throw (Exception)
    {
        int num_cells_depth = 11;
        int num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        CancerParameters::Instance()->SetCryptLength(crypt_length);
        CancerParameters::Instance()->SetCryptWidth(crypt_width);
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true,-1.0);
                
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        // Set up crypt simulation
        Meineke2001SpringSystem<2> spring_system(crypt);
        spring_system.UseCutoffPoint(sqrt(2)); // root2 is a sensible choice
        
        CryptSimulation2d simulator(crypt, &spring_system);

        simulator.SetOutputDirectory("MonolayerCutoffPointNoGhosts");
        simulator.SetEndTime(12.0);
        
        simulator.Solve();
    }  
};


#endif /*TESTCRYPTSIMULATION2DNIGHTLY_HPP_*/
