#ifndef TESTCRYPTSIMULATION2D_HPP_
#define TESTCRYPTSIMULATION2D_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "OutputFileHandler.hpp"
#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"
#include "SloughingCellKiller.hpp"
#include "WntGradient.hpp"
#include "VoronoiTessellation.cpp"
#include "CellsGenerator.hpp"



class TestCryptSimulation2d : public CxxTest::TestSuite
{

    /**
     * Compare 2 meshes to see if they are 'the same'.  Doesn't check everything,
     * but is fairly thorough.  Used for testing serialization.
     */
    template<unsigned DIM>
    void CompareMeshes(ConformingTetrahedralMesh<DIM,DIM>* pMesh1,
                       ConformingTetrahedralMesh<DIM,DIM>* pMesh2)
    {
        TS_ASSERT_EQUALS(pMesh1->GetNumAllNodes(), pMesh2->GetNumAllNodes());
        TS_ASSERT_EQUALS(pMesh1->GetNumNodes(), pMesh2->GetNumNodes());
        TS_ASSERT_EQUALS(pMesh1->GetNumBoundaryNodes(), pMesh2->GetNumBoundaryNodes());
        TS_ASSERT_EQUALS(pMesh1->GetNumCornerNodes(), pMesh2->GetNumCornerNodes());
        
        for (unsigned i=0; i<pMesh1->GetNumAllNodes(); i++)
        {
            Node<DIM> *p_node = pMesh1->GetNode(i);
            Node<DIM> *p_node2 = pMesh2->GetNode(i);
            TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
            TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());
            TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());
            for (unsigned j=0; j<DIM; j++)
            {
                TS_ASSERT_DELTA(p_node->rGetLocation()[j], p_node2->rGetLocation()[j], 1e-16);
            }
        }
        
        TS_ASSERT_EQUALS(pMesh1->GetNumElements(), pMesh2->GetNumElements());
        TS_ASSERT_EQUALS(pMesh1->GetNumAllElements(), pMesh2->GetNumAllElements());
        TS_ASSERT_EQUALS(pMesh1->GetNumBoundaryElements(), pMesh2->GetNumBoundaryElements());
        TS_ASSERT_EQUALS(pMesh1->GetNumAllBoundaryElements(), pMesh2->GetNumAllBoundaryElements());
        typename ConformingTetrahedralMesh<DIM,DIM>::ElementIterator it=pMesh1->GetElementIteratorBegin();
        typename ConformingTetrahedralMesh<DIM,DIM>::ElementIterator it2=pMesh2->GetElementIteratorBegin();
        for (;
             it != pMesh1->GetElementIteratorEnd();
             ++it, ++it2)
        {
            Element<DIM,DIM>* p_elt = *it;
            Element<DIM,DIM>* p_elt2 = *it2;
            TS_ASSERT_EQUALS(p_elt->GetNumNodes(), p_elt2->GetNumNodes());
            for (unsigned i=0; i<p_elt->GetNumNodes(); i++)
            {
                TS_ASSERT_EQUALS(p_elt->GetNodeGlobalIndex(i), p_elt2->GetNodeGlobalIndex(i));
            }
        }
    }

public:
    void Test2DCylindrical() throw (Exception)
    {        
        CancerParameters::Instance()->Reset();

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);// true = mature cells

        Crypt<2> crypt(*p_mesh, cells);               
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetEndTime(0.1);
        TS_ASSERT_THROWS_ANYTHING(simulator.SetMaxCells(10));
        simulator.SetMaxCells(500);
        TS_ASSERT_THROWS_ANYTHING(simulator.SetMaxElements(10));
        simulator.SetMaxElements(1000);
        
        simulator.SetOutputCellTypes(true);
                
        // These are for coverage and use the defaults
        simulator.SetDt(1.0/120.0);
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(false);
        simulator.SetOutputDirectory("Crypt2DCylindrical");
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);

        simulator.Solve();

        // test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<bool> ghost_cells = crypt.rGetGhostNodes();
        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
        
        TS_ASSERT_EQUALS(number_of_nodes, ghost_cells.size());
        TS_ASSERT_EQUALS(number_of_cells, cells_across*cells_up+1u);  // 6 cells in a row*12 rows + 1 birth
        TS_ASSERT_EQUALS(number_of_nodes, number_of_cells+thickness_of_ghost_layer*2*cells_across); 

        // Coverage of exceptions (after main test to avoid problems with SimulationTime).
        simulator.SetEndTime(10.0);
        simulator.SetOutputDirectory("");
        TS_ASSERT_THROWS_ANYTHING(simulator.Solve());

        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    
    void TestEdgeLengthBasedSpring() throw (Exception)
    {        
        CancerParameters::Instance()->Reset();

        // Set up the simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0,10u);
                
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 3;
       
        HoneycombMeshGenerator regular_generator(cells_across, cells_up,thickness_of_ghost_layer, false,crypt_width/cells_across);
        ConformingTetrahedralMesh<2,2>* p_regular_mesh=regular_generator.GetMesh();
        std::set<unsigned> regular_ghost_node_indices = regular_generator.GetGhostNodeIndices();      
        
        // Set up cells
        std::vector<MeinekeCryptCell> regular_cells;
        CellsGenerator<2>::GenerateForCrypt(regular_cells, *p_regular_mesh, FIXED, true);// true = mature cells

        Crypt<2> regular_crypt(*p_regular_mesh, regular_cells);               
        regular_crypt.SetGhostNodes(regular_ghost_node_indices);

        TissueSimulation<2> regular_simulator(regular_crypt);
        
        regular_simulator.SetEndTime(1.0);
        TS_ASSERT_THROWS_ANYTHING(regular_simulator.SetMaxCells(10));
        regular_simulator.SetMaxCells(500);
        TS_ASSERT_THROWS_ANYTHING(regular_simulator.SetMaxElements(10));
        regular_simulator.SetMaxElements(1000);
        
        // These are for coverage and use the defaults
        regular_simulator.SetDt(1.0/120.0);
        regular_simulator.SetReMeshRule(true);
        regular_simulator.SetNoBirth(false);
        regular_simulator.SetOutputDirectory("Crypt2DEdgeBasedSpring");
        regular_simulator.SetOutputCellTypes(true);
        
        // check that the force between nodes is correctly calculated when the spring constant is constant (!)
        regular_simulator.SetEdgeBasedSpringConstant(false);
                      
        for(Crypt<2>::SpringIterator spring_iterator=regular_crypt.SpringsBegin();
        spring_iterator!=regular_crypt.SpringsEnd();
        ++spring_iterator)
        {
            
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = regular_simulator.CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1],6.25,1e-3);
        }
        
        // check that the force between nodes is correctly calculated when the spring constant 
        // is proportional to the length of the edge between adjacent cells  
        regular_simulator.SetEdgeBasedSpringConstant(true); 
        regular_simulator.mrCrypt.CreateVoronoiTessellation();  
        
        
        for(Crypt<2>::SpringIterator spring_iterator=regular_crypt.SpringsBegin();
        spring_iterator!=regular_crypt.SpringsEnd();
        ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = regular_simulator.CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1],4.34027778,1e-3);
        }
        
        // choose two interior neighbour nodes
        c_vector<double, 2> force = regular_simulator.CalculateForceBetweenNodes(20u,21u);
        
        TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1],4.34027778,1e-3);
        
        // now move node 21 a bit and check that the force calculation changes correctly
        c_vector<double,2> shift;
        shift[0] = 0.01;
        shift[1] = 0.0;                 
        ChastePoint<2> new_point(regular_simulator.rGetCrypt().rGetMesh().GetNode(21u)->rGetLocation() + shift);
        regular_simulator.rGetCrypt().rGetMesh().SetNode(21u, new_point, false);
        
        // check that the new force between nodes is correctly calculated
        regular_simulator.mrCrypt.CreateVoronoiTessellation();  
        c_vector<double, 2> new_force = regular_simulator.CalculateForceBetweenNodes(20u,21u);
        
        // force calculation: shift is along x-axis so we should have
        // new_edge_length = (5/6 + shift[0])*tan(0.5*arctan(5*sqrt(3)/(5 + 12*shift[0]))),
        // force^2 = mu^2 * (new_edge_length*sqrt(3))^2 * (1 - 5/6 - shift[0])^2
        TS_ASSERT_DELTA(new_force[0]*new_force[0] + new_force[1]*new_force[1], 3.83479824,1e-3);
    
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();    
    }
    
    void TestEdgeBasedSpringsOnPeriodicMesh() throw (Exception)
    {     
        // Test on a periodic mesh
        CancerParameters::Instance()->Reset();

        // Set up the simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0,10u);
                
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 3;
       
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);// true = mature cells

        Crypt<2> crypt(*p_mesh, cells);               
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetEndTime(1.0);
        TS_ASSERT_THROWS_ANYTHING(simulator.SetMaxCells(10));
        simulator.SetMaxCells(500);
        TS_ASSERT_THROWS_ANYTHING(simulator.SetMaxElements(10));
        simulator.SetMaxElements(1000);
        
        // These are for coverage and use the defaults
        simulator.SetDt(1.0/120.0);
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(false);
        simulator.SetOutputDirectory("Crypt2DCylindricalWithEdgeSprings");
        simulator.SetOutputCellTypes(true);
        

        // check that the force between nodes is correctly calculated when the spring constant is constant (!)
        simulator.SetEdgeBasedSpringConstant(false);
                      
        for(Crypt<2>::SpringIterator spring_iterator=crypt.SpringsBegin();
        spring_iterator!=crypt.SpringsEnd();
        ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = simulator.CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
                        
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1],6.25,1e-3);
        }
        
        // check that the force between nodes is correctly calculated when the spring constant 
        // is proportional to the length of the edge between adjacenet cells  
        simulator.SetEdgeBasedSpringConstant(true); 
        simulator.mrCrypt.CreateVoronoiTessellation();  
        for(Crypt<2>::SpringIterator spring_iterator=crypt.SpringsBegin();
        spring_iterator!=crypt.SpringsEnd();
        ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = simulator.CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1],4.34027778,1e-3);
        }              
            
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();                
    }
    
    
    void Test2DCylindricalMultipleDivisions() throw (Exception)
    {   
        // create a log of this test.
        LogFile* p_log_file = LogFile::Instance();
        p_log_file->Set(2,"Crypt2DCylindricalMultipleDivisions");
             
        unsigned cells_across = 6;
        unsigned cells_up = 8;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);// true = mature cells
        
        for(unsigned i=0; i< cells.size();i++)
        {
            cells[i].SetBirthTime(-11.5);
        }
        
        
        Crypt<2> crypt(*p_mesh, cells);               
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true); 
        simulator.AddCellKiller(p_sloughing_cell_killer); 
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        // These are for coverage and use the defaults
        simulator.SetDt(1.0/120.0);
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(false);
        simulator.SetOutputDirectory("Crypt2DCylindricalMultipleDivisions");
        
        simulator.SetEndTime(0.6);
        simulator.Solve();
        //Find the height of the current crypt
        double height_after_division=p_mesh->GetWidth(1);
        simulator.SetEndTime(0.8);
        simulator.Solve();
        //Find the height of the current crypt
        double height_after_relaxation=p_mesh->GetWidth(1);
         
        TS_ASSERT_LESS_THAN(height_after_division, height_after_relaxation);
        

        simulator.SetEndTime(2.0);
        simulator.Solve();
        
        //All fully diffs has sloughed off
        for (Crypt<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
             TS_ASSERT(cell_iter->GetCellType() != DIFFERENTIATED);
        }
        
        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        
        // close the log file opened in this test
        LogFile::Close();
    }

    // This is a rubbish test - all cells start at birthTime = 0.
    // So bizarrely the crypt shrinks as the rest lengths are shortened! 
    // But at least it uses Wnt cell cycle and runs reasonably quickly...
    // For a better test with more randomly distributed cell ages see the Nightly test pack.
    void TestWithWntDependentCells() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 0;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells 
        std::vector<MeinekeCryptCell> cells;                      
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, WNT, false);
        
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);  
        
        WntGradient::Instance()->SetType(LINEAR);             
        WntGradient::Instance()->SetCrypt(crypt);
        
        TissueSimulation<2> simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DPeriodicWnt");
        
        // Set length of simulation here
        simulator.SetEndTime(0.3);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);   
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);
        
        simulator.Solve();
        
        std::vector<double> node_35_location = simulator.GetNodeLocation(35);
        
        TS_ASSERT_DELTA(node_35_location[0], 5.5000 , 1e-4);
        // Old version of this test had cells with age zero, therefore small spring lengths.
        // Variable spring lengths now only associated with cell division.
        TS_ASSERT_DELTA(node_35_location[1], 4.33013 , 1e-4);
          
        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
    }
    
    // A better check that the loaded mesh is the same as that saved
    void TestMeshSurvivesSaveLoad() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        
        // Set up a simulation
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        std::vector<MeinekeCryptCell> cells;
                
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, WNT, false);
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
        
        WntGradient::Instance()->SetType(LINEAR);
        WntGradient::Instance()->SetCrypt(crypt);
        //wnt_gradient.SetCrypt(crypt);
        
        TissueSimulation<2> simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DMeshArchive");
        simulator.SetEndTime(0.1);
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);        
        
        // Memory leak (unconditional jump) without the following line.
        // The archiver assumes that a Solve has been called and simulation time has been set up properly.
        // In this test it hasn't so we need this to avoid memory leak.
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(0.1, 100);

        // Save
        simulator.Save();
        
        // Load
        TissueSimulation<2>* p_simulator;
        p_simulator = TissueSimulation<2>::Load("Crypt2DMeshArchive", 0.0);

        // Create an identical mesh for comparison purposes
        HoneycombMeshGenerator generator2(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh2=generator2.GetCylindricalMesh();
        
        // Compare
        CompareMeshes(p_mesh2, &(p_simulator->rGetCrypt().rGetMesh()));
        
        delete p_simulator;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
    }

    void TestStandardResultForArchivingTestsBelow() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000); 
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
        
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        // We have a Wnt Gradient - but not Wnt dependent cells
        // so that the test runs quickly, but we test archiving of it!
        WntGradient::Instance()->SetType(LINEAR);
        WntGradient::Instance()->SetCrypt(crypt);
        
        TissueSimulation<2> simulator(crypt);

        simulator.SetOutputDirectory("Crypt2DPeriodicStandardResult");
        
        // Our full end time is 0.25
        simulator.SetEndTime(0.25);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);

        simulator.Solve();
        
        // These cells just divided and have been gradually moving apart.
        // These results are from time 0.25, which is also tested below
        // after a save and a load. (To check archiving of mDivisionPairs)
        std::vector<double> node_28_location = simulator.GetNodeLocation(28);
        TS_ASSERT_DELTA(node_28_location[0], 4.2123 , 1e-4);
        TS_ASSERT_DELTA(node_28_location[1], 0.0 , 1e-4);
        std::vector<double> node_120_location = simulator.GetNodeLocation(120);
        TS_ASSERT_DELTA(node_120_location[0], 3.7968 , 1e-4);
        TS_ASSERT_DELTA(node_120_location[1], 0.1050 , 1e-4);
                
        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        
        // test the Wnt gradient result
        MeinekeCryptCell* p_cell = &(crypt.rGetCellAtNodeIndex(28));
        TS_ASSERT_DELTA(WntGradient::Instance()->GetWntLevel(p_cell), 1.0, 1e-9);
        p_cell = &(crypt.rGetCellAtNodeIndex(120));
        TS_ASSERT_DELTA(WntGradient::Instance()->GetWntLevel(p_cell), 0.9898, 1e-4);
        WntGradient::Destroy();
        
    }

    // Testing Save 
    void TestSave() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
        
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
        
        WntGradient::Instance()->SetType(LINEAR);
        WntGradient::Instance()->SetCrypt(crypt);

        TissueSimulation<2> simulator(crypt);

        simulator.SetOutputDirectory("Crypt2DPeriodicSaveAndLoad");
        
        // Our full end time is 0.25, here we run until 0.1 then load and run more below.
        simulator.SetEndTime(0.1);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);

        simulator.Solve();
        
        // save the results..
        simulator.Save();
        
        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
    }
    

    // Testing Load (based on previous two tests)
    void TestLoad() throw (Exception) 
    {
        CancerParameters::Instance()->Reset();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);

        // Load the simulation from the TestSave method above and
        // run it from 0.1 to 0.2
        TissueSimulation<2>* p_simulator1;
        
        WntGradient::Instance();   // Make sure there is no existing Wnt Gradient before load.
        WntGradient::Destroy();
        
        p_simulator1 = TissueSimulation<2>::Load("Crypt2DPeriodicSaveAndLoad", 0.1);
        
        p_simulator1->SetEndTime(0.2);
        p_simulator1->Solve();
        
        // save that then reload
        // and run from 0.2 to 0.25.
        NodeMap map(0) ;
        p_simulator1->rGetCrypt().rGetMesh().ReMesh(map);
        p_simulator1->Save();
        
        TissueSimulation<2>* p_simulator2 = TissueSimulation<2>::Load("Crypt2DPeriodicSaveAndLoad", 0.2);
        
        CompareMeshes(&(p_simulator1->rGetCrypt().rGetMesh()),
                      &(p_simulator2->rGetCrypt().rGetMesh()));
        
        p_simulator2->SetEndTime(0.25);
        p_simulator2->Solve();
        
        
        // These cells just divided and have been gradually moving apart.
        // These results are from time 0.25 in the StandardResult test above.
        std::vector<double> node_28_location = p_simulator2->GetNodeLocation(28);
        TS_ASSERT_DELTA(node_28_location[0], 4.2123 , 1e-4);
        TS_ASSERT_DELTA(node_28_location[1], 0.0 , 1e-4);
        std::vector<double> node_120_location = p_simulator2->GetNodeLocation(120);
        TS_ASSERT_DELTA(node_120_location[0], 3.7968 , 1e-4);
        TS_ASSERT_DELTA(node_120_location[1], 0.1050 , 1e-4);
        
        // test Wnt Gradient was set up correctly
        TS_ASSERT_EQUALS(WntGradient::Instance()->IsGradientSetUp(),true);
        // test the Wnt gradient result
        MeinekeCryptCell* p_cell = &(p_simulator2->rGetCrypt().rGetCellAtNodeIndex(28));
        TS_ASSERT_DELTA(WntGradient::Instance()->GetWntLevel(p_cell), 1.0, 1e-9);
        p_cell = &(p_simulator2->rGetCrypt().rGetCellAtNodeIndex(120));
        TS_ASSERT_DELTA(WntGradient::Instance()->GetWntLevel(p_cell), 0.9898, 1e-4);
        
        delete p_simulator1;
        delete p_simulator2;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        
        WntGradient::Destroy();
    }
    
    /* 
     * Cells are compressed and are trying to spread out. So the cells 
     * at the bottom try to push across y=0 but are prevented from doing so.
     * This is not covered in previous tests because cells are all at age = 0
     * and try to come together to simulate birth.
     * 
     * It is potentially an expensive test computationally, because all 
     * Wnt cells have to run cell cycle models for a large time 
     * to be 'mature' cells which won't shrink together. 
     * Limited this by using only four cells of minimum age.
     */
    void TestWntCellsCannotMoveAcrossYEqualsZeroAndVoronoiWriter() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        
        unsigned cells_across = 2;
        unsigned cells_up = 2;
        double crypt_width = 0.5;
        unsigned thickness_of_ghost_layer = 1;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, false, crypt_width/cells_across);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
                
        // Set up cells
        std::vector<MeinekeCryptCell> cells;        
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, WNT, true);
        
//        cells[0].SetBirthTime(-1.0);   // Make cell cycle models do minimum work
//        cells[1].SetBirthTime(-1.0);
//        cells[1].SetMutationState(LABELLED);
//        cells[2].SetBirthTime(-1.0);
//        cells[2].SetMutationState(APC_ONE_HIT);
//        cells[3].SetBirthTime(-1.0);
//        cells[3].SetMutationState(BETA_CATENIN_ONE_HIT);

        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
        
        Crypt<2>::Iterator cell_iterator = crypt.Begin();
        cell_iterator->SetBirthTime(-1.0);   // Make cell cycle models do minimum work
        ++cell_iterator;
        cell_iterator->SetBirthTime(-1.0);
        cell_iterator->SetMutationState(LABELLED);
        ++cell_iterator;
        cell_iterator->SetBirthTime(-1.0);
        cell_iterator->SetMutationState(APC_ONE_HIT);
        ++cell_iterator;
        cell_iterator->SetBirthTime(-1.0);
        cell_iterator->SetMutationState(BETA_CATENIN_ONE_HIT);
                
        WntGradient::Instance()->SetType(LINEAR);
        WntGradient::Instance()->SetCrypt(crypt);
                
        TissueSimulation<2> simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DWntMatureCells");
        // If you want to visualize this use the 'notcylindrical' option
        // (it is too small for it to figure out what's happening on its own).
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        simulator.SetEndTime(0.01);
        simulator.SetOutputCellTypes(true);   
        // cover the write voronoi data method
        simulator.SetWriteVoronoiData(true, false);     
        simulator.Solve();
        
        // Check that nothing has moved below y=0
        for (Crypt<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            TS_ASSERT_LESS_THAN(-1e-15,p_mesh->GetNode(cell_iter->GetNodeIndex())->rGetLocation()[1]);
        }
        
        c_vector<unsigned,5> cellTypeCount = simulator.GetCellTypeCount();
        TS_ASSERT_EQUALS(cellTypeCount[0],3u); // see ticket:500
        TS_ASSERT_EQUALS(cellTypeCount[1],1u);
        TS_ASSERT_EQUALS(cellTypeCount[2],1u);
        TS_ASSERT_EQUALS(cellTypeCount[3],0u);  // No APC two hit, one of all the rest.
        TS_ASSERT_EQUALS(cellTypeCount[4],1u);
        
        // check writing of voronoi data
        OutputFileHandler handler("Crypt2DWntMatureCells",false);
        std::string results_file = handler.GetTestOutputDirectory() + "VoronoiAreaAndPerimeter.dat";
        TS_ASSERT_EQUALS(system(("cmp " + results_file + " cancer/test/data/Crypt2DWntMatureCells/VoronoiAreaAndPerimeter.dat").c_str()), 0);
            
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
    }
    
    // This is a strange test -- all cells divide within a quick time, it gives
    // good testing of the periodic boundaries though... [comment no longer valid?]
    void TestWithTysonNovakCells() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();      
        // There is no limit on transit cells in T&N
        p_params->SetMaxTransitGenerations(1000);
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, TYSONNOVAK, true);
        
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
        
        TissueSimulation<2> simulator(crypt); 
               
        AbstractCellKiller<2>* p_cell_killer = new SloughingCellKiller(&crypt);
        simulator.AddCellKiller(p_cell_killer);
                
        simulator.SetOutputDirectory("Crypt2DPeriodicTysonNovak");
        simulator.SetEndTime(0.05);
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        simulator.SetDt(0.001);
        
        // Test that labelling a few cells doesn't make any difference to the simulation
        // and therefore log them in the visualizer files for the next test to check.
        simulator.rGetCrypt().rGetCellAtNodeIndex(57).SetMutationState(LABELLED);
        simulator.rGetCrypt().rGetCellAtNodeIndex(56).SetMutationState(APC_ONE_HIT);
        simulator.rGetCrypt().rGetCellAtNodeIndex(51).SetMutationState(APC_TWO_HIT);
        simulator.rGetCrypt().rGetCellAtNodeIndex(63).SetMutationState(BETA_CATENIN_ONE_HIT);
        simulator.SetOutputCellTypes(true);             
        simulator.Solve();
        
        // test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<bool> ghost_cells = crypt.rGetGhostNodes();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
        
        TS_ASSERT_EQUALS(number_of_nodes,ghost_cells.size());
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 81u);
        TS_ASSERT_EQUALS(number_of_nodes, 123u);

        delete p_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    /*
     * This test compares the visualizer output from the previous test with a known file.
     * 
     * Note - if the previous test is changed we need to update the file this test refers to. 
     */
    void TestVisualizerOutput() throw (Exception)
    {
        // work out where the previous test wrote its files
        OutputFileHandler handler("Crypt2DPeriodicTysonNovak",false);
        std::string results_dir = handler.GetTestOutputDirectory() + "results_from_time_0/vis_results";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/results.vizelements cancer/test/data/Crypt2DPeriodicTysonNovak_vis/results.vizelements").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/results.viznodes cancer/test/data/Crypt2DPeriodicTysonNovak_vis/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/results.vizsetup cancer/test/data/Crypt2DPeriodicTysonNovak_vis/results.vizsetup").c_str()), 0);
    }
   
    
    void TestPrivateFunctionsOf2DCryptSimulation() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
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
        CellsGenerator<2>::GenerateForCrypt(cells, mesh, FIXED, false, 0.0, 3.0, 6.5, 8.0);
        
        cells[60].SetBirthTime(-50.0);
        
        Crypt<2> crypt(mesh, cells);
        
        TissueSimulation<2> simulator(crypt);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);
        
        unsigned num_deaths = simulator.DoCellRemoval();
        unsigned num_births = simulator.DoCellBirth();
                                                                
        TS_ASSERT_EQUALS(num_births, 1u);
        TS_ASSERT_EQUALS(num_deaths,11u);
               
        delete p_sloughing_cell_killer;       
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    
    void TestPrivateFunctionsOf2DCryptSimulationOnHoneycombMesh() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        /*
         ************************************************************************
         ************************************************************************ 
         *     Set up a simulation class to run the individual tests on.
         ************************************************************************
         ************************************************************************ 
         */
        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 3;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer,false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        unsigned num_cells = p_mesh->GetNumAllNodes();
        
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator::Instance();
                
        // Set up cells by iterating through the mesh nodes
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            double birth_time;
            CryptCellType cell_type;
            unsigned generation;
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
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
                  
            WntCellCycleModel* p_model = new WntCellCycleModel();
            
            MeinekeCryptCell cell(cell_type, mutation_state, generation, p_model);
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        WntGradient::Instance()->SetType(LINEAR);
        WntGradient::Instance()->SetCrypt(crypt);
        
        TissueSimulation<2> simulator(crypt);
        
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(400);

       /*
        ************************************************************************
        ************************************************************************ 
        *  Test Calculate Velocities on each node
        ************************************************************************
        ************************************************************************ 
        */
                
        std::vector<c_vector<double, 2> > velocities_on_each_node(p_mesh->GetNumAllNodes());
        
        velocities_on_each_node = simulator.CalculateVelocitiesOfEachNode();
 
        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            std::set<unsigned>::iterator iter = ghost_node_indices.find(i);
            bool is_a_ghost_node = (iter!=ghost_node_indices.end());

            if (!is_a_ghost_node)
            {
                TS_ASSERT_DELTA(velocities_on_each_node[i][0], 0.0, 1e-4);
                TS_ASSERT_DELTA(velocities_on_each_node[i][1], 0.0, 1e-4);
            }
        }
        
        // Move a node along the x-axis and calculate the force exerted on a neighbour
        c_vector<double,2> old_point = p_mesh->GetNode(59)->rGetLocation();
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = old_point[0]+0.5;
        new_point.rGetLocation()[1] = old_point[1];
  
        p_mesh->SetNode(59, new_point, false);
        velocities_on_each_node = simulator.CalculateVelocitiesOfEachNode();
        
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
        ChastePoint<2> new_point2;
        new_point2.rGetLocation()[0] = new_point[0] + 0.01;
        new_point2.rGetLocation()[1] = new_point[1] + 0.01 ;
        
        unsigned elem_index = p_mesh->GetContainingElementIndex(new_point2,false);
        Element<2,2>* p_element = p_mesh->GetElement(elem_index);
        
        force_on_spring = simulator.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),p_element->GetNodeGlobalIndex(0));
        TS_ASSERT_DELTA(force_on_spring[0], 0.5*p_params->GetSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);
        
        

       /*
        ************************************************************************
        ************************************************************************ 
        *  Test UpdateNodePositions
        ************************************************************************
        ************************************************************************ 
        */
        
        ChastePoint<2> point_of_node60 = p_mesh->GetNode(60)->rGetLocation();
        
        simulator.SetDt(0.01);
        simulator.UpdateNodePositions(velocities_on_each_node);
        
        TS_ASSERT_DELTA(p_mesh->GetNode(60)->rGetLocation()[0],point_of_node60.rGetLocation()[0]+force_on_spring[0]/p_params->GetDampingConstantMutant() *0.01, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(60)->rGetLocation()[1],point_of_node60.rGetLocation()[1], 1e-4);
        
       /*
        ************************************************************************
        ************************************************************************ 
        * Test UpdateCellTypes is being done
        ************************************************************************
        ************************************************************************ 
        */    
        
        for (Crypt<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            CryptCellType cell_type;
            cell_type = cell_iter->GetCellType();
            if (!cell_type==STEM)
            {
                //std::cout << "Cell type = " << cell_type << std::endl;
                WntCellCycleModel *p_this_model = static_cast<WntCellCycleModel*>(cell_iter->GetCellCycleModel());
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
        
        /*
         ******************************************
         ******************************************
         *  test force with cutoff point
         ******************************************
         ******************************************
         */
        double dist = norm_2( p_mesh->GetVectorFromAtoB(p_element->GetNode(0)->rGetLocation(), p_element->GetNode(1)->rGetLocation()) );   

        simulator.UseCutoffPoint(dist-0.1);
        
        force_on_spring = simulator.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),p_element->GetNodeGlobalIndex(0));
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);
        //\todo The above test is not guaranteed repeatible.  It appears to 
        //give different answers on different runs
        
        //Here's a double check that the geometry is the same:
        TS_ASSERT_DELTA(dist, 0.7607, 1e-4);
        TS_ASSERT_DELTA(p_element->GetNode(0)->rGetLocation()[0], 4.2767, 1e-4);
        TS_ASSERT_DELTA(p_element->GetNode(0)->rGetLocation()[1], 0.8660, 1e-4);
        TS_ASSERT_DELTA(p_element->GetNode(1)->rGetLocation()[0], 5.0375, 1e-4);
        TS_ASSERT_DELTA(p_element->GetNode(1)->rGetLocation()[1], 0.8660, 1e-4);
                

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
    }
    
    void TestCalculateDividingCellCentreLocationsConfMesh() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        CancerParameters::Instance()->SetDivisionRestingSpringLength(0.9);//Only coverage
        CancerParameters::Instance()->SetDivisionSeparation(0.1);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Make a parent node
        c_vector<double ,2> location;
        location[0]=1.0;
        location[1]=1.0;
        Node<2>* p_node = new Node<2>(0u,location, false);

        ConformingTetrahedralMesh<2,2> conf_mesh;
        conf_mesh.AddNode(p_node);
        
        // Set up cells
        std::vector<MeinekeCryptCell> conf_cells;
        CellsGenerator<2>::GenerateForCrypt(conf_cells, conf_mesh, TYSONNOVAK, true);

        Crypt<2> conf_crypt(conf_mesh, conf_cells);

        Crypt<2>::Iterator conf_iter = conf_crypt.Begin();

        TissueSimulation<2> simulator(conf_crypt);
                     
        c_vector<double, 2> daughter_location = simulator.CalculateDividingCellCentreLocations(conf_iter);
        c_vector<double, 2> new_parent_location = conf_mesh.GetNode(0)->rGetLocation();
        c_vector<double, 2> parent_to_daughter = conf_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);
        TS_ASSERT_DELTA(norm_2(parent_to_daughter), 
            CancerParameters::Instance()->GetDivisionSeparation(),
            1e-7);

        SimulationTime::Destroy();
    }
    

    void TestCalculateDividingCellCentreLocationsConfMeshStemCell() throw (Exception)
    {
        CancerParameters::Instance()->Reset();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Make a parent node
        c_vector<double ,2> location;
        location[0]=1.0;
        location[1]=0.0;                                    // <- y=0
        Node<2>* p_node = new Node<2>(0u,location, false);
        ConformingTetrahedralMesh<2,2> conf_mesh;
        conf_mesh.AddNode(p_node);
        
        // Set up cells
        std::vector<MeinekeCryptCell> conf_cells;
        CellsGenerator<2>::GenerateForCrypt(conf_cells, conf_mesh, TYSONNOVAK, true);        
        Crypt<2> conf_crypt(conf_mesh, conf_cells);

        Crypt<2>::Iterator conf_iter = conf_crypt.Begin();

        TissueSimulation<2> simulator(conf_crypt);
        
        // repeat two times for coverage
        // need vector from parent to daughter to have both +ve and -ve y component
        // different branches will execute to make sure daughter stays in crypt ie. +ve y component
        for (unsigned repetitions=0; repetitions<=1; repetitions++)
        {                
            c_vector<double, 2> daughter_location = simulator.CalculateDividingCellCentreLocations(conf_iter);
            c_vector<double, 2> new_parent_location = conf_mesh.GetNode(0)->rGetLocation();
            c_vector<double, 2> parent_to_daughter = conf_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);
            // The parent stem cell should stay where it is and the daughter be introduced at positive y.
            
            TS_ASSERT_DELTA(new_parent_location[0], location[0], 1e-7);
            TS_ASSERT_DELTA(new_parent_location[1], location[1], 1e-7);
            TS_ASSERT(daughter_location[1]>=location[1]);
            TS_ASSERT_DELTA(norm_2(parent_to_daughter), 
                1.0*CancerParameters::Instance()->GetDivisionSeparation(),
                1e-7);
       }

        SimulationTime::Destroy();
    }

    void TestCalculateDividingCellCentreLocationsCylindricalMesh() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Make a mesh
        c_vector<double ,2> location;
        location[0]=1.0;
        location[1]=1.0;
        Node<2>* p_node = new Node<2>(0u,location, false);
        Cylindrical2dMesh cyl_mesh(6.0);
        cyl_mesh.AddNode(p_node);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cyl_cells;
        CellsGenerator<2>::GenerateForCrypt(cyl_cells, cyl_mesh, TYSONNOVAK, true);        
        Crypt<2> cyl_crypt(cyl_mesh, cyl_cells);

        Crypt<2>::Iterator cyl_iter = cyl_crypt.Begin();

        TissueSimulation<2> simulator(cyl_crypt);                
        c_vector<double, 2> daughter_location = simulator.CalculateDividingCellCentreLocations(cyl_iter);
        c_vector<double, 2> new_parent_location = cyl_mesh.GetNode(0)->rGetLocation();
        c_vector<double, 2> parent_to_daughter = cyl_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);
        TS_ASSERT_DELTA(norm_2(parent_to_daughter), 
            CancerParameters::Instance()->GetDivisionSeparation(),
            1e-7);

        SimulationTime::Destroy();
    }

    void TestCalculateDividingCellCentreLocationsCylindricalMeshStemCell() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Make a mesh
        c_vector<double ,2> location;
        location[0]=1.0;
        location[1]=0.0;                                    // <- y=0
        Node<2>* p_node = new Node<2>(0u,location, false);
        Cylindrical2dMesh cyl_mesh(6.0);
        cyl_mesh.AddNode(p_node);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cyl_cells;
        CellsGenerator<2>::GenerateForCrypt(cyl_cells, cyl_mesh, TYSONNOVAK, true);        
        Crypt<2> cyl_crypt(cyl_mesh, cyl_cells);

        Crypt<2>::Iterator cyl_iter = cyl_crypt.Begin();

        TissueSimulation<2> simulator(cyl_crypt);                
        c_vector<double, 2> daughter_location = simulator.CalculateDividingCellCentreLocations(cyl_iter);
        c_vector<double, 2> new_parent_location = cyl_mesh.GetNode(0)->rGetLocation();
        c_vector<double, 2> parent_to_daughter = cyl_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);
        
        // The parent stem cell should stay where it is and the daughter be introduced at positive y.
        TS_ASSERT_DELTA(new_parent_location[0], location[0], 1e-7);
        TS_ASSERT_DELTA(new_parent_location[1], location[1], 1e-7);
        TS_ASSERT(daughter_location[1]>=location[1]);
        TS_ASSERT_DELTA(norm_2(parent_to_daughter), 
            CancerParameters::Instance()->GetDivisionSeparation(),
            1e-7);

        SimulationTime::Destroy();
    }
    
    // short test which sets mNoBirth for coverage
    void TestNoBirth() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        
        std::string output_directory = "Crypt2DCylindricalNoBirth";        
        unsigned cells_across = 2;
        unsigned cells_up = 3;
        double crypt_width = 2.0;
        unsigned thickness_of_ghost_layer = 0;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);// true = mature cells

        Crypt<2> crypt(*p_mesh, cells);               
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);

        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(2.0); // long enough for a cell to be born were SetNoBirth not called
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        // These are for coverage and use the defaults
        simulator.SetDt(1.0/120.0);
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(true);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);

        simulator.Solve();

        // test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
                
        TS_ASSERT_EQUALS(number_of_cells, cells_across*cells_up); 
        TS_ASSERT_EQUALS(number_of_nodes, number_of_cells+thickness_of_ghost_layer*2*cells_across); 

        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    

    // Death on a non-periodic mesh
    // Massive amount of random death eventually leading to every cell being killed off..
    // Note that birth does occur too.
    void TestRandomDeathOnNonPeriodicCrypt() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        
        unsigned cells_across = 2;
        unsigned cells_up = 1;
        unsigned thickness_of_ghost_layer = 1;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
              
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DRandomDeathNonPeriodic");
        simulator.SetEndTime(0.5);
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);

        AbstractCellKiller<2>* p_random_cell_killer = new RandomCellKiller<2>(&crypt, 0.1);
        simulator.AddCellKiller(p_random_cell_killer);

        simulator.Solve();
        
        // there should be no cells left after this amount of time
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 0u);
    
        delete p_random_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    void TestUsingNonFlatBottomSurfaceNonPeriodic()
    {
        CancerParameters::Instance()->Reset();
        
        unsigned cells_across = 20;
        unsigned cells_up = 4;
        unsigned thickness_of_ghost_layer = 2;
        double crypt_width = 20.0;
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);
              
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DNonFlatBottomSurface");
        simulator.SetEndTime(0.5);
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);

        simulator.UseNonFlatBottomSurface();

        simulator.Solve();
        
        double frequency = floor(crypt_width/4.0)+1.0;
        TS_ASSERT_DELTA(p_mesh->GetNode(48)->rGetLocation()[0],8.0, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(48)->rGetLocation()[1],  0.05*(sin(2.0*frequency*M_PI*8.0/crypt_width) + 1.0 ), 1e-4);
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTCRYPTSIMULATION2D_HPP_*/
