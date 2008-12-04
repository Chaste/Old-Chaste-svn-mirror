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
#ifndef TESTCRYPTSIMULATION2D_HPP_
#define TESTCRYPTSIMULATION2D_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CryptSimulation2d.hpp"
#include "FixedCellCycleModelCellsGenerator.hpp"
#include "StochasticCellCycleModelCellsGenerator.hpp"
#include "WntCellCycleModelCellsGenerator.hpp"
#include "TysonNovakCellCycleModelCellsGenerator.hpp"
#include "IngeWntSwatCellCycleModelCellsGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractCancerTestSuite.hpp"
#include "../../global/test/NumericFileComparison.hpp"
#include "CancerEventHandler.hpp"

class TestCryptSimulation2d : public AbstractCancerTestSuite
{
private:

    /**
     * Compare two meshes to see if they are 'the same'.  Doesn't check everything,
     * but is fairly thorough.  Used for testing serialization.
     */
    template<unsigned DIM>
    void CompareMeshes(MutableMesh<DIM,DIM>* pMesh1,
                       MutableMesh<DIM,DIM>* pMesh2)
    {
        TS_ASSERT_EQUALS(pMesh1->GetNumAllNodes(), pMesh2->GetNumAllNodes());
        TS_ASSERT_EQUALS(pMesh1->GetNumNodes(), pMesh2->GetNumNodes());
        TS_ASSERT_EQUALS(pMesh1->GetNumBoundaryNodes(), pMesh2->GetNumBoundaryNodes());

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
        typename MutableMesh<DIM,DIM>::ElementIterator it=pMesh1->GetElementIteratorBegin();
        typename MutableMesh<DIM,DIM>::ElementIterator it2=pMesh2->GetElementIteratorBegin();
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

    void TestUpdatePositions() throw (Exception)
    {
        HoneycombMeshGenerator generator(3, 3, 1, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(tissue);

        std::vector<c_vector<double, 2> > old_posns(p_mesh->GetNumNodes());
        std::vector<c_vector<double, 2> > velocities_on_each_node(p_mesh->GetNumNodes());

        // Make some velocities up..
        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            old_posns[i][0] = p_mesh->GetNode(i)->rGetLocation()[0];
            old_posns[i][1] = p_mesh->GetNode(i)->rGetLocation()[1];

            velocities_on_each_node[i][0] = i*0.01;
            velocities_on_each_node[i][1] = 2*i*0.01;
       }

        simulator.SetDt(0.01);
        simulator.UpdateNodePositions(velocities_on_each_node);

        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            std::set<unsigned>::iterator iter = ghost_node_indices.find(i);
            bool is_a_ghost_node = (iter!=ghost_node_indices.end());
            Node<2>* p_node = p_mesh->GetNode(i);
            if (!is_a_ghost_node)
            {
                if(old_posns[i][1]==0) // stem
                {
                    // No wnt so shouldn't have been moved
                    TS_ASSERT_DELTA(p_node->rGetLocation()[0], old_posns[i][0], 1e-9);
                    TS_ASSERT_DELTA(p_node->rGetLocation()[1], old_posns[i][1], 1e-9);
                }
                else
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[0], old_posns[i][0] +   i*0.01*0.01, 1e-9);
                    TS_ASSERT_DELTA(p_node->rGetLocation()[1], old_posns[i][1] + 2*i*0.01*0.01, 1e-9);
                }
            }
        }
    }

    void Test2DCylindrical() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);// true = mature cells

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(crypt);

        simulator.SetEndTime(0.1);

        simulator.SetOutputCellMutationStates(true);

        // These are for coverage and use the defaults
        simulator.SetDt(1.0/120.0);
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(false);
        simulator.SetOutputDirectory("Crypt2DCylindrical");

        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        simulator.Solve();

        // Test we have the same number of cells and nodes at the end of each time
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
        CancerEventHandler::Reset(); // otherwise event handler left in bad state after throw
    }


    void Test2DCylindricalMultipleDivisions() throw (Exception)
    {
        // Create a log of this test.
        LogFile* p_log_file = LogFile::Instance();
        p_log_file->Set(2,"Crypt2DCylindricalMultipleDivisions");

        unsigned cells_across = 6;
        unsigned cells_up = 8;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);// true = mature cells

        for (unsigned i=0; i<cells.size();i++)
        {
            cells[i].SetBirthTime(-11.5);
        }

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(crypt);

        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        simulator.SetOutputDirectory("Crypt2DCylindricalMultipleDivisions");
        simulator.SetEndTime(0.6);
        simulator.Solve();

        // Find the height of the current crypt
        double height_after_division=p_mesh->GetWidth(1);
        simulator.SetEndTime(0.8);
        simulator.Solve();

        // Find the height of the current crypt
        double height_after_relaxation=p_mesh->GetWidth(1);

        TS_ASSERT_LESS_THAN(height_after_division, height_after_relaxation);

        simulator.SetEndTime(2.0);
        simulator.Solve();

        // All fully differentiated cells have sloughed off
        for (MeshBasedTissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
             TS_ASSERT(cell_iter->GetCellType() != DIFFERENTIATED);
        }

        // Close the log file opened in this test
        LogFile::Close();
    }

    // This is a rubbish test - all cells start at birthTime = 0.
    // So bizarrely the crypt shrinks as the rest lengths are shortened!
    // But at least it uses Wnt cell cycle and runs reasonably quickly...
    // For a better test with more randomly distributed cell ages see the Nightly test pack.
    void TestWithWntDependentCells() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 0;

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        WntCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, false);

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(crypt);

        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DPeriodicWnt");

        // Set length of simulation here
        simulator.SetEndTime(0.3);

        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        simulator.Solve();

        std::vector<double> node_35_location = simulator.GetNodeLocation(35);

        TS_ASSERT_DELTA(node_35_location[0], 5.5000 , 1e-4);

        // Old version of this test had cells with age zero, therefore small spring lengths.
        // Variable spring lengths now only associated with cell division.
        TS_ASSERT_DELTA(node_35_location[1], 4.33013 , 1e-4);

        WntConcentration::Destroy();
    }

    // A better check that the loaded mesh is the same as that saved
    void TestMeshSurvivesSaveLoad() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Set up a simulation
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        std::vector<TissueCell> cells;

        WntCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, false);
        
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(crypt);

        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DMeshArchive");
        simulator.SetEndTime(0.1);

        // Memory leak (unconditional jump) without the following line.
        // The archiver assumes that a Solve has been called and simulation time has been set up properly.
        // In this test it hasn't so we need this to avoid memory leak.
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.1, 100);

        // Save
        simulator.Save();

        // Load
        CryptSimulation2d* p_simulator;
        p_simulator = CryptSimulation2d::Load("Crypt2DMeshArchive", 0.0);

        // Create an identical mesh for comparison purposes
        HoneycombMeshGenerator generator2(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh2=generator2.GetCylindricalMesh();

        // Compare
        MutableMesh<2,2>& r_mesh = (static_cast<MeshBasedTissue<2>*>(&(p_simulator->rGetTissue())))->rGetMesh();
        CompareMeshes(p_mesh2, &r_mesh);

        delete p_simulator;
        WntConcentration::Destroy();
    }

    // A check that save and load works when a Voronoi tessellation is involved
    void TestMeshSurvivesSaveLoadWithVoroniTessellation() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Set up a simulation
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        std::vector<TissueCell> cells;

        WntCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, false);
        
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(crypt);

        Meineke2001SpringSystem<2>* p_meineke_spring_system = new Meineke2001SpringSystem<2>(crypt);
        p_meineke_spring_system->SetAreaBasedViscosity(true);
        p_meineke_spring_system->SetEdgeBasedSpringConstant(true);

        CryptSimulation2d simulator(crypt, p_meineke_spring_system, false, true);
        simulator.SetOutputDirectory("Crypt2DMeshArchive2");
        simulator.SetEndTime(0.1);

        simulator.Solve();

        // Save
        simulator.Save();

        // Load
        CryptSimulation2d* p_simulator;
        p_simulator = CryptSimulation2d::Load("Crypt2DMeshArchive2", 0.1);
        p_simulator->SetEndTime(0.15);
        p_simulator->Solve();

        // Tidy up
        delete p_simulator;
        delete p_meineke_spring_system;
        WntConcentration::Destroy();
    }


    void TestStandardResultForArchivingTestsBelow() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        // We have a Wnt Gradient - but not Wnt dependent cells
        // so that the test runs quickly, but we test archiving of it!
        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(crypt);

        CryptSimulation2d simulator(crypt);

        simulator.SetOutputDirectory("Crypt2DPeriodicStandardResult");

        // Our full end time is 0.25
        simulator.SetEndTime(0.25);

        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);
        simulator.SetOutputCellCyclePhases(true);// For coverage only...
        simulator.Solve();

        // These cells just divided and have been gradually moving apart.
        // These results are from time 0.25, which is also tested below
        // after a save and a load. (See #420, #479.)
        std::vector<double> node_28_location = simulator.GetNodeLocation(28);
        TS_ASSERT_DELTA(node_28_location[0], 3.7875 , 1e-4);
        TS_ASSERT_DELTA(node_28_location[1], 0.0 , 1e-4);
        std::vector<double> node_120_location = simulator.GetNodeLocation(120);
        TS_ASSERT_DELTA(node_120_location[0], 4.2035 , 1e-4);
        TS_ASSERT_DELTA(node_120_location[1], 0.1033 , 1e-4);

        // Test the Wnt gradient result
        TissueCell* p_cell = &(crypt.rGetCellAtNodeIndex(28));
        TS_ASSERT_DELTA(WntConcentration::Instance()->GetWntLevel(p_cell), 1.0, 1e-9);
        p_cell = &(crypt.rGetCellAtNodeIndex(120));
        TS_ASSERT_DELTA(WntConcentration::Instance()->GetWntLevel(p_cell), 0.9900, 1e-4);
        WntConcentration::Destroy();

        // Check writing of voronoi data
        OutputFileHandler handler("Crypt2DPeriodicStandardResult",false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellcyclephases.dat";

        NumericFileComparison comp(results_file,"cancer/test/data/CellCyclePhaseOutput/cellcyclephases.dat");
        TS_ASSERT(comp.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_file + " cancer/test/data/CellCyclePhaseOutput/cellcyclephases.dat").c_str()), 0);
    }

    // Testing Save
    void TestSave() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(crypt);

        CryptSimulation2d simulator(crypt);

        simulator.SetOutputDirectory("Crypt2DPeriodicSaveAndLoad");

        // Our full end time is 0.25, here we run until 0.1 then load and run more below.
        simulator.SetEndTime(0.1);

        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        simulator.Solve();

        // Save the results..
        simulator.Save();

        WntConcentration::Destroy();
    }


    // Testing Load (based on previous two tests)
    void TestLoad() throw (Exception)
    {
        // Load the simulation from the TestSave method above and
        // run it from 0.1 to 0.2
        CryptSimulation2d* p_simulator1;

        WntConcentration::Instance();   // Make sure there is no existing Wnt Gradient before load.
        WntConcentration::Destroy();

        p_simulator1 = CryptSimulation2d::Load("Crypt2DPeriodicSaveAndLoad", 0.1);

        p_simulator1->SetEndTime(0.2);
        p_simulator1->Solve();

        //Save that then reload
        // and run from 0.2 to 0.25.
        NodeMap map(0);

        MutableMesh<2,2>& r_mesh1 = (static_cast<MeshBasedTissue<2>*>(&(p_simulator1->rGetTissue())))->rGetMesh();
        r_mesh1.ReMesh(map);
        p_simulator1->Save();

        CryptSimulation2d* p_simulator2 = CryptSimulation2d::Load("Crypt2DPeriodicSaveAndLoad", 0.2);

        MutableMesh<2,2>& r_mesh2 = (static_cast<MeshBasedTissue<2>*>(&(p_simulator2->rGetTissue())))->rGetMesh();

        CompareMeshes(&r_mesh1, &r_mesh2);

        p_simulator2->SetEndTime(0.25);
        p_simulator2->Solve();

        // These cells just divided and have been gradually moving apart.
        // These results are from time 0.25 in the StandardResult test above.
        std::vector<double> node_28_location = p_simulator2->GetNodeLocation(28);
        TS_ASSERT_DELTA(node_28_location[0], 3.7875 , 1e-4);
        TS_ASSERT_DELTA(node_28_location[1], 0.0 , 1e-4);
        std::vector<double> node_120_location = p_simulator2->GetNodeLocation(120);
        TS_ASSERT_DELTA(node_120_location[0], 4.2035 , 1e-4);
        TS_ASSERT_DELTA(node_120_location[1], 0.1033 , 1e-4);

        // Test Wnt Gradient was set up correctly
        TS_ASSERT_EQUALS(WntConcentration::Instance()->IsWntSetUp(),true);

        // Test the Wnt gradient result
        TissueCell* p_cell = &(p_simulator2->rGetTissue().rGetCellAtNodeIndex(28));
        TS_ASSERT_DELTA(WntConcentration::Instance()->GetWntLevel(p_cell), 1.0, 1e-9);
        p_cell = &(p_simulator2->rGetTissue().rGetCellAtNodeIndex(120));
        TS_ASSERT_DELTA(WntConcentration::Instance()->GetWntLevel(p_cell), 0.9900, 1e-4);

        delete p_simulator1;
        delete p_simulator2;
        WntConcentration::Destroy();
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
        unsigned cells_across = 2;
        unsigned cells_up = 2;
        double crypt_width = 0.5;
        unsigned thickness_of_ghost_layer = 1;

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, false, crypt_width/cells_across);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        WntCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        // Cover the write Voronoi data method
        crypt.SetWriteVoronoiData(true, false);

        MeshBasedTissue<2>::Iterator cell_iterator = crypt.Begin();
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

        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(crypt);

        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DWntMatureCells");

        // If you want to visualize this use the 'notcylindrical' option
        // (it is too small for it to figure out what's happening on its own).

        simulator.SetEndTime(0.01);

        // cover exceptions
        TS_ASSERT_THROWS_ANYTHING(simulator.GetCellMutationStateCount());
        simulator.SetOutputCellMutationStates(true);

        TS_ASSERT_THROWS_ANYTHING(simulator.GetCellTypeCount());
        simulator.SetOutputCellTypes(true);

        TS_ASSERT_THROWS_ANYTHING(simulator.GetCellCyclePhaseCount());

        simulator.Solve();

        // Check that nothing has moved below y=0
        for (MeshBasedTissue<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            TS_ASSERT_LESS_THAN(-1e-15,p_mesh->GetNode(cell_iter->GetNodeIndex())->rGetLocation()[1]);
        }

        c_vector<unsigned,5> cell_mutation_state_count = simulator.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_state_count[0], 1u);
        TS_ASSERT_EQUALS(cell_mutation_state_count[1], 1u);
        TS_ASSERT_EQUALS(cell_mutation_state_count[2], 1u);
        TS_ASSERT_EQUALS(cell_mutation_state_count[3], 0u);  // No APC two hit, one of all the rest.
        TS_ASSERT_EQUALS(cell_mutation_state_count[4], 1u);

        c_vector<unsigned,5> cell_type_count = simulator.GetCellTypeCount();
        TS_ASSERT_EQUALS(cell_type_count[0], 0u);
        TS_ASSERT_EQUALS(cell_type_count[1], 4u);
        TS_ASSERT_EQUALS(cell_type_count[2], 0u);
        TS_ASSERT_EQUALS(cell_type_count[3], 0u);

        // Check writing of voronoi data
        OutputFileHandler handler("Crypt2DWntMatureCells",false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizvoronoi";

        NumericFileComparison comp(results_file,"cancer/test/data/Crypt2DWntMatureCells/VoronoiAreaAndPerimeter.dat");
        TS_ASSERT(comp.CompareFiles(2e-6));
        //Cover writing logged cell
        crypt.SetWriteVoronoiData(true, true);
        simulator.SetEndTime(0.01 + 1./120.);
        // Set the first cell to be logged
        crypt.Begin()->SetLogged();

        simulator.Solve();

        WntConcentration::Destroy();
    }

    // This is a strange test -- all cells divide within a quick time, it gives
    // good testing of the periodic boundaries though... [comment no longer valid?]
    void TestWithTysonNovakCells() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        TysonNovakCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(crypt);

        SloughingCellKiller sloughing_cell_killer(&crypt);
        simulator.AddCellKiller(&sloughing_cell_killer);

        simulator.SetOutputDirectory("Crypt2DPeriodicTysonNovak");
        simulator.SetEndTime(0.05);
        simulator.SetDt(0.001);

        // Test that labelling a few cells doesn't make any difference to the simulation
        // and therefore log them in the visualizer files for the next test to check.
        simulator.rGetTissue().rGetCellAtNodeIndex(57).SetMutationState(LABELLED);
        simulator.rGetTissue().rGetCellAtNodeIndex(56).SetMutationState(APC_ONE_HIT);
        simulator.rGetTissue().rGetCellAtNodeIndex(51).SetMutationState(APC_TWO_HIT);
        simulator.rGetTissue().rGetCellAtNodeIndex(63).SetMutationState(BETA_CATENIN_ONE_HIT);
        simulator.SetOutputCellMutationStates(true);
        simulator.Solve();

        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<bool> ghost_cells = crypt.rGetGhostNodes();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();

        TS_ASSERT_EQUALS(number_of_nodes,ghost_cells.size());
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 75u);
        TS_ASSERT_EQUALS(number_of_nodes, 123u);
    }

    /*
     * This test compares the visualizer output from the previous test with a known file.
     *
     * The results of this should be a yellow crypt with a line of
     * blue stem cells at the base, and four labelled cells of different
     * colours in the centre. A couple of cells divide, the crypt stays
     * periodic and a couple of cells swap sides.
     *
     * Note - if the previous test is changed we need to update the file this test refers to.
     */
    void TestVisualizerOutput() throw (Exception)
    {
        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt2DPeriodicTysonNovak",false);
        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0";

        NumericFileComparison comp_ele(results_dir + "/results.vizelements","cancer/test/data/Crypt2DPeriodicTysonNovak_vis/results.vizelements");
        TS_ASSERT(comp_ele.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizelements cancer/test/data/Crypt2DPeriodicTysonNovak_vis/results.vizelements").c_str()), 0);

        NumericFileComparison comp_nodes(results_dir + "/results.viznodes","cancer/test/data/Crypt2DPeriodicTysonNovak_vis/results.viznodes");
        TS_ASSERT(comp_nodes.CompareFiles(1e-15));
        //TS_ASSERT_EQUALS(system(("ndiff -abserr 1e-14 " + results_dir + "/results.viznodes cancer/test/data/Crypt2DPeriodicTysonNovak_vis/results.viznodes").c_str()), 0);

        //NumericFileComparison comp_setup(results_dir + "/results.vizsetup","cancer/test/data/Crypt2DPeriodicTysonNovak_vis/results.vizsetup");
        //TS_ASSERT(comp_setup.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizsetup cancer/test/data/Crypt2DPeriodicTysonNovak_vis/results.vizsetup").c_str()), 0);
    }


    void TestAddCellKiller() throw (Exception)
    {
        double crypt_length = 9.3;
        double crypt_width = 10.0;

        CancerParameters::Instance()->SetCryptLength(crypt_length);
        CancerParameters::Instance()->SetCryptWidth(crypt_width);

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, mesh, false, 0.0, 3.0, 6.5, 8.0);

        cells[60].SetBirthTime(-50.0);

        MeshBasedTissueWithGhostNodes<2> crypt(mesh, cells);

        CryptSimulation2d simulator(crypt);

        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        unsigned num_deaths = simulator.DoCellRemoval();
        unsigned num_births = simulator.DoCellBirth();

        TS_ASSERT_EQUALS(num_births, 1u);
        TS_ASSERT_EQUALS(num_deaths,11u);
    }

    void TestCalculateDividingCellCentreLocationsConfMesh() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        CancerParameters::Instance()->SetDivisionRestingSpringLength(0.9);//Only coverage
        CancerParameters::Instance()->SetDivisionSeparation(0.1);

        // Make a parent node
        c_vector<double,2> location;
        location[0] = 1.0;
        location[1] = 1.0;
        Node<2>* p_node = new Node<2>(0u,location, false);

        MutableMesh<2,2> conf_mesh;
        conf_mesh.AddNode(p_node);

        // Set up cells
        std::vector<TissueCell> conf_cells;
        TysonNovakCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(conf_cells, conf_mesh, true);

        MeshBasedTissue<2> conf_crypt(conf_mesh, conf_cells);

        MeshBasedTissue<2>::Iterator conf_iter = conf_crypt.Begin();

        CryptSimulation2d simulator(conf_crypt);

        c_vector<double, 2> daughter_location = simulator.CalculateDividingCellCentreLocations(conf_iter);
        c_vector<double, 2> new_parent_location = conf_mesh.GetNode(0)->rGetLocation();
        c_vector<double, 2> parent_to_daughter = conf_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);
        TS_ASSERT_DELTA(norm_2(parent_to_daughter),
            CancerParameters::Instance()->GetDivisionSeparation(),
            1e-7);
    }


    void TestCalculateDividingCellCentreLocationsConfMeshStemCell() throw (Exception)
    {
        // Make a parent node
        c_vector<double,2> location;
        location[0] = 1.0;
        location[1] = 0.0; // <- y=0
        Node<2>* p_node = new Node<2>(0u,location, false);
        MutableMesh<2,2> conf_mesh;
        conf_mesh.AddNode(p_node);

        // Set up cells
        std::vector<TissueCell> conf_cells;
        TysonNovakCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(conf_cells, conf_mesh, true);
        MeshBasedTissue<2> conf_crypt(conf_mesh, conf_cells);

        MeshBasedTissue<2>::Iterator conf_iter = conf_crypt.Begin();

        CryptSimulation2d simulator(conf_crypt);

        // Repeat two times for coverage
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
    }

    void TestCalculateDividingCellCentreLocationsCylindricalMesh() throw (Exception)
    {
        // Make a mesh
        c_vector<double,2> location;
        location[0] = 1.0;
        location[1] = 1.0;
        Node<2>* p_node = new Node<2>(0u,location, false);
        Cylindrical2dMesh cyl_mesh(6.0);
        cyl_mesh.AddNode(p_node);

        // Set up cells
        std::vector<TissueCell> cyl_cells;
        TysonNovakCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cyl_cells, cyl_mesh, true);
        MeshBasedTissue<2> cyl_crypt(cyl_mesh, cyl_cells);

        MeshBasedTissue<2>::Iterator cyl_iter = cyl_crypt.Begin();

        CryptSimulation2d simulator(cyl_crypt);
        c_vector<double, 2> daughter_location = simulator.CalculateDividingCellCentreLocations(cyl_iter);
        c_vector<double, 2> new_parent_location = cyl_mesh.GetNode(0)->rGetLocation();
        c_vector<double, 2> parent_to_daughter = cyl_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);
        TS_ASSERT_DELTA(norm_2(parent_to_daughter),
            CancerParameters::Instance()->GetDivisionSeparation(),
            1e-7);
    }

    void TestCalculateDividingCellCentreLocationsCylindricalMeshStemCell() throw (Exception)
    {
        // Make a mesh
        c_vector<double,2> location;
        location[0] = 1.0;
        location[1] = 0.0; // <- y=0
        Node<2>* p_node = new Node<2>(0u,location, false);
        Cylindrical2dMesh cyl_mesh(6.0);
        cyl_mesh.AddNode(p_node);

        // Set up cells
        std::vector<TissueCell> cyl_cells;
        TysonNovakCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cyl_cells, cyl_mesh, true);
        
        MeshBasedTissue<2> cyl_crypt(cyl_mesh, cyl_cells);

        MeshBasedTissue<2>::Iterator cyl_iter = cyl_crypt.Begin();

        CryptSimulation2d simulator(cyl_crypt);
        c_vector<double,2> daughter_location = simulator.CalculateDividingCellCentreLocations(cyl_iter);
        c_vector<double,2> new_parent_location = cyl_mesh.GetNode(0)->rGetLocation();
        c_vector<double,2> parent_to_daughter = cyl_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);

        // The parent stem cell should stay where it is and the daughter be introduced at positive y.
        TS_ASSERT_DELTA(new_parent_location[0], location[0], 1e-7);
        TS_ASSERT_DELTA(new_parent_location[1], location[1], 1e-7);
        TS_ASSERT(daughter_location[1]>=location[1]);
        TS_ASSERT_DELTA(norm_2(parent_to_daughter),
            CancerParameters::Instance()->GetDivisionSeparation(),
            1e-7);
    }

    // Short test which sets mNoBirth for coverage
    void TestNoBirth() throw (Exception)
    {
        std::string output_directory = "Crypt2DCylindricalNoBirth";
        unsigned cells_across = 2;
        unsigned cells_up = 3;
        double crypt_width = 2.0;
        unsigned thickness_of_ghost_layer = 0;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);// true = mature cells

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(crypt);

        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(2.0); // long enough for a cell to be born were SetNoBirth not called

        // These are for coverage and use the defaults
        simulator.SetDt(1.0/120.0);
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(true);

        SloughingCellKiller sloughing_cell_killer(&crypt, true);
        simulator.AddCellKiller(&sloughing_cell_killer);

        simulator.Solve();

        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();

        TS_ASSERT_EQUALS(number_of_cells, cells_across*cells_up);
        TS_ASSERT_EQUALS(number_of_nodes, number_of_cells+thickness_of_ghost_layer*2*cells_across);
    }


    // Death on a non-periodic mesh
    // Massive amount of random death eventually leading to every cell being killed off..
    // Note that birth does occur too.
    void TestRandomDeathOnNonPeriodicCrypt() throw (Exception)
    {
        unsigned cells_across = 2;
        unsigned cells_up = 1;
        unsigned thickness_of_ghost_layer = 1;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(crypt);

        simulator.SetOutputDirectory("Crypt2DRandomDeathNonPeriodic");
        simulator.SetEndTime(0.6);

        RandomCellKiller<2> random_cell_killer(&crypt, 0.1);
        simulator.AddCellKiller(&random_cell_killer);

        simulator.Solve();

        // There should be no cells left after this amount of time
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 0u);
    }

    void TestUsingJiggledBottomSurface()
    {
        HoneycombMeshGenerator generator(4, 4, 0, true, 1.0);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true);

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        CryptSimulation2d simulator(crypt);

        simulator.SetOutputDirectory("Crypt2DJiggledBottomCells");
        simulator.SetEndTime(0.01);
        simulator.UseJiggledBottomCells();

        //Move the first cell (which should be on y=0) down a bit
        MeshBasedTissue<2>::Iterator cell_iter = crypt.Begin();
        assert(cell_iter.rGetLocation()[1] == 0.0);

        // Move the cell (can't use the iterator for this as it is const)
        crypt.rGetMesh().GetNode(0)->rGetModifiableLocation()[1] = -0.1;
        assert(cell_iter.rGetLocation()[1] < 0.0);

        simulator.Solve();

        //The cell should have been pulled up, but not above y=0. However it should
        // then been moved to above y=0 by the jiggling
        TS_ASSERT_LESS_THAN(0.0, cell_iter.rGetLocation()[1]);
    }

    void TestWriteBetaCatenin() throw (Exception)
    {
        HoneycombMeshGenerator generator(5, 4, 1);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        IngeWntSwatCellCycleModelCellsGenerator<2> cells_generator(1u);
        cells_generator.GenerateForCrypt(cells, *p_mesh, false);

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(crypt);

        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("CryptBetaCatenin");
        simulator.SetEndTime(0.01);

        simulator.Solve();

        // Check writing of beta-catenin data
        OutputFileHandler handler("CryptBetaCatenin",false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizbCat";
        std::string results_setup_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizsetup";

        NumericFileComparison comp_bcat(results_file,"cancer/test/data/CryptBetaCatenin/results.vizbCat");
        TS_ASSERT(comp_bcat.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_file + " cancer/test/data/CryptBetaCatenin/results.vizbCat").c_str()), 0);

        TS_ASSERT_EQUALS(system(("diff " + results_setup_file + " cancer/test/data/CryptBetaCatenin/results.vizsetup").c_str()), 0);

        WntConcentration::Destroy();
    }

    void TestApoptosisSpringLengths() throw (Exception)
    {
        unsigned num_cells_depth = 2;
        unsigned num_cells_width = 2;
        double crypt_length = num_cells_depth-0.0;
        double crypt_width = num_cells_width-0.0;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        CancerParameters* p_params = CancerParameters::Instance();
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);

        // Set up cells
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(TRANSIT, HEALTHY, new FixedCellCycleModel());
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*(p_params->GetTransitCellG1Duration()
                                               +p_params->GetSG2MDuration());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, ghost_node_indices);

        TissueSimulation<2> simulator(tissue);

        simulator.SetOutputDirectory("2dSpheroidApoptosis");
        simulator.SetEndTime(1.0);

        CancerParameters::Instance()->SetApoptosisTime(2.0);
        tissue.rGetCellAtNodeIndex(14).StartApoptosis();
        tissue.rGetCellAtNodeIndex(15).StartApoptosis();
        simulator.SetNoBirth(true);

        simulator.Solve();

        /* We track the locations of two dying cells (a and b) and two
         * live cells adjacent to them (c and d)
         *
         * All cells begin distance 1 apart.
         *
         * a and b move together to leave a gap of 0.
         * a and c (and b and d) move to a distance of 0.5 apart.
         */

        c_vector<double, 2> a_location = tissue.rGetMesh().GetNode(14)->rGetLocation();
        c_vector<double, 2> b_location = tissue.rGetMesh().GetNode(15)->rGetLocation();
        c_vector<double, 2> c_location = tissue.rGetMesh().GetNode(20)->rGetLocation();
        c_vector<double, 2> d_location = tissue.rGetMesh().GetNode(21)->rGetLocation();

        double a_b_separation = sqrt((a_location[0]-b_location[0])*(a_location[0]-b_location[0]) +
                                (a_location[1]-b_location[1])*(a_location[1]-b_location[1]));
        double a_c_separation = sqrt((a_location[0]-c_location[0])*(a_location[0]-c_location[0]) +
                                (a_location[1]-c_location[1])*(a_location[1]-c_location[1]));
        double c_d_separation = sqrt((d_location[0]-c_location[0])*(d_location[0]-c_location[0]) +
                                (d_location[1]-c_location[1])*(d_location[1]-c_location[1]));

        TS_ASSERT_DELTA(a_b_separation , 0.5, 1e-1);
        TS_ASSERT_DELTA(a_c_separation , 0.75, 1e-1);
        TS_ASSERT_DELTA(c_d_separation , 1.0, 1e-1);

        simulator.SetEndTime(1.99);
        simulator.Solve();

        a_location = tissue.rGetMesh().GetNode(14)->rGetLocation();
        b_location = tissue.rGetMesh().GetNode(15)->rGetLocation();
        c_location = tissue.rGetMesh().GetNode(20)->rGetLocation();
        d_location = tissue.rGetMesh().GetNode(21)->rGetLocation();

        a_b_separation = sqrt((a_location[0]-b_location[0])*(a_location[0]-b_location[0]) +
                         (a_location[1]-b_location[1])*(a_location[1]-b_location[1]));
        a_c_separation = sqrt((a_location[0]-c_location[0])*(a_location[0]-c_location[0]) +
                         (a_location[1]-c_location[1])*(a_location[1]-c_location[1]));
        c_d_separation = sqrt((d_location[0]-c_location[0])*(d_location[0]-c_location[0]) +
                         (d_location[1]-c_location[1])*(d_location[1]-c_location[1]));

        TS_ASSERT_DELTA(a_b_separation , 0.01, 1e-1);
        TS_ASSERT_DELTA(a_c_separation , 0.5, 1e-1);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 4u);

        simulator.SetEndTime(2.01);
        simulator.Solve();

        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 2u);
    }

    void TestAncestorCryptSimulations() throw (Exception)
    {
        std::string output_directory = "AncestorCrypt";

        unsigned cells_across = 13;
        unsigned cells_up = 25;
        double crypt_width = 12.1;
        unsigned thickness_of_ghost_layer = 3;

        CancerParameters* p_params = CancerParameters::Instance();

        p_params->SetDampingConstantNormal(1.0); // normally 1

        // Do not give mutant cells any different movement properties to normal ones
        p_params->SetDampingConstantMutant(p_params->GetDampingConstantNormal());
        p_params->SetSpringStiffness(30.0); //normally 15.0;

        std::set<unsigned> ghost_node_indices;
        ghost_node_indices.clear();

        std::vector<TissueCell> cells;

        double time_of_each_run;
        std::vector<bool> labelled;

        HoneycombMeshGenerator generator = HoneycombMeshGenerator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        ghost_node_indices = generator.GetGhostNodeIndices();

        Cylindrical2dMesh* p_mesh;
        p_mesh = generator.GetCylindricalMesh();

        // Set up cells
        StochasticCellCycleModelCellsGenerator<2>cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, true,
                                            0.3,2.0,3.0,4.0,true);
        // Set up crypt
        MeshBasedTissueWithGhostNodes<2>* p_crypt = new MeshBasedTissueWithGhostNodes<2>(*p_mesh, cells, ghost_node_indices);

        p_crypt->SetBottomCellAncestors();


        // Set up crypt simulation
        CryptSimulation2d simulator(*p_crypt, NULL, false, false);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputCellAncestors(true);

        // Set simulation to output cell types
        simulator.SetOutputCellMutationStates(true);

        // Set length of simulation here
        time_of_each_run = 10.0*simulator.GetDt(); // for each run
        simulator.SetEndTime(time_of_each_run);

        // Set up cell killer
        SloughingCellKiller killer(&(simulator.rGetTissue()),0.01);
        simulator.AddCellKiller(&killer);

        simulator.UseJiggledBottomCells();

        // Run for a bit
        simulator.Solve();

        WntConcentration::Destroy();

        delete p_crypt;

        // ... and checking visualization of labelled cells against previous run
        OutputFileHandler handler("AncestorCrypt",false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizAncestors";
        NumericFileComparison comp_nodes(results_file1,"cancer/test/data/AncestorCrypt/results.viznodes");
        TS_ASSERT(comp_nodes.CompareFiles());
        NumericFileComparison comp_ans(results_file2,"cancer/test/data/AncestorCrypt/results.vizAncestors");
        TS_ASSERT(comp_ans.CompareFiles());

        TS_ASSERT_EQUALS(system(("diff " + results_file1 + " cancer/test/data/AncestorCrypt/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_file2 + " cancer/test/data/AncestorCrypt/results.vizAncestors").c_str()), 0);

        delete p_params;
    }

};

#endif /*TESTCRYPTSIMULATION2D_HPP_*/
