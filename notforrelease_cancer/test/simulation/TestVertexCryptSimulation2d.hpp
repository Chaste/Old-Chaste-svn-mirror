/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef TESTVERTEXCRYPTSIMULATION2D_HPP_
#define TESTVERTEXCRYPTSIMULATION2D_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cancer headers
#include "TissueSimulationArchiver.hpp"

#include "VertexCryptSimulation2d.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "NagaiHondaForce.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "SloughingCellKiller.hpp"
#include "AbstractCancerTestSuite.hpp"
#include "CancerEventHandler.hpp"

class TestVertexCryptSimulation2d : public AbstractCancerTestSuite
{
private:

    /**
     * Compare two meshes to see if they are 'the same'.  Doesn't check everything,
     * but is fairly thorough.  Used for testing serialization.
     */
    template<unsigned DIM>
    void CompareMeshes(VertexMesh<DIM,DIM>* pMesh1, VertexMesh<DIM,DIM>* pMesh2)
    {
        TS_ASSERT_EQUALS(pMesh1->GetNumNodes(), pMesh2->GetNumNodes());

        for (unsigned i=0; i<pMesh1->GetNumNodes(); i++)
        {
            Node<DIM>* p_node1 = pMesh1->GetNode(i);
            Node<DIM>* p_node2 = pMesh2->GetNode(i);

            TS_ASSERT_EQUALS(p_node1->IsDeleted(), p_node2->IsDeleted());
            TS_ASSERT_EQUALS(p_node1->GetIndex(), p_node2->GetIndex());
            TS_ASSERT_EQUALS(p_node1->IsBoundaryNode(), p_node2->IsBoundaryNode());

            for (unsigned j=0; j<DIM; j++)
            {
                TS_ASSERT_DELTA(p_node1->rGetLocation()[j], p_node2->rGetLocation()[j], 1e-4);
            }
        }

        TS_ASSERT_EQUALS(pMesh1->GetNumElements(), pMesh2->GetNumElements());

		for (typename VertexMesh<DIM,DIM>::VertexElementIterator iter = pMesh1->GetElementIteratorBegin();
             iter != pMesh1->GetElementIteratorEnd();
             ++iter)
		{
            unsigned elem_index = iter->GetIndex();
            VertexElement<DIM,DIM>* p_elt2 = pMesh2->GetElement(elem_index);
            TS_ASSERT_EQUALS(iter->GetNumNodes(), p_elt2->GetNumNodes());

            for (unsigned j=0; j<iter->GetNumNodes(); j++)
            {
                TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(j), p_elt2->GetNodeGlobalIndex(j));
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
        // Create mesh
        Cylindrical2dVertexMesh mesh(6, 6, 0.01, 2.0);

        // Set parameters
        TissueConfig::Instance()->SetMaxTransitGenerations(UINT_MAX);

        // Create cells
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<mesh.GetNumElements(); elem_index++)
        {
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                ( TissueConfig::Instance()->GetTransitCellG1Duration()
                                    + TissueConfig::Instance()->GetSG2MDuration() );

            TissueCell cell(TRANSIT, HEALTHY, new StochasticDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(tissue, force_collection);

        std::vector<c_vector<double, 2> > old_node_locations(mesh.GetNumNodes());
        std::vector<c_vector<double, 2> > forces(mesh.GetNumNodes());

        // Make up some forces
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            old_node_locations[i][0] = mesh.GetNode(i)->rGetLocation()[0];
            old_node_locations[i][1] = mesh.GetNode(i)->rGetLocation()[1];

            forces[i][0] = i*0.01;
            forces[i][1] = 2*i*0.01;
       }

        simulator.SetDt(0.01);
        simulator.UpdateNodePositions(forces);

        for (unsigned node_index=0; node_index<simulator.rGetTissue().GetNumNodes(); node_index++)
        {
            c_vector<double, 2> node_location = simulator.rGetTissue().GetNode(node_index)->rGetLocation();

            TS_ASSERT_DELTA(node_location[0], old_node_locations[node_index][0] +   node_index*0.01*0.01, 1e-9);
            TS_ASSERT_DELTA(node_location[1], old_node_locations[node_index][1] + 2*node_index*0.01*0.01, 1e-9);
        }
    }


    void TestCryptWithNoBirth() throw (Exception)
    {
        // Create mesh
        Cylindrical2dVertexMesh mesh(6, 12, 0.01, 2.0);

        // Create cells
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<mesh.GetNumElements(); elem_index++)
        {
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                ( TissueConfig::Instance()->GetTransitCellG1Duration()
                                    + TissueConfig::Instance()->GetSG2MDuration() );

            TissueCell cell(DIFFERENTIATED, HEALTHY, new StochasticDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> crypt(mesh, cells);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetEndTime(1.0);
        simulator.SetOutputDirectory("TestVertexCryptWithNoBirth");

        SloughingCellKiller<2> sloughing_cell_killer(&crypt, false);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }

    void TestCryptWithBirth() throw (Exception)
    {
        // Create mesh
        Cylindrical2dVertexMesh mesh(4, 6, 0.01, 2.0, true);

        // Create cells
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<mesh.GetNumElements(); elem_index++)
        {
            double birth_time = - RandomNumberGenerator::Instance()->ranf()*
                                 ( TissueConfig::Instance()->GetTransitCellG1Duration()
                                    + TissueConfig::Instance()->GetSG2MDuration() );

            CellType cell_type;

            // Cell 1 should divide at time t=0.5
            if (elem_index==0)
            {
                birth_time = -23.5;
                cell_type = STEM;
            }
            // Cells 2 3 and 4 should divide at later times 
            else if ((elem_index==1)||(elem_index==2)||(elem_index==3))
            {
                birth_time = -15.5 - 2.0*(double)elem_index;
                cell_type = STEM;
            }
            else
            {
                cell_type = DIFFERENTIATED;
            }

            TissueCell cell(cell_type, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> crypt(mesh, cells);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(1.0);
        simulator.SetOutputDirectory("TestVertexCryptWithBirth");
        
        // Make crypt shorter for sloughing 
        TissueConfig::Instance()->SetCryptLength(5.0);
                
        SloughingCellKiller<2> sloughing_cell_killer(&crypt);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }
    
    void noTestCryptSimulationLong() throw (Exception)
    {
        // Create mesh
        Cylindrical2dVertexMesh mesh(4, 8, 0.01, 2.0, true);

        // Create cells
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<mesh.GetNumElements(); elem_index++)
        {
            double birth_time = - RandomNumberGenerator::Instance()->ranf()*
                                 ( TissueConfig::Instance()->GetTransitCellG1Duration()
                                    + TissueConfig::Instance()->GetSG2MDuration() );

            CellType cell_type;

            // Cells 0 1 2 and 3 are stem cells
            if (elem_index<4)
            {
                birth_time = - 2.0*(double)elem_index;
                cell_type = STEM;
            }
            else
            {
                cell_type = DIFFERENTIATED;
            }

            TissueCell cell(cell_type, HEALTHY, new StochasticDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> crypt(mesh, cells);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetSamplingTimestepMultiple(2);
        simulator.SetEndTime(100);
        simulator.SetOutputDirectory("TestVertexCryptLong");

        // Modified parameters to make cells equilibriate 
        TissueConfig::Instance()->SetAreaBasedDampingConstantParameter(1.0);
        TissueConfig::Instance()->SetDeformationEnergyParameter(10.0);
        TissueConfig::Instance()->SetMembraneSurfaceEnergyParameter(5.0);
        TissueConfig::Instance()->SetMaxTransitGenerations(2);
        

        // Make crypt shorter for sloughing 
        TissueConfig::Instance()->SetCryptLength(10.0);
                
        SloughingCellKiller<2> sloughing_cell_killer(&crypt);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }


    void TestMeshSurvivesSaveLoad() throw (Exception)
    {
        // Create mesh
        Cylindrical2dVertexMesh mesh(6, 12, 0.01, 2.0);

        // Create cells
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<mesh.GetNumElements(); elem_index++)
        {
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                ( TissueConfig::Instance()->GetTransitCellG1Duration()
                                    + TissueConfig::Instance()->GetSG2MDuration() );

            TissueCell cell(STEM, HEALTHY, new StochasticDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> crypt(mesh, cells);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("VertexCrypt2DMeshArchive");
        simulator.SetEndTime(0.1);

        // Memory leak (unconditional jump) without the following line.
        // The archiver assumes that a Solve has been called and simulation time has been set up properly.
        // In this test it hasn't so we need this to avoid memory leak.
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.1, 100);

        // Save
        TissueSimulationArchiver<2, VertexCryptSimulation2d>::Save(&simulator);

        // Load
        VertexCryptSimulation2d* p_simulator;
        p_simulator = TissueSimulationArchiver<2, VertexCryptSimulation2d>::Load("VertexCrypt2DMeshArchive", 0.0);

        // Create an identical mesh for comparison purposes
        Cylindrical2dVertexMesh mesh2(6, 12, 0.01, 2.0);

        // Compare meshes
        VertexMesh<2,2>& r_mesh = (static_cast<VertexBasedTissue<2>*>(&(p_simulator->rGetTissue())))->rGetMesh();
        CompareMeshes(&mesh2, &r_mesh);

        // Tidy up
        delete p_simulator;
    }

    /// \todo Add more tests (see #923)

};

#endif /*TESTVERTEXCRYPTSIMULATION2D_HPP_*/
