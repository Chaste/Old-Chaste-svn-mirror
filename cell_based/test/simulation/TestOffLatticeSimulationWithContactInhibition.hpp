/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTOFFLATTICESIMULATIONWITHCONTACTINHIBITION_HPP_
#define TESTOFFLATTICESIMULATIONWITHCONTACTINHIBITION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "ContactInhibitionOffLatticeSimulation.hpp"
#include "WildTypeCellMutationState.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SmartPointers.hpp"
#include "OutputFileHandler.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SimulationTime.hpp"
#include "CellLabel.hpp"
#include "OutputFileHandler.hpp"
#include "MutableMesh.hpp"
#include "Warnings.hpp"

class TestContactInhibitionOffLatticeSimulation : public CxxTest::TestSuite
{
public:

    void TestBoxContactInhibitionCellCycleModel()
      {
          // Set up SimulationTime
          SimulationTime* p_simulation_time = SimulationTime::Instance();
          p_simulation_time->SetStartTime(0.0);

          // Create a simple mesh
          HoneycombMeshGenerator generator(2, 2, 0);
          MutableMesh<2,2>* p_mesh = generator.GetMesh();

          // Create cell state
          MAKE_PTR(WildTypeCellMutationState, p_state);
          std::vector<CellPtr> cells;

          for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
          {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetCellProliferativeType(STEM);
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-10.0);
            p_cycle_model->SetQuiescentVolumeFraction(0.7);
            p_cycle_model->SetEquilibriumVolume(1.0);
            p_cycle_model->SetStemCellG1Duration(0.1);
            p_cycle_model->SetTransitCellG1Duration(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);
          }

          // Create a cell population
          MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

          // Create a force law
          MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
          p_force->SetCutOffLength(1.5);

          // Create a singleton class to store the volume of the cells and initialize it
          CellwiseData<2>* p_data = CellwiseData<2>::Instance();
          p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
          p_data->SetCellPopulation(&cell_population);

          for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
          {
              p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
          }

          // Create a contact inhibition simulator
          ContactInhibitionOffLatticeSimulation<2> simulator(cell_population);
          simulator.SetOutputDirectory("TestOffLatticeWithContactInhibition");
          simulator.SetEndTime(1.0);
          simulator.AddForce(p_force);

          // Run simulation
          TS_ASSERT_THROWS_NOTHING(simulator.Solve());

          // Test that the volumes of the cells are correct in CellWiseData
          cell_population.CreateVoronoiTessellation();

          for (MeshBasedCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                   cell_iter != cell_population.End();
                   ++cell_iter)
            {
                unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
                TS_ASSERT_DELTA(cell_population.GetVolumeOfVoronoiElement(node_index), p_data->GetValue(*cell_iter,0),1e-4);
            }

          // Tidy up
          SimulationTime::Destroy();
          RandomNumberGenerator::Destroy();
          CellwiseData<2>::Destroy();
      }

    void TestArchiving() throw (Exception)
        {
			EXIT_IF_PARALLEL;

			// Set up SimulationTime
			SimulationTime* p_simulation_time = SimulationTime::Instance();
			p_simulation_time->SetStartTime(0.0);

			// Create a simple mesh
			HoneycombMeshGenerator generator(2, 2, 0);
			MutableMesh<2,2>* p_mesh = generator.GetMesh();

			// Create cell state
			MAKE_PTR(WildTypeCellMutationState, p_state);
			std::vector<CellPtr> cells;

			for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
			{
				ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
				p_cycle_model->SetCellProliferativeType(STEM);
				p_cycle_model->SetDimension(2);
				p_cycle_model->SetBirthTime(-1.0);
				p_cycle_model->SetQuiescentVolumeFraction(0.7);
				p_cycle_model->SetEquilibriumVolume(1.0);
				p_cycle_model->SetStemCellG1Duration(0.1);
				p_cycle_model->SetTransitCellG1Duration(0.1);

				CellPtr p_cell(new Cell(p_state, p_cycle_model));
				p_cell->InitialiseCellCycleModel();
				cells.push_back(p_cell);
			}

			// Create a cell population
			MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

			// Create a force law
			MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
			p_force->SetCutOffLength(1.5);

			// Create a singleton class to store the volume of the cells and initialize it
			CellwiseData<2>* p_data = CellwiseData<2>::Instance();
			p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
			p_data->SetCellPopulation(&cell_population);

			for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
			{
				p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
			}

			// Create a contact inhibition simulator
			ContactInhibitionOffLatticeSimulation<2> simulator(cell_population);
			simulator.SetOutputDirectory("TestContactInhibitionOffLatticeSimulationSaveAndLoad");
			double end_time=0.01;
			simulator.SetEndTime(end_time);
			simulator.AddForce(p_force);

			// Run simulation
			simulator.Solve();

            CellBasedSimulationArchiver<2, ContactInhibitionOffLatticeSimulation<2> >::Save(&simulator);

            TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 4u);
            TS_ASSERT_EQUALS((static_cast<MeshBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation())))->GetNumRealCells(), 4u);

            TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.01, 1e-9);
            CellPtr p_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(3);
            TS_ASSERT_DELTA(p_cell->GetAge(), 1.01, 1e-4);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Load simulation
            ContactInhibitionOffLatticeSimulation<2>* p_simulator
                = CellBasedSimulationArchiver<2, ContactInhibitionOffLatticeSimulation<2> >::Load("TestContactInhibitionOffLatticeSimulationSaveAndLoad", end_time);

            p_simulator->SetEndTime(0.2);

            TS_ASSERT_EQUALS(p_simulator->rGetCellPopulation().GetNumRealCells(), 4u);
            TS_ASSERT_EQUALS((static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation())))->GetNumRealCells(), 4u);

            TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.01, 1e-9);
            CellPtr p_cell2 = p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(3);
            TS_ASSERT_DELTA(p_cell2->GetAge(), 1.01, 1e-4);

            // Run simulation
            TS_ASSERT_THROWS_NOTHING(p_simulator->Solve());

            // Tidy up
            delete p_simulator;

            // Test Warnings
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
            Warnings::QuietDestroy();
        }
};

#endif /*TESTOFFLATTICESIMULATIONWITHCONTACTINHIBITION_HPP_*/
