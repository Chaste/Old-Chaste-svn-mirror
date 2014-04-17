/*

Copyright (c) 2005-2014, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTFARHADIFARTYPEMODIFIER_HPP_
#define TESTFARHADIFARTYPEMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include <ctime>
#include "SmartPointers.hpp"

#include "FarhadifarTypeModifier.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "NagaiHondaForce.hpp"
#include "CellBasedEventHandler.hpp"

// This test is only run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestFarhadifarTypeModifier : public AbstractCellBasedTestSuite
{
public:
    void TestSetupSolveException() throw (Exception)
    {
        // First set up SimulationTime (this is usually handled by a simulation object)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a FarhadifarTypeModifier
        MAKE_PTR(FarhadifarTypeModifier<2>, p_modifier);

        // Create a cell population whose type should not be used with a FarhadifarTypeModifier
        HoneycombMeshGenerator generator(4, 4, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> population(*p_mesh, cells);

        // Test that the correct exception is thrown if we try to call UpdateTargetAreas() on the population
        TS_ASSERT_THROWS_THIS(p_modifier->SetupSolve(population, "unused_argument"),
                "AbstractTargetAreaModifiers are to be used with a VertexBasedCellPopulation only");
        CellBasedEventHandler::Reset(); // Otherwise logging has been started but not stopped due to exception above.

    }

   void TestFarhadifarTypeModifierMethods() throw (Exception)
    {

        // create our modifier
        MAKE_PTR(FarhadifarTypeModifier<2>,p_growth_modifier);

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3*0.25, 3);

        // Create mesh
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        //MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

//        cells[0]->SetCellProliferativeType(p_diff_type);
//        cells[4]->SetCellProliferativeType(p_diff_type);

        // the cells have varying ages
        for (unsigned i=0; i<cells.size(); i++)
        {
            double birth_time;
            if (i > 3)
            {
                birth_time = 0 - (20. + 0.8*(i-4.));
            }
            else
            {
                birth_time = 0.0 - 2*i;
            }
            cells[i]->SetBirthTime(birth_time);
        }

        // Create a cell population

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.InitialiseCells(); // this method must be called explicitly as there is no simulation

        // Check UpdateTargetAreaOfCell() method
        for (VertexMesh<2,2>::VertexElementIterator iter = p_mesh->GetElementIteratorBegin();
                iter != p_mesh->GetElementIteratorEnd();
                ++iter)
        {
            unsigned elem_index = iter->GetIndex();
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(elem_index);

            // get the mature cell target area
            double expected_area = p_growth_modifier->GetReferenceTargetArea();

            // have the growth modifier update the cell target area
            p_growth_modifier->UpdateTargetAreaOfCell(p_cell);

            if (elem_index>3)
            {
                expected_area *= (1 + (0.8*(elem_index-4.))/4.0);
            }

            double actual_area = p_cell->GetCellData()->GetItem("target area");

            TS_ASSERT_DELTA(actual_area, expected_area, 1e-12);
        }

        // check that the growth rule for individual cells works together
        // with the UpdateTargetAreas method

        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);
        CellPtr p_cell_1 = cell_population.GetCellUsingLocationIndex(1);
        CellPtr p_cell_3 = cell_population.GetCellUsingLocationIndex(3);

        // Make cell 1 and 3 undergo apoptosis
        p_cell_1->StartApoptosis();
        p_cell_3->StartApoptosis();

        // update the target areas
        p_growth_modifier->UpdateTargetAreas(cell_population);

        double actual_area_0 = p_cell_0->GetCellData()->GetItem("target area");
        double actual_area_1 = p_cell_1->GetCellData()->GetItem("target area");
        double actual_area_3 = p_cell_3->GetCellData()->GetItem("target area");

        double expected_area_0 = 1.0;
        double expected_area_1 = p_growth_modifier->GetReferenceTargetArea();
        double expected_area_3 = p_growth_modifier->GetReferenceTargetArea();

        TS_ASSERT_DELTA(actual_area_0, expected_area_0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1, expected_area_1, 1e-12);
        TS_ASSERT_DELTA(actual_area_3, expected_area_3, 1e-12);

        // Run for one time step
        p_simulation_time->IncrementTimeOneStep();
        p_growth_modifier->UpdateTargetAreas(cell_population);

        double actual_area_0_after_dt = p_cell_0->GetCellData()->GetItem("target area");
        double actual_area_1_after_dt = p_cell_1->GetCellData()->GetItem("target area");
        double actual_area_3_after_dt = p_cell_3->GetCellData()->GetItem("target area");

        // The target areas of cells 1 and 3 should have halved
        expected_area_0 = ( p_growth_modifier->GetReferenceTargetArea() );

        TS_ASSERT_DELTA(actual_area_0_after_dt, expected_area_0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1_after_dt, 0.5*expected_area_1, 1e-12);
        TS_ASSERT_DELTA(actual_area_3_after_dt, 0.5*expected_area_3, 1e-12);

        // Make cell 0 undergo apoptosis
        p_cell_0->StartApoptosis();

        // Now run on for a further time step
        p_simulation_time->IncrementTimeOneStep();
        p_growth_modifier->UpdateTargetAreas(cell_population);

        double actual_area_0_after_2dt = p_cell_0->GetCellData()->GetItem("target area");
        double actual_area_1_after_2dt = p_cell_1->GetCellData()->GetItem("target area");
        double actual_area_3_after_2dt = p_cell_3->GetCellData()->GetItem("target area");

        // Cells 1 and 3 should now have zero target area and the target area of cell 0 should have halved
        TS_ASSERT_DELTA(actual_area_0_after_2dt, 0.5*expected_area_0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1_after_2dt, 0.0, 1e-12);
        TS_ASSERT_DELTA(actual_area_3_after_2dt, 0.0, 1e-12);

        // Now run on for even further, for coverage
        p_simulation_time->IncrementTimeOneStep();
        p_growth_modifier->UpdateTargetAreas(cell_population);

        double actual_area_0_after_3dt = p_cell_0->GetCellData()->GetItem("target area");
        double actual_area_1_after_3dt = p_cell_1->GetCellData()->GetItem("target area");
        double actual_area_3_after_3dt = p_cell_3->GetCellData()->GetItem("target area");

        // All apoptotic cells should now have zero target area
        TS_ASSERT_DELTA(actual_area_0_after_3dt, 0.0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1_after_3dt, 0.0, 1e-12);
        TS_ASSERT_DELTA(actual_area_3_after_3dt, 0.0, 1e-12);
    }

    void TestTargetAreaOfDaughterCells()
    {
        // Plan: initialise a cell and have it divide at a fixed time. If cell divides at t, check that:
        // target area of mother cell at t - dt = mature target area
        // target area of daughter cells at t = half of that
        // target area of daughter cells at t + dt = half of mature target area + little increment

        // Create a simple 2D MutableVertexMesh with only one cell
        HoneycombVertexMeshGenerator generator(1, 1);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cell.
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();

        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_transit_type);
        double birth_time = -11.0; // The cell cycle duration is 12
        p_cell->SetBirthTime(birth_time);
        cells.push_back(p_cell);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestTargetAreaOfDaughterCells");
        simulator.SetEndTime(0.997);
        //dt=0.002

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // The FarhadifarTypeModifier
        MAKE_PTR(FarhadifarTypeModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        // Run simulation

        simulator.Solve();

        unsigned num_cells_before_division = simulator.rGetCellPopulation().GetNumRealCells();
        // We should only have one cell now
        TS_ASSERT_EQUALS(num_cells_before_division,1u);

        // This is the cell from before, let's see what its target area is
        double target_area_before_division = p_cell->GetCellData()->GetItem("target area");
        double expected_area = p_growth_modifier->GetReferenceTargetArea() +
                                  (p_cell->GetAge() - 8.0)/(p_cell->GetCellCycleModel()->GetG2Duration());
        TS_ASSERT_DELTA(target_area_before_division,expected_area,1e-9);

        // Now we adjust the end time and run the simulation a bit further
        simulator.SetEndTime(1.001);
        simulator.Solve();

        unsigned num_cells_at_division = simulator.rGetCellPopulation().GetNumRealCells();
        // We should have two cells now
        TS_ASSERT_EQUALS(num_cells_at_division,2u);

        // Iterate over the cells, checking their target areas
        for (VertexBasedCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            double target_area_at_division = cell_iter->GetCellData()->GetItem("target area");
            TS_ASSERT_DELTA(target_area_at_division,1.0,1e-9);
        }

        // Now we do the same thing again
        simulator.SetEndTime(1.003);
        simulator.Solve();

        unsigned num_cells_after_division = simulator.rGetCellPopulation().GetNumRealCells();

        // We should still have two cells
        TS_ASSERT_EQUALS(num_cells_after_division,2u);

        for (VertexBasedCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            double target_area_after_division = cell_iter->GetCellData()->GetItem("target area");

            TS_ASSERT_DELTA(target_area_after_division,1.0,1e-9);
        }
    }


    void TestFarhadifarTypeModifierArchiving()
    {
        //create a file for archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "growth_modifier.arch";

        // separate scope to write the archive
        {
            // initialise a growth modifier and set a non-standard mature target area
            MAKE_PTR(FarhadifarTypeModifier<2>, p_growth_modifier);
            p_growth_modifier->SetReferenceTargetArea(14.3);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the growth modifier to the archive
            boost::shared_ptr<AbstractCellBasedSimulationModifier<2,2> > p_abstract_simulation_modifier = p_growth_modifier;
            output_arch << p_abstract_simulation_modifier;
        }

        // separate scope to read the archive
        {
            //MAKE_PTR((AbstractCellBasedSimulationModifier<2,2>) , p_growth_modifier);
            boost::shared_ptr<AbstractCellBasedSimulationModifier<2,2> > p_growth_modifier(new AbstractCellBasedSimulationModifier<2,2>);
            // Restore the modifier
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_growth_modifier;

            // get a pointer of type growth modifier
            boost::shared_ptr<FarhadifarTypeModifier<2> > p_real_growth_modifier =
                    boost::dynamic_pointer_cast<FarhadifarTypeModifier<2> >(p_growth_modifier);
            TS_ASSERT(p_real_growth_modifier != NULL);

            // see whether we read out the correct target area
            double mature_target_area = p_real_growth_modifier->GetReferenceTargetArea();
            TS_ASSERT_DELTA(mature_target_area, 14.3, 1e-9);
        }
    }
};

#endif /*TESTFARHADIFARTYPEMODIFIER_HPP_*/
