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
#ifndef TESTTISSUESIMULATION_HPP_
#define TESTTISSUESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <ctime>
#include <cmath>

#include "TissueSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "StochasticWntCellCycleModelCellsGenerator.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "RandomCellKiller.hpp"
#include "AbstractCancerTestSuite.hpp"
#include "FixedDurationGenerationBasedCellCycleModelCellsGenerator.hpp"
#include "MeshBasedTissueWithGhostNodes.hpp"
#include "../../global/test/NumericFileComparison.hpp"
#include "CancerEventHandler.hpp"


// Simple subclass of TissueSimulation which just overloads StoppingEventHasOccurred
// for testing the stopping event functionality..
class TissueSimulationWithMyStoppingEvent : public TissueSimulation<2>
{
private:
    // define a stopping event with says stop if t>3.14
    bool StoppingEventHasOccurred()
    {
        return  (SimulationTime::Instance()->GetTime() > 3.1415);
    }

public:
    TissueSimulationWithMyStoppingEvent(AbstractTissue<2>& rTissue,
                                        std::vector<AbstractForce<2>* > forceCollection)
      : TissueSimulation<2>(rTissue, forceCollection)
   {
   }
};


/**
 *  Note: Most tests of TissueSimulation are in TestCryptSimulation2d
 */
class TestTissueSimulation : public AbstractCancerTestSuite
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

    void TestOutputStatistics() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        StochasticWntCellCycleModelCellsGenerator<2> cell_generator;
        cell_generator.GenerateForCrypt(cells,* p_mesh, location_indices, true);

        // Set up tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        TissueConfig::Instance()->SetOutputTissueAreas(true); // record the spheroid radius and apoptotic radius

        // Set up Wnt Gradient
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetTissue(tissue);

        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("TissueSimulationWritingProteins");
        simulator.SetEndTime(0.5);

        TS_ASSERT_DELTA(simulator.GetDt(), 1.0/120.0, 1e-12);

        TissueConfig::Instance()->SetOutputCellVariables(true);
        TissueConfig::Instance()->SetOutputCellCyclePhases(true);
        TissueConfig::Instance()->SetOutputCellAges(true);

        // Run tissue simulation
        TS_ASSERT_EQUALS(simulator.GetOutputDirectory(), "TissueSimulationWritingProteins");
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        OutputFileHandler handler("TissueSimulationWritingProteins", false);

        std::string cell_variables_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellvariables.dat";
        NumericFileComparison comp_cell_variables(cell_variables_file, "cancer/test/data/TissueSimulationWritingProteins/cellvariables.dat");
        TS_ASSERT(comp_cell_variables.CompareFiles(1e-2));

        std::string cell_cycle_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellcyclephases.dat";
        NumericFileComparison comp_cell_cycle(cell_cycle_file, "cancer/test/data/TissueSimulationWritingProteins/cellcyclephases.dat");
        TS_ASSERT(comp_cell_cycle.CompareFiles(1e-2));

        std::string cell_ages_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellages.dat";
        NumericFileComparison comp_cell_ages(cell_ages_file, "cancer/test/data/TissueSimulationWritingProteins/cellages.dat");
        TS_ASSERT(comp_cell_ages.CompareFiles(1e-2));

        // Tidy up
        WntConcentration<2>::Destroy();
    }


    /**
     * Test a tissue simulation with a cell killer.
     *
     * In this test, we solve a tissue simulation without ghost nodes and
     * check that the numbers of nodes and cells match at the end of the
     * simulation.
     */
    void TestTissueSimulationWithCellDeath() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each node. Give each cell a random birth time.
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                (TissueConfig::Instance()->GetStemCellG1Duration()
                                    + TissueConfig::Instance()->GetSG2MDuration() );
            TissueCell cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create a tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);

        // Create a force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("TestTissueSimulationWithCellDeath");
        simulator.SetEndTime(0.5);

        // Add cell killer
        RandomCellKiller<2> random_cell_killer(&tissue, 0.997877574);
        simulator.AddCellKiller(&random_cell_killer);

        // For coverage of an exception.
        simulator.SetUpdateTissueRule(false);
        TS_ASSERT_THROWS_THIS(simulator.Solve(),"Tissue has had births or deaths but mUpdateTissue is set to false, please set it to true.");
        CancerEventHandler::Reset(); // Otherwise logging has been started but not stopped due to exception above.

        simulator.SetUpdateTissueRule(true);
        simulator.Solve();

        // Check that the number of nodes is equal to the number of cells
        TS_ASSERT_EQUALS(simulator.rGetTissue().GetNumNodes(), simulator.rGetTissue().GetNumRealCells());

        // For coverage of these 'Get' functions
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 17u);
    }


    void TestTissueSimulationWithStoppingEvent() throw (Exception)
    {
        HoneycombMeshGenerator generator(2, 2, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells, one for each node. Give each cell a random birth time.
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells,* p_mesh, location_indices, true);

        // Create a tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);

        // Create a force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation WITH the stopping event
        TissueSimulationWithMyStoppingEvent simulator(tissue, force_collection);
        simulator.SetOutputDirectory("TestTissueSimWithStoppingEvent");

        // Set the end time to 10.0 - the stopping event is, however, t>3.1415.
        simulator.SetEndTime(10.0);

        // Run tissue simulation
        simulator.Solve();

        double time = SimulationTime::Instance()->GetTime();
        TS_ASSERT_DELTA(time, 3.1415, 1e-1); // big tol, doesn't matter, just want t~3.14 and t!=10
        // t should be strictly greater than the 3.1415
        TS_ASSERT_LESS_THAN(3.1415, time);
    }


    void TestApoptosisSpringLengths() throw (Exception)
    {
        unsigned num_cells_depth = 2;
        unsigned num_cells_width = 2;
        double crypt_length = num_cells_depth-0.0;
        double crypt_width = num_cells_width-0.0;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        TissueConfig* p_params = TissueConfig::Instance();
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);

        // Set up cells
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            TissueCell cell(TRANSIT, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*(p_params->GetTransitCellG1Duration()
                                               +p_params->GetSG2MDuration());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, location_indices);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from tissue and force law
        TissueSimulation<2> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("2dSpheroidApoptosis");
        simulator.SetEndTime(1.0);

        TissueConfig::Instance()->SetApoptosisTime(2.0);
        tissue.rGetCellUsingLocationIndex(14).StartApoptosis();
        tissue.rGetCellUsingLocationIndex(15).StartApoptosis();
        simulator.SetNoBirth(true);

        // Run tissue simulation
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

        TS_ASSERT_DELTA(a_b_separation, 0.5, 1e-1);
        TS_ASSERT_DELTA(a_c_separation, 0.75, 1e-1);
        TS_ASSERT_DELTA(c_d_separation, 1.0, 1e-1);

        // Reset end time and run tissue simulation
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

        TS_ASSERT_DELTA(a_b_separation, 0.01, 1e-1);
        TS_ASSERT_DELTA(a_c_separation, 0.5, 1e-1);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 4u);

        // Reset end time and run tissue simulation
        simulator.SetEndTime(2.01);
        simulator.Solve();

        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 2u);
    }

};
#endif /*TESTTISSUESIMULATION_HPP_*/
