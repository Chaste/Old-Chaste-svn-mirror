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
#ifndef TESTTISSUESIMULATION_HPP_
#define TESTTISSUESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <stdio.h>
#include <ctime>
#include <math.h>

#include "TissueSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "FixedCellCycleModel.hpp"
#include "CryptProjectionSpringSystem.hpp"
#include "RandomCellKiller.hpp"
#include "RadialSloughingCellKiller.hpp"
#include "AbstractCancerTestSuite.hpp"


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
        RandomNumberGenerator *p_gen=RandomNumberGenerator::Instance();
        CancerParameters::Instance()->SetHepaOneParameters();

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new TysonNovakCellCycleModel());
            double birth_time = -1.0*p_gen->ranf();
            cell.SetLocationIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Set up tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        tissue.SetWriteTissueAreas(true); // record the spheroid radius and necrotic radius

        Meineke2001SpringSystem<2> spring_system(tissue);
        spring_system.UseCutoffPoint(1.5);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, &spring_system);
        simulator.SetOutputDirectory("TissueSimulationWritingProteins");
        simulator.SetEndTime(0.5);
        simulator.SetOutputCellVariables(true);
        simulator.SetOutputCellCyclePhases(true);

        // Run tissue simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        OutputFileHandler handler("TissueSimulationWritingProteins",false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellvariables.dat";
        TS_ASSERT_EQUALS(system(("diff " + results_file + " cancer/test/data/TissueSimulationWritingProteins/cellvariables.dat").c_str()), 0);
    }

    /**
     *  Test a tissue simulation with a non-Meineke spring system.
     *
     *  This test consists of a standard crypt projection model simulation with a
     *  radial sloughing cell killer, a crypt projection cell cycle model that
     *  depends on a radial Wnt gradient, and the crypt projection model spring
     *  system, and store the results for use in later archiving tests.
     */
    void TestTissueSimulationWithCryptProjectionSpringSystem() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->SetWntStemThreshold(0.95);

        double a = 0.2;
        double b = 2.0;
        p_params->SetCryptProjectionParameterA(a);
        p_params->SetCryptProjectionParameterB(b);

        // Set up mesh
        int num_cells_depth = 20;
        int num_cells_width = 20;
        unsigned thickness_of_ghost_layer = 3;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        c_vector<double,2> width_extremes = p_mesh->GetWidthExtremes(0u);
        c_vector<double,2> height_extremes = p_mesh->GetWidthExtremes(1u);

        double width_of_mesh = (num_cells_width/(num_cells_width+2.0*thickness_of_ghost_layer))*(width_extremes[1] - width_extremes[0]);
        double height_of_mesh = (num_cells_depth/(num_cells_depth+2.0*thickness_of_ghost_layer))*(height_extremes[1] - height_extremes[0]);

        p_mesh->Translate(-width_of_mesh/2,-height_of_mesh/2);

        // To start off with, set up all cells to be of type TRANSIT
        std::vector<TissueCell> cells;

        std::cout << "num nodes = " << p_mesh->GetNumNodes() << "\n" << std::flush;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(TRANSIT, HEALTHY, new SimpleWntCellCycleModel());
            cell.InitialiseCellCycleModel();
            double birth_time = - RandomNumberGenerator::Instance()->ranf()*
                                  ( p_params->GetTransitCellG1Duration()
                                   +p_params->GetSG2MDuration());
            cell.SetLocationIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Make a tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);

        // Set up the Wnt gradient
        WntConcentration::Instance()->SetType(RADIAL);
        WntConcentration::Instance()->SetTissue(crypt);

        // Create the spring system
        CryptProjectionSpringSystem spring_system(crypt);

        // Make a tissue simulation
        TissueSimulation<2> crypt_projection_simulator(crypt, &spring_system, false, false);

        // Create a radial cell killer and pass it in to the tissue simulation
        c_vector<double,2> centre = zero_vector<double>(2);
        double crypt_radius = pow(CancerParameters::Instance()->GetCryptLength()/a, 1.0/b);

        RadialSloughingCellKiller killer(&crypt, centre, crypt_radius);
        crypt_projection_simulator.AddCellKiller(&killer);

        // Set up the simulation
        crypt_projection_simulator.SetOutputDirectory("CryptProjectionSimulation");
        crypt_projection_simulator.SetEndTime(0.25);

        // Run the simulation
        TS_ASSERT_THROWS_NOTHING(crypt_projection_simulator.Solve());

        // These cells just divided and have been gradually moving apart.
        // These results are from time 0.25.
        std::vector<double> node_302_location = crypt_projection_simulator.GetNodeLocation(302);
        std::vector<double> node_506_location = crypt_projection_simulator.GetNodeLocation(506);
        c_vector<double, 2> distance_between;
        distance_between(0) = node_506_location[0]-node_302_location[0];
        distance_between(1) = node_506_location[1]-node_302_location[1];
        TS_ASSERT_DELTA(norm_2(distance_between), 0.7029, 1e-3);

        // Test the Wnt gradient result
        TissueCell* p_cell = &(crypt.rGetCellAtNodeIndex(302));
        TS_ASSERT_DELTA(WntConcentration::Instance()->GetWntLevel(p_cell), 0.999, 1e-3);
        p_cell = &(crypt.rGetCellAtNodeIndex(506));
        TS_ASSERT_DELTA(WntConcentration::Instance()->GetWntLevel(p_cell), 0.989, 1e-3);

        // Tidy up
        WntConcentration::Destroy();
    }

    /**
     *  Test a tissue simulation with a cell killer.
     *
     *  In this test, we solve a tissue simulation without ghost nodes and
     *  check that the numbers of nodes and cells match at the end of the
     *  simulation.
     *
     *  This test currently fails so is commented out.
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
                                (CancerParameters::Instance()->GetStemCellG1Duration()
                                    + CancerParameters::Instance()->GetSG2MDuration() );
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            cell.SetLocationIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create a tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);

        // Create a mechanics system
        Meineke2001SpringSystem<2> spring_system(tissue);
        spring_system.UseCutoffPoint(1.5);

        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, &spring_system);
        simulator.SetOutputDirectory("TestTissueSimulationWithCellDeath");
        simulator.SetEndTime(0.5);

        // Add cell killer
        RandomCellKiller<2> random_cell_killer(&tissue, 0.05);
        simulator.AddCellKiller(&random_cell_killer);

        // Run tissue simulation
        simulator.Solve();

        // Check that the number of nodes is equal to the number of cells
        TS_ASSERT_EQUALS(simulator.rGetTissue().GetNumNodes(), simulator.rGetTissue().GetNumRealCells());

        // For coverage of these 'Get' functions
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 17u);
    }

};
#endif /*TESTTISSUESIMULATION_HPP_*/
