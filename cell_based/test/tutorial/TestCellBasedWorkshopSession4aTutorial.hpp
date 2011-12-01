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

#ifndef TESTCELLBASEDWORKSHOPSESSION4ATUTORIAL_HPP_
#define TESTCELLBASEDWORKSHOPSESSION4ATUTORIAL_HPP_

#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "RepulsionForce.hpp"
#include "NagaiHondaForce.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "SmartPointers.hpp"

class TestCellBasedWorkshopSession4aTutorial : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    void TestVertexBasedMonolayer() throw (Exception)
    {
        TS_FAIL("This test suite is currently too long for the continuous test pack.");
    	HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), TRANSIT);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Demo0");
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(20.0);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();
    }

    // Use a node-based cell population and demonstrate usage of cell killers
    void TestNodeBasedMonolayer() throw (Exception)
    {
        // Create mesh
        HoneycombMeshGenerator generator(2, 2); //**Changed**//
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh(); //**Changed**//
        NodesOnlyMesh<2> mesh; //**Changed**//
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh); //**Changed**//

        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), TRANSIT);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);//**Changed**//
        cell_population.SetMechanicsCutOffLength(1.5); //**Changed**// //to speed up simulations

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Demo1"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(40.0); //**Changed**//

        MAKE_PTR(RepulsionForce<2>, p_force); //**Changed**//
        simulator.AddForce(p_force);

        // Create cell killer and pass in to simulation
//        RandomCellKiller<2> cell_killer(&cell_population, 0.01);
//        simulator.AddCellKiller(&cell_killer);

        simulator.Solve();
    }

    // Use a mesh-based cell population
    void TestMeshBasedMonolayer() throw (Exception)
    {
        HoneycombMeshGenerator generator(2, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), TRANSIT);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells); //**Changed**//

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Demo2"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(40.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force); //**Changed**//
        simulator.AddForce(p_force);

        simulator.Solve();
   }

    // Add ghost nodes
    void TestMeshBasedMonolayerWithGhostNodes() throw (Exception)
    {
        HoneycombMeshGenerator generator(2, 2, 2); //**Changed**//
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();//**Changed**//
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), TRANSIT); //**Changed**//

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices); //**Changed**//

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Demo3"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(40.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();
    }

    // Make the simulation periodic
    void TestMeshBasedMonolayerPeriodic() throw (Exception)
    {
        CylindricalHoneycombMeshGenerator generator(5, 2, 2); //**Changed**//
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh(); //**Changed**//

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), TRANSIT);

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Demo4"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(40.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();
    }

    // Add boundary at x=0
    void TestMeshBasedMonolayerPeriodicSolidBottomBoundary() throw (Exception)
    {
        CylindricalHoneycombMeshGenerator generator(5, 2, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), TRANSIT);

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        OffLatticeSimulation<2> simulator(cell_population); //**Changed**//
        simulator.SetOutputDirectory("Demo5"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(40.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();
    }
};

#endif /*TESTCELLBASEDWORKSHOPSESSION4ATUTORIAL_HPP_*/
