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
#ifndef TESTTISSUESIMULATION3D_HPP_
#define TESTTISSUESIMULATION3D_HPP_

#include <cxxtest/TestSuite.h>

#include "TissueSimulation.hpp"
#include "FixedCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestTissueSimulation3d : public AbstractCancerTestSuite
{
private:

    MutableMesh<3,3> Make3dMesh(unsigned width=3, unsigned height=3, unsigned depth=3)
    {
        MutableMesh<3,3> mesh;
        mesh.ConstructCuboid(width,height,depth,true);
        TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);

        return mesh;
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
    void TestDoCellBirth() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CellType cell_type;
            unsigned generation;
            cell_type = STEM;
            generation = 0;
            TissueCell cell(cell_type, HEALTHY, new FixedCellCycleModel());
            cell.GetCellCycleModel()->SetGeneration(generation);
            cell.SetLocationIndex(i);
            if ( i == 50u)
            {
                cell.SetBirthTime(-50.0 );
            }

            cells.push_back(cell);
        }

        MeshBasedTissue<3> tissue(mesh,cells);
        TissueSimulation<3> simulator(tissue);

        unsigned num_births = simulator.DoCellBirth();

        TS_ASSERT_EQUALS(num_births, 1u);
    }

    void TestBirthOccursDuringSolve() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element");

        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<mesh.GetNumNodes()-1; i++)
        {
            cell.SetLocationIndex(i);
            cell.SetBirthTime(0.0);
            cells.push_back(cell);
        }

        // Setting last cell to undergo cell birth.
        cell.SetLocationIndex(mesh.GetNumNodes()-1);
        cell.SetBirthTime(-50.0);
        cells.push_back(cell);

        MeshBasedTissue<3> tissue(mesh,cells);
        TissueSimulation<3> simulator(tissue);

        TrianglesMeshWriter<3,3> mesh_writer1("Test3DCellBirth","StartMesh");
        mesh_writer1.WriteFilesUsingMesh(mesh);

        simulator.SetOutputDirectory("Test3DCellBirth");
        simulator.SetEndTime(1.0);

        simulator.Solve();

        TS_ASSERT_EQUALS(mesh.GetNumNodes(),5u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(),3u);

        TrianglesMeshWriter<3,3> mesh_writer2("Test3DCellBirth","EndMesh",false);
        mesh_writer2.WriteFilesUsingMesh(mesh);
    }

    void TestSolveMethodSpheroidSimulation3D() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<3,3> mesh_writer("TestSolveMethodSpheroidSimulation3DMesh","StartMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);

        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CellType cell_type;
            unsigned generation;
            cell_type = STEM;
            generation = 0;
            TissueCell cell(cell_type, HEALTHY, new FixedCellCycleModel());
            cell.GetCellCycleModel()->SetGeneration(generation);
            cell.SetLocationIndex(i);
            cell.SetBirthTime(-RandomNumberGenerator::Instance()->ranf()*
                               ( CancerParameters::Instance()->GetStemCellG1Duration()
                                 + CancerParameters::Instance()->GetSG2MDuration()   ));
            cells.push_back(cell);
        }

        MeshBasedTissue<3> tissue(mesh,cells);
        TissueSimulation<3> simulator(tissue);
        simulator.SetOutputDirectory("TestSolveMethodSpheroidSimulation3D");

        // Test SetSamplingTimestepMultiple method
        TS_ASSERT_EQUALS(simulator.mSamplingTimestepMultiple, 1u);
        simulator.SetSamplingTimestepMultiple(2);
        TS_ASSERT_EQUALS(simulator.mSamplingTimestepMultiple, 2u);

        simulator.SetEndTime(0.1);
        simulator.Solve();

        TrianglesMeshWriter<3,3> mesh_writer2("TestSolveMethodSpheroidSimulation3DMesh","EndMesh",false);
        mesh_writer2.WriteFilesUsingMesh(mesh);
    }


    void TestGhostNodesSpheroidSimulation3DandSave() throw (Exception)
    {
        unsigned width = 3;
        unsigned height = 3;
        unsigned depth = 3;

        MutableMesh<3,3> mesh = Make3dMesh(width,height,depth);
        TrianglesMeshWriter<3,3> mesh_writer("TestGhostNodesSpheroidSimulation3D","StartMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);

        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<TissueCell> cells;

        c_vector<double, 3> spheroid_centre;
        spheroid_centre[0] = 0.5*((double) width);
        spheroid_centre[1] = 0.5*((double) height);
        spheroid_centre[2] = 0.5*((double) depth);

        std::set<unsigned> ghost_node_indices;

        for (unsigned i=0; i<num_cells; i++)
        {
            CellType cell_type;
            unsigned generation;

            c_vector<double, 3> node_location = mesh.GetNode(i)->rGetLocation();

            unsigned min_spatial_dimension;
            if (width <= height && width <= depth)
            {
                min_spatial_dimension = width;
            }
            else
            {
                if (height <= depth)
                {
                    min_spatial_dimension = height;
                }
                else
                {
                    min_spatial_dimension = depth;
                }
            }
            if ( norm_2(node_location - spheroid_centre) > 0.5*sqrt(3)*1.01*((double) min_spatial_dimension)/3.0 )
            {
                ghost_node_indices.insert(i);
            }

            cell_type = STEM;
            generation = 0;
            TissueCell cell(cell_type, HEALTHY, new FixedCellCycleModel());
            cell.GetCellCycleModel()->SetGeneration(generation);
            cell.SetLocationIndex(i);
            cell.SetBirthTime(-RandomNumberGenerator::Instance()->ranf()*
                                (  CancerParameters::Instance()->GetStemCellG1Duration() +
                                   CancerParameters::Instance()->GetSG2MDuration()  ));
            cells.push_back(cell);
        }

        TS_ASSERT(ghost_node_indices.size() < num_cells);
        TS_ASSERT(ghost_node_indices.size() > 0);
        TS_ASSERT_EQUALS(ghost_node_indices.size(), 56u);

        // Test Save with a MeshBasedTissueWithGhostNodes
        MeshBasedTissueWithGhostNodes<3> tissue(mesh, cells, ghost_node_indices);
        TissueSimulation<3> simulator(tissue);
        simulator.SetOutputDirectory("TestGhostNodesSpheroidSimulation3D");
        simulator.SetEndTime(0.1);
        simulator.Solve();
        simulator.Save();

        // To generate results for below test
        // std::cout << mesh.GetNode(23u)->rGetLocation()[2] << std::endl << std::flush;

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        // Test Save with a MeshBasedTissue - one cell born during this.
        MeshBasedTissue<3> tissue2(mesh, cells);
        TissueSimulation<3> simulator2(tissue2);
        simulator2.SetOutputDirectory("TestGhostNodesSpheroidSimulation3DNoGhosts");
        simulator2.SetEndTime(0.1);
        simulator2.Solve();
        simulator2.Save();

        // To generate results for below test
        // std::cout << mesh.GetNode(23u)->rGetLocation()[2] << std::endl << std::flush;

    }

    void TestLoadOf3DSimulation() throw (Exception)
    {
        {   // With ghost nodes - 56 ghosts 8 real cells.
            TissueSimulation<3>* p_simulator = TissueSimulation<3>::Load("TestGhostNodesSpheroidSimulation3D", 0.1);
            unsigned num_cells = p_simulator->rGetTissue().GetNumRealCells();

            TS_ASSERT_EQUALS(num_cells, 8u);
            TS_ASSERT_DELTA(SimulationTime::Instance()->GetDimensionalisedTime(), 0.1, 1e-9);
            TS_ASSERT_DELTA(p_simulator->rGetTissue().GetLocationOfCell(p_simulator->rGetTissue().rGetCellUsingLocationIndex(23u))[2] , 0.911736, 1e-6);

            delete p_simulator;
        }

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        {   // Without ghost nodes - all 65 are real cells.
            TissueSimulation<3>* p_simulator = TissueSimulation<3>::Load("TestGhostNodesSpheroidSimulation3DNoGhosts", 0.1);
            unsigned num_cells = p_simulator->rGetTissue().GetNumRealCells();

            TS_ASSERT_EQUALS(num_cells, 65u);
            TS_ASSERT_DELTA(SimulationTime::Instance()->GetDimensionalisedTime(), 0.1, 1e-9);
            TS_ASSERT_DELTA(p_simulator->rGetTissue().GetLocationOfCell(p_simulator->rGetTissue().rGetCellUsingLocationIndex(23u))[2] , 1.13958, 1e-6);

            delete p_simulator;
        }
    }
};

#endif /*TESTTISSUESIMULATION3D_HPP_*/

