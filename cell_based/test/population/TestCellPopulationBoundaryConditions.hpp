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

#ifndef TESTCELLPOPULATIONBOUNDARYCONDITIONS_HPP_
#define TESTCELLPOPULATIONBOUNDARYCONDITIONS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "SphereGeometryBoundaryCondition.hpp"
#include "CellwiseData.hpp"
#include "TrianglesMeshReader.hpp"
#include "WildTypeCellMutationState.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "ArchiveOpener.hpp"
#include "ArchiveLocationInfo.hpp"

/**
 * This class contains tests for methods on classes inheriting from AbstractCellPopulationBoundaryCondition.
 */
class TestCellPopulationBoundaryConditions : public AbstractCellBasedTestSuite
{
public:

    void TestPlaneBoundaryConditionWithNodeBasedCellPopulation() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell population boundary condition
        c_vector<double,2> point = zero_vector<double>(2);
        point(0) = 2.0;
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        PlaneBoundaryCondition<2> boundary_condition(&cell_population, point, normal);

        TS_ASSERT_EQUALS(boundary_condition.GetIdentifier(), "PlaneBoundaryCondition-2");

        // Impose boundary condition
        std::vector<c_vector<double,2> > old_locations;
        old_locations.reserve(cell_population.GetNumNodes());
        for (std::list<CellPtr>::iterator cell_iter = cell_population.rGetCells().begin();
             cell_iter != cell_population.rGetCells().end();
             ++cell_iter)
        {
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            old_locations.push_back(p_node->rGetLocation());
        }

        boundary_condition.ImposeBoundaryCondition(old_locations);

        // Test that all nodes satisfy the boundary condition
        for (std::list<CellPtr>::iterator cell_iter = cell_population.rGetCells().begin();
             cell_iter != cell_population.rGetCells().end();
             ++cell_iter)
        {
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location = p_node->rGetLocation();
            if (old_locations[p_node->GetIndex()][0] < 2.0)
            {
                TS_ASSERT_DELTA(2.0, location[0], 1e-6);
                TS_ASSERT_DELTA(location[1], old_locations[p_node->GetIndex()][1], 1e-6);
            }
            else
            {
                TS_ASSERT_DELTA(location[0], old_locations[p_node->GetIndex()][0], 1e-6);
                TS_ASSERT_DELTA(location[1], old_locations[p_node->GetIndex()][1], 1e-6);
            }
        }

        // Test VerifyBoundaryCondition() method
        TS_ASSERT_EQUALS(boundary_condition.VerifyBoundaryCondition(), true);

        // For coverage, test VerifyBoundaryCondition() method in the case DIM != 2
        PlaneBoundaryCondition<3> plane_boundary_condition_3d(NULL, zero_vector<double>(3), unit_vector<double>(3,2));
        TS_ASSERT_THROWS_THIS(plane_boundary_condition_3d.VerifyBoundaryCondition(),
                              "PlaneBoundaryCondition is not yet implemented in 1D or 3D");
    }

    void TestPlaneBoundaryConditionWithVertexBasedCellPopulation() throw(Exception)
    {
        // Create mesh
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);


        // Set up cell population boundary condition x>1
        c_vector<double,2> point = zero_vector<double>(2);
        point(0) = 1.0;
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        PlaneBoundaryCondition<2> boundary_condition(&cell_population, point, normal);

        TS_ASSERT_EQUALS(boundary_condition.GetIdentifier(), "PlaneBoundaryCondition-2");

        // Impose boundary condition
        std::vector<c_vector<double,2> > old_locations;
        old_locations.reserve(cell_population.GetNumNodes());

        for (MutableVertexMesh<2,2>::NodeIterator node_iter = cell_population.rGetMesh().GetNodeIteratorBegin();
                node_iter != cell_population.rGetMesh().GetNodeIteratorEnd();
                ++node_iter)
        {
            old_locations.push_back(node_iter->rGetLocation());
        }

        boundary_condition.ImposeBoundaryCondition(old_locations);

        // Test that all nodes satisfy the boundary condition
        for (MutableVertexMesh<2,2>::NodeIterator node_iter = cell_population.rGetMesh().GetNodeIteratorBegin();
                node_iter != cell_population.rGetMesh().GetNodeIteratorEnd();
                ++node_iter)
        {
            c_vector<double, 2> location = node_iter->rGetLocation();
            if (old_locations[node_iter->GetIndex()][0] < 1.0)
            {
                TS_ASSERT_DELTA(1.0, location[0], 1e-6);
                TS_ASSERT_DELTA(location[1], old_locations[node_iter->GetIndex()][1], 1e-6);
            }
            else
            {
                TS_ASSERT_DELTA(location[0], old_locations[node_iter->GetIndex()][0], 1e-6);
                TS_ASSERT_DELTA(location[1], old_locations[node_iter->GetIndex()][1], 1e-6);
            }
        }

        // Test VerifyBoundaryCondition() method
        TS_ASSERT_EQUALS(boundary_condition.VerifyBoundaryCondition(), true);

        // For coverage, test VerifyBoundaryCondition() method in the case DIM != 2
        PlaneBoundaryCondition<3> plane_boundary_condition_3d(NULL, zero_vector<double>(3), unit_vector<double>(3,2));
        TS_ASSERT_THROWS_THIS(plane_boundary_condition_3d.VerifyBoundaryCondition(),
                              "PlaneBoundaryCondition is not yet implemented in 1D or 3D");
    }

//    void TestPlaneBoundaryConditionExceptions() throw(Exception)
//    {
//        // Create a simple 2D PottsMesh
//        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
//        PottsMesh<2>* p_mesh = generator.GetMesh();
//
//        // Create cells
//        std::vector<CellPtr> cells;
//        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), DIFFERENTIATED);
//
//        // Create cell population
//        PottsBasedCellPopulation<2> potts_cell_population(*p_mesh, cells);
//
//        // Attempt to set up cell population boundary condition
//        c_vector<double,2> point = zero_vector<double>(2);
//        c_vector<double,2> normal = zero_vector<double>(2);
//        normal(0) = 1.0;
//        TS_ASSERT_THROWS_THIS(PlaneBoundaryCondition<2> plane_boundary_condition(&potts_cell_population, point, normal),
//            "PlaneBoundaryCondition require a subclass of AbstractOffLatticeCellPopulation.");
//    }

    void TestSphereGeometryBoundaryCondition() throw (Exception)
    {
        // We first test that the correct exception is thrown in 1D
        TrianglesMeshReader<1,1> mesh_reader_1d("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> generating_mesh_1d;
        generating_mesh_1d.ConstructFromMeshReader(mesh_reader_1d);

        NodesOnlyMesh<1> mesh_1d;
        mesh_1d.ConstructNodesWithoutMesh(generating_mesh_1d);

        std::vector<CellPtr> cells_1d;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 1> cells_generator_1d;
        cells_generator_1d.GenerateBasic(cells_1d, mesh_1d.GetNumNodes());

        NodeBasedCellPopulation<1> population_1d(mesh_1d, cells_1d);

        c_vector<double,1> centre_1d = zero_vector<double>(1);

        TS_ASSERT_THROWS_THIS(SphereGeometryBoundaryCondition<1> bc_1d(&population_1d, centre_1d, 1.0),
            "This boundary condition is not implemented in 1D.");

        // Next we test that the correct exception is thrown if not using a NodeBasedCellPopulation
        TrianglesMeshReader<2,2> mesh_reader_2d("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh_2d;
        mesh_2d.ConstructFromMeshReader(mesh_reader_2d);

        std::vector<CellPtr> cells_2d;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator_2d;
        cells_generator_2d.GenerateBasic(cells_2d, mesh_2d.GetNumNodes());

        MeshBasedCellPopulation<2> population_2d(mesh_2d, cells_2d);

        c_vector<double,2> centre_2d = zero_vector<double>(2);

        TS_ASSERT_THROWS_THIS(SphereGeometryBoundaryCondition<2> bc_2d(&population_2d, centre_2d, 1.0),
            "A NodeBasedCellPopulation must be used with this boundary condition object.");

        // We now test the methods of this class
        TrianglesMeshReader<3,3> mesh_reader_3d("mesh/test/data/cube_136_elements");
        TetrahedralMesh<3,3> generating_mesh_3d;
        generating_mesh_3d.ConstructFromMeshReader(mesh_reader_3d);

        NodesOnlyMesh<3> mesh_3d;
        mesh_3d.ConstructNodesWithoutMesh(generating_mesh_3d);

        std::vector<CellPtr> cells_3d;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator_3d;
        cells_generator_3d.GenerateBasic(cells_3d, mesh_3d.GetNumNodes());

        NodeBasedCellPopulation<3> population_3d(mesh_3d, cells_3d);

        c_vector<double,3> centre_3d = zero_vector<double>(3);

        SphereGeometryBoundaryCondition<3> bc_3d(&population_3d, centre_3d, 0.4, 1e-4);

        // Test that member variables were initialised correctly
        TS_ASSERT_DELTA(bc_3d.rGetCentreOfSphere()[0], centre_3d(0), 1e-4);
        TS_ASSERT_DELTA(bc_3d.rGetCentreOfSphere()[1], centre_3d(1), 1e-4);
        TS_ASSERT_DELTA(bc_3d.GetRadiusOfSphere(), 0.4, 1e-4);

        TS_ASSERT_EQUALS(bc_3d.VerifyBoundaryCondition(), false);

        // Store the location of each node prior to imposing the boundary condition
        std::vector<c_vector<double,3> > old_locations;
        old_locations.reserve(population_3d.GetNumNodes());
        for (std::list<CellPtr>::iterator cell_iter = population_3d.rGetCells().begin();
             cell_iter != population_3d.rGetCells().end();
             ++cell_iter)
        {
            c_vector<double,3> location = population_3d.GetLocationOfCellCentre(*cell_iter);
            old_locations.push_back(location);
        }

        bc_3d.ImposeBoundaryCondition(old_locations);

        // Test that the boundary condition was imposed correctly
        TS_ASSERT_EQUALS(bc_3d.VerifyBoundaryCondition(), true);
    }

    void TestArchivingOfPlaneBoundaryCondition() throw (Exception)
    {
        // Set up singleton classes
        OutputFileHandler handler("archive", false); // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "PlaneBoundaryCondition.arch";

        {
            // Create an output archive
            PlaneBoundaryCondition<2> boundary_condition(NULL, zero_vector<double>(2), unit_vector<double>(2,1));

            TS_ASSERT_DELTA(boundary_condition.rGetPointOnPlane()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(boundary_condition.rGetPointOnPlane()[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(boundary_condition.rGetNormalToPlane()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(boundary_condition.rGetNormalToPlane()[1], 1.0, 1e-6);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            AbstractCellPopulationBoundaryCondition<2>* const p_boundary_condition = &boundary_condition;
            output_arch << p_boundary_condition;
        }

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCellPopulationBoundaryCondition<2>* p_boundary_condition;

            // Restore from the archive
            input_arch >> p_boundary_condition;

            // Test we have restored the plane geometry correctly
            TS_ASSERT_DELTA(static_cast<PlaneBoundaryCondition<2>*>(p_boundary_condition)->rGetPointOnPlane()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(static_cast<PlaneBoundaryCondition<2>*>(p_boundary_condition)->rGetPointOnPlane()[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(static_cast<PlaneBoundaryCondition<2>*>(p_boundary_condition)->rGetNormalToPlane()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(static_cast<PlaneBoundaryCondition<2>*>(p_boundary_condition)->rGetNormalToPlane()[1], 1.0, 1e-6);

            // Tidy up
            delete p_boundary_condition;
       }
    }

    void TestArchivingOfSphereGeometryBoundaryCondition() throw (Exception)
    {
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,
1);

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);
        ArchiveLocationInfo::SetMeshFilename("mesh");

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> population(mesh, cells);

        c_vector<double,2> centre = zero_vector<double>(2);
        centre(0) = 0.5;
        centre(1) = 0.7;

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "SphereGeometryBoundaryCondition.arch";
        ArchiveLocationInfo::SetMeshFilename("SphereGeometryBoundaryCondition");

        {
            SphereGeometryBoundaryCondition<2> bc(&population, centre, 0.56, 1e-3);

            // Create an output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Serialize via pointer
            AbstractCellPopulationBoundaryCondition<2>* const p_bc = &bc;
            (*p_arch) << p_bc;
        }

        {
            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCellPopulationBoundaryCondition<2>* p_bc;

            // Restore from the archive
            (*p_arch) >> p_bc;

            // Test we have restored the object correctly
            TS_ASSERT_DELTA(static_cast<SphereGeometryBoundaryCondition<2>*>(p_bc)->rGetCentreOfSphere()[0], 0.5, 1e-6);
            TS_ASSERT_DELTA(static_cast<SphereGeometryBoundaryCondition<2>*>(p_bc)->rGetCentreOfSphere()[1], 0.7, 1e-6);
            TS_ASSERT_DELTA(static_cast<SphereGeometryBoundaryCondition<2>*>(p_bc)->GetRadiusOfSphere(), 0.56, 1e-6);

            // Tidy up
            delete p_bc->mpCellPopulation;
            delete p_bc;
       }
    }

    void TestCellBoundaryConditionsOutputParameters()
    {
        std::string output_directory = "TestCellBoundaryConditionsOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with PlaneBoundaryCondition
        PlaneBoundaryCondition<2> plane_boundary_condition(NULL, zero_vector<double>(2), unit_vector<double>(2,1));
        TS_ASSERT_EQUALS(plane_boundary_condition.GetIdentifier(), "PlaneBoundaryCondition-2");

        out_stream plane_boundary_condition_parameter_file = output_file_handler.OpenOutputFile("plane_results.parameters");
        plane_boundary_condition.OutputCellPopulationBoundaryConditionParameters(plane_boundary_condition_parameter_file);
        plane_boundary_condition_parameter_file->close();

        std::string plane_boundary_condition_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + plane_boundary_condition_results_dir + "plane_results.parameters cell_based/test/data/TestCellBoundaryConditionsOutputParameters/plane_results.parameters").c_str()), 0);

        // Test with SphereGeometryBoundaryCondition
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> population(mesh, cells);
        c_vector<double,2> centre = zero_vector<double>(2);

        SphereGeometryBoundaryCondition<2> sphere_boundary_condition(&population, centre, 0.56, 1e-3);
        TS_ASSERT_EQUALS(sphere_boundary_condition.GetIdentifier(), "SphereGeometryBoundaryCondition-2");

        out_stream sphere_boundary_condition_parameter_file = output_file_handler.OpenOutputFile("sphere_results.parameters");
        sphere_boundary_condition.OutputCellPopulationBoundaryConditionParameters(sphere_boundary_condition_parameter_file);
        sphere_boundary_condition_parameter_file->close();

        std::string sphere_boundary_condition_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + sphere_boundary_condition_results_dir + "sphere_results.parameters cell_based/test/data/TestCellBoundaryConditionsOutputParameters/sphere_results.parameters").c_str()), 0);

        // Test OutputCellPopulationBoundaryConditionInfo() method
        out_stream plane_boundary_condition_info_file = output_file_handler.OpenOutputFile("plane_results.info");
        plane_boundary_condition.OutputCellPopulationBoundaryConditionInfo(plane_boundary_condition_info_file);
        plane_boundary_condition_info_file->close();

        TS_ASSERT_EQUALS(system(("diff " + plane_boundary_condition_results_dir + "plane_results.info cell_based/test/data/TestCellBoundaryConditionsOutputParameters/plane_results.info").c_str()), 0);
    }
};

#endif /*TESTCELLPOPULATIONBOUNDARYCONDITIONS_HPP_*/
