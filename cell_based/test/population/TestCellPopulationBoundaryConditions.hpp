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
#include "PlaneBoundaryCondition.hpp"
#include "CellwiseData.hpp"
#include "TrianglesMeshReader.hpp"
#include "WildTypeCellMutationState.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "NodeBasedCellPopulation.hpp"

/**
 * This class contains tests for methods on classes
 * inheriting from AbstractCellPopulationBoundaryConditions.
 */
class TestCellPopulationBoundaryConditions : public AbstractCellBasedTestSuite
{
public:

    void TestPlaneBoundaryConditionWithNodeBasedCellPopulation() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        TS_ASSERT(p_mesh->GetNumNodes()>0);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell population BCS
        c_vector<double,2> point = zero_vector<double>(2);
        point(0) = 2.0;
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        PlaneBoundaryCondition<2> boundary_condition(&cell_population, point, normal);

        TS_ASSERT_EQUALS(boundary_condition.GetIdentifier(), "PlaneBoundaryCondition-2");

        // Impose BCS
        std::vector<c_vector<double,2> > old_locations;
        old_locations.reserve(cell_population.GetNumNodes());
        for (std::list<CellPtr>::iterator cell_iter = cell_population.rGetCells().begin();
             cell_iter != cell_population.rGetCells().end();
             ++cell_iter)
        {
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            old_locations.push_back(p_node->rGetLocation());
        }

        boundary_condition.ImposeBoundaryConditions(old_locations);

        // Test that all nodes satisfy the boundary condition
        for (std::list<CellPtr>::iterator cell_iter = cell_population.rGetCells().begin();
             cell_iter != cell_population.rGetCells().end();
             ++cell_iter)
        {
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location = p_node->rGetLocation();
            if (old_locations[p_node->GetIndex()][1]<2.0)
            {
                TS_ASSERT_LESS_THAN_EQUALS(2.0, location[0]);
                TS_ASSERT_DELTA(location[1],old_locations[p_node->GetIndex()][1], 1e-6);
            }
            else
            {
                TS_ASSERT_DELTA(location[0], old_locations[p_node->GetIndex()][0], 1e-6);
                TS_ASSERT_DELTA(location[1], old_locations[p_node->GetIndex()][1], 1e-6);
            }
        }

        TS_ASSERT(boundary_condition.VerifyBoundaryConditions());
    }

    void TestArchivingOfPlaneBoundaryCondition() throw (Exception)
	{
	    // Set up singleton classes
	    OutputFileHandler handler("archive", false);    // don't erase contents of folder
	    std::string archive_filename = handler.GetOutputDirectoryFullPath() + "single_boundary_conditon.arch";

	    {
	    	// Create an output archive
	        PlaneBoundaryCondition<2> boundary_condition(NULL, zero_vector<double>(2), unit_vector<double>(2,1));

		    std::ofstream ofs(archive_filename.c_str());
		    boost::archive::text_oarchive output_arch(ofs);

		    // Serialize via pointer
		    PlaneBoundaryCondition<2>* const p_boundary_condition = &boundary_condition;
		    output_arch << p_boundary_condition;

            TS_ASSERT_DELTA(p_boundary_condition->rGetPointOnPlane()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(p_boundary_condition->rGetPointOnPlane()[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(p_boundary_condition->rGetNormalToPlane()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(p_boundary_condition->rGetNormalToPlane()[1], 1.0, 1e-6);
	    }

	    {
		    // Create an input archive
		    std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
		    boost::archive::text_iarchive input_arch(ifs);

		    PlaneBoundaryCondition<2>* p_boundary_condition;

		    // Restore from the archive
		    input_arch >> p_boundary_condition;

		    // Test we have restored the plane geometry correctly
            TS_ASSERT_DELTA(p_boundary_condition->rGetPointOnPlane()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(p_boundary_condition->rGetPointOnPlane()[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(p_boundary_condition->rGetNormalToPlane()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(p_boundary_condition->rGetNormalToPlane()[1], 1.0, 1e-6);

		    delete p_boundary_condition;
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
    }
};

#endif /*TESTCELLPOPULATIONBOUNDARYCONDITIONS_HPP_*/
