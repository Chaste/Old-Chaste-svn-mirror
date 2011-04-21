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
#ifndef TESTFAILINGDIVISION3D_HPP_
#define TESTFAILINGDIVISION3D_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "TrianglesMeshReader.hpp"
#include "CellBasedSimulation.hpp"
#include "TrianglesMeshWriter.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"

class TestFailingDivision3d : public AbstractCellBasedTestSuite
{
private:
    
    MutableMesh<3,3>* Make3dMesh(unsigned width=3, unsigned height=3, unsigned depth=3)
    {
        MutableMesh<3,3>* p_mesh = new MutableMesh<3,3>;
        p_mesh->ConstructCuboid(width, height, depth);
        TrianglesMeshWriter<3,3> mesh_writer("", "3dSpringMesh");
        mesh_writer.WriteFilesUsingMesh(*p_mesh);

        return p_mesh;
    }
    
public:

    /* This test will pass if you run it for 10 hours, but fails for longer times when a cell tries to divide
     *
     */
    void TestFailsWhenCellDivides() throw (Exception)
    {
        unsigned width = 4;
        unsigned height = 4;
        unsigned depth = 1;

        MutableMesh<3,3>* p_mesh = Make3dMesh(width, height, depth);

        // Set up cells by iterating through the mesh nodes
        unsigned num_nodes = p_mesh->GetNumAllNodes();

        unsigned num_epithelial_cells = (width+1)*(height+1);

    	boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

        std::vector<CellPtr> cells;
		for (unsigned i=0; i<num_nodes-num_epithelial_cells; i++)
		{
			FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            p_model->SetDimension(3);

            CellPtr p_differentiated_cell(new Cell(p_state, p_model));
			cells.push_back(p_differentiated_cell);
        }

        for (unsigned i=num_nodes-num_epithelial_cells; i<num_nodes; i++)
        {
			FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(TRANSIT);
            p_model->SetDimension(3);

            CellPtr p_epithelial_cell(new Cell(p_state, p_model));
            p_epithelial_cell->InitialiseCellCycleModel();
            
            p_epithelial_cell->SetBirthTime(-10.0);

			cells.push_back(p_epithelial_cell);
        }

        // Test Save with a MeshBasedCellPopulationWithGhostNodes
        MeshBasedCellPopulation<3> cell_population(*p_mesh, cells);

        CellBasedSimulation<3> simulator(cell_population);

        simulator.SetOutputDirectory("TestFailsWhenCellDivides");
        simulator.SetEndTime(20.0);			// Keep at 10 hours max and no cell divisions will occur. Cell divisions don't seem to work yet.

        // Create a force law and pass it to the simulation
        GeneralisedLinearSpringForce<3> linear_force;
        linear_force.SetCutOffLength(1.5);
        simulator.AddForce(&linear_force);

        simulator.Solve();
    }

};

#endif /* TESTFAILINGDIVISION3D_HPP_ */
