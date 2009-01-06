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

#ifndef TESTVERTEXBASEDTISSUE_HPP_
#define TESTVERTEXBASEDTISSUE_HPP_

#include <cxxtest/TestSuite.h>

#include "VertexBasedTissue.hpp"
#include "FixedCellCycleModel.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestBasedTissue : public AbstractCancerTestSuite
{

public:

    // Test construction, accessors and iterator
    void TestCreateSmallVertexBasedTissue() throw(Exception)
    {
        // Create a simple 2D VertexMesh
        VertexMesh<2,2> mesh(5,3); // columns then rows
        
        // Set up cells, one for each VertexElement. Give each cell 
        // a birth time of -elem_index, so its age is elem_index
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TissueCell cell(DIFFERENTIATED, HEALTHY, new FixedCellCycleModel());
            double birth_time = 0.0-i;
            cell.SetLocationIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);
        
        // Test we have the correct number of cells and elements
        TS_ASSERT_EQUALS(tissue.rGetMesh().GetNumElements(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), cells.size());

        unsigned counter = 0;

        // Test VertexBasedTissue::Iterator
        for (VertexBasedTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            // Test operator* and that cells are in sync
            TS_ASSERT_EQUALS((*cell_iter).GetLocationIndex(), counter);

            // Test operator-> and that cells are in sync
            TS_ASSERT_DELTA(cell_iter->GetAge(), (double)counter, 1e-12);

            // Test GetElement on the iterator
            TS_ASSERT_EQUALS(cell_iter.GetElement()->GetIndex(), mesh.GetElement(counter)->GetIndex());

            // Test iter.GetElement()->GetIndex() is consistent with cell.GetLocationIndex()
            TS_ASSERT_EQUALS((*cell_iter).GetLocationIndex(), cell_iter.GetElement()->GetIndex());

            counter++;
        }

        // Test we have gone through all cells in the for loop
        TS_ASSERT_EQUALS(counter, tissue.GetNumRealCells());
        
        // Test GetNumNodes() method
        TS_ASSERT_EQUALS(tissue.GetNumNodes(), mesh.GetNumNodes());
    }
    
    void TestValidateVertexBasedTissue()
    {
        // Create a simple vertex-based mesh
        VertexMesh<2,2> mesh(4,6); // columns then rows

        // Set up cells, one for each VertexElement. Give each cell 
        // a birth time of -elem_index, so its age is elem_index
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TissueCell cell(DIFFERENTIATED, HEALTHY, new FixedCellCycleModel());
            double birth_time = 0.0-i;
            cell.SetLocationIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        cells[0].SetLocationIndex(1);

        // This test fails as there is no cell to element 0
        TS_ASSERT_THROWS_ANYTHING(VertexBasedTissue<2> tissue(mesh, cells));
    }
    
    void TestRemoveDeadCellsAndUpdate()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);
        
        // Create a simple vertex-based mesh
        VertexMesh<2,2> mesh(4,6); // columns then rows

        // Set up cells, one for each VertexElement. Give each cell 
        // a birth time of -elem_index, so its age is elem_index
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TissueCell cell(DIFFERENTIATED, HEALTHY, new FixedCellCycleModel());
            double birth_time = 0.0-i;
            cell.SetLocationIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        cells[5].StartApoptosis();

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 24u);        
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 24u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 24u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 68u);

        p_simulation_time->IncrementTimeOneStep();
        
        // Remove dead cells
        unsigned num_cells_removed = tissue.RemoveDeadCells();
        
        /// \todo Currently RemoveDeadCells() does nothing, and is only 
        //        in a test for coverage. Cell death will be implemented 
        //        in #853.
        TS_ASSERT_EQUALS(num_cells_removed, 0u);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 24u);        
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 24u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 24u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 68u);
                
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), cells.size()); // the tissue now copies cells

        tissue.Update();
        
        // Finally, check the cells node indices have updated

        // We expect the cell node indices to be {0,11,...,23}
        std::set<unsigned> expected_node_indices;
        for (unsigned i=0; i<tissue.GetNumRealCells(); i++)
        {
            expected_node_indices.insert(i);
        }

        // Get actual cell node indices
        std::set<unsigned> node_indices;

        for (VertexBasedTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            // Record node index corresponding to cell
            unsigned node_index = tissue.GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
            node_indices.insert(node_index);
        }

        TS_ASSERT_EQUALS(node_indices, expected_node_indices);        
    }

    
};


#endif /*TESTVERTEXBASEDTISSUE_HPP_*/
