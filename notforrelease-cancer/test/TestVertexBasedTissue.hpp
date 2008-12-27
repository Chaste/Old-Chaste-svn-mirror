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

//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>

#include "VertexBasedTissue.hpp"
//#include "MeinekeInteractionForce.hpp"
//#include "HoneycombMeshGenerator.hpp"
#include "FixedCellCycleModel.hpp"
//#include "OutputFileHandler.hpp"
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
        
        // Create the tissue
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

        // Fails as no cell or ghost correponding to node 0
        TS_ASSERT_THROWS_ANYTHING(VertexBasedTissue<2> tissue(mesh, cells));
    }
    
};


#endif /*TESTVERTEXBASEDTISSUE_HPP_*/
