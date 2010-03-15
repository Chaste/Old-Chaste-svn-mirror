/*

Copyright (C) University of Oxford, 2005-2010

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

#ifndef CELLSGENERATOR_HPP_
#define CELLSGENERATOR_HPP_

#include "AbstractCellsGenerator.hpp"

template<class CELL_CYCLE_MODEL, unsigned DIM>
class CellsGenerator
{
public:

    CellsGenerator()
    {}

    ~CellsGenerator()
    {}

    void GenerateBasic(std::vector<TissueCell>& rCells,
                       unsigned numCells,
                       const std::vector<unsigned> locationIndices=std::vector<unsigned>())
    {
        rCells.clear();
        unsigned num_cells;
        if (!locationIndices.empty())
        {
            num_cells = locationIndices.size();
        }
        else
        {
            num_cells = numCells;
        }
        rCells.reserve(num_cells);

        // Create cells
        for (unsigned i=0; i<num_cells; i++)
        {
            CELL_CYCLE_MODEL* p_cell_cycle_model = new CELL_CYCLE_MODEL();
            TissueCell cell(STEM, HEALTHY, p_cell_cycle_model);

            double birth_time;
            if (!locationIndices.empty())
            {
                birth_time = 0.0 - locationIndices[i];
            }
            else
            {
                birth_time = 0.0 - i;
            }
            cell.SetBirthTime(birth_time);
            rCells.push_back(cell);
        }
    }
};



#endif /* CELLSGENERATOR_HPP_ */
