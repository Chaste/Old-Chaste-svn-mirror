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


#include <vector>
#include "TissueCell.hpp"
#include "WildTypeCellMutationState.hpp"

/**
 * A helper class for generating a vector of cells for a given mesh.
 *
 * It is templated over types of cell cycle model and spatial dimension.
 */
template<class CELL_CYCLE_MODEL, unsigned DIM>
class CellsGenerator
{
public:

    /**
     * Fills a vector of cells with a specified cell cycle model, to match
     * a given number of cells. Gives them birth times of 0 for node 0,
     * -1 for node 1, -2 for node 2 etc...
     *
     * @param rCells  An empty vector of cells to fill up.
     * @param numCells  The number of cells to generate.
     * @param locationIndices is used when a birth-time hint is needed for individual cell.
     *             Defaults to an empty vector -- otherwise must be of length numCells
     *
     */
    void GenerateBasic(std::vector<TissueCellPtr>& rCells,
                       unsigned numCells,
                       const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Fills a vector of cells with birth times to match a given vector of location indices.
     *
     * @param rCells  An empty vector of cells to fill up.
     * @param locationIndices  The indices of the tissue to assign real cells to.
     */
    void GenerateGivenLocationIndices(std::vector<TissueCellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices);

};

template<class CELL_CYCLE_MODEL, unsigned DIM>
void CellsGenerator<CELL_CYCLE_MODEL,DIM>::GenerateBasic(std::vector<TissueCellPtr>& rCells,
                                                         unsigned numCells,
                                                         const std::vector<unsigned> locationIndices)
{
    CellPropertyRegistry::Instance()->Clear();
    rCells.clear();
    if (!locationIndices.empty())
    {
        // If location indices is given, then it needs to match the number of output cells
        if (numCells != locationIndices.size())
        {
            EXCEPTION("The size of the locationIndices vector must match the required number of output cells");
        }
    }
    rCells.reserve(numCells);

    // Create cells
    for (unsigned i=0; i<numCells; i++)
    {
        CELL_CYCLE_MODEL* p_cell_cycle_model = new CELL_CYCLE_MODEL;
        p_cell_cycle_model->SetDimension(DIM);
        p_cell_cycle_model->SetCellProliferativeType(STEM);

        boost::shared_ptr<AbstractCellMutationState> p_state = boost::dynamic_pointer_cast<AbstractCellMutationState>(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        TissueCellPtr p_cell(new TissueCell(p_state, p_cell_cycle_model));

        double birth_time;
        if (!locationIndices.empty())
        {
            birth_time = 0.0 - locationIndices[i];
        }
        else
        {
            birth_time = 0.0 - i;
        }
        p_cell->SetBirthTime(birth_time);
        rCells.push_back(p_cell);
    }
}

template<class CELL_CYCLE_MODEL, unsigned DIM>
void CellsGenerator<CELL_CYCLE_MODEL,DIM>::GenerateGivenLocationIndices(std::vector<TissueCellPtr>& rCells,
                                                                        const std::vector<unsigned> locationIndices)
{
    assert(!locationIndices.empty());

    unsigned num_cells = locationIndices.size();

    rCells.clear();
    rCells.reserve(num_cells);
    CellPropertyRegistry::Instance()->Clear();

    for (unsigned i=0; i<num_cells; i++)
    {
        CELL_CYCLE_MODEL* p_cell_cycle_model = new CELL_CYCLE_MODEL;
        p_cell_cycle_model->SetDimension(DIM);
        p_cell_cycle_model->SetCellProliferativeType(STEM);

        boost::shared_ptr<AbstractCellMutationState> p_state = boost::dynamic_pointer_cast<AbstractCellMutationState>(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        TissueCellPtr p_cell(new TissueCell(p_state, p_cell_cycle_model));

        double birth_time = 0.0 - locationIndices[i];
        p_cell->SetBirthTime(birth_time);
        rCells.push_back(p_cell);
    }
}

#endif /* CELLSGENERATOR_HPP_ */
