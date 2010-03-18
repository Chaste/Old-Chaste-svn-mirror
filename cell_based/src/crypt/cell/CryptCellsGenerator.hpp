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
#ifndef CRYPTCELLSGENERATOR_HPP_
#define CRYPTCELLSGENERATOR_HPP_

#include "CellsGenerator.hpp"

/**
 * A subclass of CellsGenerator that generates
 * cells for crypt simulations.
 *
 * It is templated over types of cell cycle model.
 */
template<class CELL_CYCLE_MODEL, unsigned DIM>
class CryptCellsGenerator : public CellsGenerator<CELL_CYCLE_MODEL,DIM>
{
public:
    /**
     * Generates cells of a specified cell cycle type under the correct
     * crypt conditions and gives random ages if required,
     * or gives them an age of 0.0 - creates least work for solver startup.
     *
     * @param rCells  An empty cells vector for this function to fill up
     * @param rMesh  The crypt mesh (can be cylindrical)
     * @param locationIndices the node indices corresponding to real cells
     * @param randomBirthTimes  Whether to assign the cells random birth times
     *    (this can be expensive computationally with ODE models)
     * @param y0  below this line cells are generation 0 (defaults to 0.3)
     * @param y1  below this line cells are generation 1 (defaults to 2.0)
     * @param y2  below this line cells are generation 2 (defaults to 3.0)
     * @param y3  below this line cells are generation 3 (defaults to 4.0)
     * @param initialiseCells  whether to initialise the cell cycle models as each
     *   cell is created
     */
    void Generate(std::vector<TissueCell>& rCells,
                  TetrahedralMesh<2,2>& rMesh,
                  const std::vector<unsigned> locationIndices,
                  bool randomBirthTimes,
                  double y0 = 0.3,
                  double y1 = 2.0,
                  double y2 = 3.0,
                  double y3 = 4.0,
                  bool initialiseCells = false);
};


template<class CELL_CYCLE_MODEL, unsigned DIM>
void CryptCellsGenerator<CELL_CYCLE_MODEL,DIM>::Generate(
									  std::vector<TissueCell>& rCells,
									  TetrahedralMesh<2,2>& rMesh,
									  const std::vector<unsigned> locationIndices,
									  bool randomBirthTimes,
									  double y0,
									  double y1,
									  double y2,
									  double y3,
									  bool initialiseCells)
{
    #define COVERAGE_IGNORE
    assert(DIM==2);
    #undef COVERAGE_IGNORE

    RandomNumberGenerator* p_random_num_gen = RandomNumberGenerator::Instance();

    unsigned num_cells = locationIndices.empty() ? rMesh.GetNumNodes() : locationIndices.size();

    rCells.clear();
    rCells.reserve(num_cells);

    for (unsigned i=0; i<rMesh.GetNumNodes(); i++)
    {
    	CellProliferativeType cell_type;
        unsigned generation;

        double y = 0.0;
        if (!locationIndices.empty())
        {
            if ( std::find(locationIndices.begin(), locationIndices.end(), i) != locationIndices.end() )
            {
                y = rMesh.GetNode(i)->GetPoint().rGetLocation()[1];
            }
        }
        else
        {
            y = rMesh.GetNode(i)->GetPoint().rGetLocation()[1];
        }

        CELL_CYCLE_MODEL* p_cell_cycle_model = new CELL_CYCLE_MODEL;
        p_cell_cycle_model->SetDimension(DIM);

        double typical_transit_cycle_time = p_cell_cycle_model->GetAverageTransitCellCycleTime();
        double typical_stem_cycle_time = p_cell_cycle_model->GetAverageStemCellCycleTime();

        double birth_time = 0.0;
        if (randomBirthTimes)
        {
            birth_time = -p_random_num_gen->ranf();
        }

        if (y <= y0)
        {
            cell_type = STEM;
            generation = 0;
            birth_time *= typical_stem_cycle_time; // hours
        }
        else if (y < y1)
        {
            cell_type = TRANSIT;
            generation = 1;
            birth_time *= typical_transit_cycle_time; // hours
        }
        else if (y < y2)
        {
            cell_type = TRANSIT;
            generation = 2;
            birth_time *= typical_transit_cycle_time; // hours
        }
        else if (y < y3)
        {
            cell_type = TRANSIT;
            generation = 3;
            birth_time *= typical_transit_cycle_time; // hours
        }
        else
        {
            cell_type = this->CellsCanDifferentiate() ? DIFFERENTIATED : TRANSIT;
            generation = 4;
            birth_time *= typical_transit_cycle_time; // hours
        }

        if (dynamic_cast<AbstractSimpleGenerationBasedCellCycleModel*>(p_cell_cycle_model))
        {
        	dynamic_cast<AbstractSimpleGenerationBasedCellCycleModel*>(p_cell_cycle_model)->SetGeneration(generation);
        }

        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        TissueCell cell(cell_type, p_state, p_cell_cycle_model);
        if (initialiseCells)
        {
            cell.InitialiseCellCycleModel();
        }

        cell.SetBirthTime(birth_time);

        if (!locationIndices.empty())
        {
            if ( std::find(locationIndices.begin(), locationIndices.end(), i) != locationIndices.end() )
            {
                rCells.push_back(cell);
            }
        }
        else
        {
            rCells.push_back(cell);
        }
    }
}

#endif /* CRYPTCELLSGENERATOR_HPP_ */
