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

//Shoot me now, and put me out of my misery
#include <boost/mpl/integral_c.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>

#include <vector>
#include "TissueCell.hpp"


#include "IngeWntSwatCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "SingleOdeWntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"

/**
 * A helper class for generating a vector of cells for a given mesh.
 * \todo write a generator for meshes
 *
 * It is templated over types of cell cycle model.
 */
template<class CELL_CYCLE_MODEL, unsigned DIM>
class CellsGenerator
{
public:

    /**
     * Whether cells are able to fully differentiate.
     * Defaults to false unless overridden.
     */
    bool CellsCanDifferentiate();

    /**
     * Fills a vector of cells with a specified cell cycle model, to match
     * a given number of cells. Gives them birth times of 0 for node 0,
     * -1 for node 1, -2 for node 2 etc...
     *
     * @param rCells  An empty vector of cells to fill up.
     * @param numCells  The number of cells to generate.
     * @param locationIndices is used when a birth-time hint is needed for individual cell.
     * 			Defaults to an empty vector -- otherwise must be of length numCells
     *
     */
    void GenerateBasic(std::vector<TissueCell>& rCells,
                       unsigned numCells,
                       const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Fills a vector of cells with birth times to match a given vector of location indices.
     *
     * @param rCells  An empty vector of cells to fill up.
     * @param locationIndices  The indices of the tissue to assign real cells to.
     */
    void GenerateGivenLocationIndices(std::vector<TissueCell>& rCells,
    		                          const std::vector<unsigned> locationIndices);

};

template<class CELL_CYCLE_MODEL, unsigned DIM>
void CellsGenerator<CELL_CYCLE_MODEL,DIM>::GenerateBasic(std::vector<TissueCell>& rCells,
											             unsigned numCells,
											             const std::vector<unsigned> locationIndices)
{
	rCells.clear();
	if (!locationIndices.empty())
	{
		//If location indices is given, then it needs to match the number of output cells
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

		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		TissueCell cell(STEM, p_state, p_cell_cycle_model);

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

template<class CELL_CYCLE_MODEL, unsigned DIM>
void CellsGenerator<CELL_CYCLE_MODEL,DIM>::GenerateGivenLocationIndices(std::vector<TissueCell>& rCells,
                                                                        const std::vector<unsigned> locationIndices)
{
    assert(!locationIndices.empty());

    unsigned num_cells = locationIndices.size();

    rCells.clear();
    rCells.reserve(num_cells);

    for (unsigned i=0; i<num_cells; i++)
    {
    	CELL_CYCLE_MODEL* p_cell_cycle_model = new CELL_CYCLE_MODEL;
        p_cell_cycle_model->SetDimension(DIM);

    	boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

    	TissueCell cell(STEM, p_state, p_cell_cycle_model);

    	double birth_time = 0.0 - locationIndices[i];
        cell.SetBirthTime(birth_time);
        rCells.push_back(cell);
    }
}

template<class CELL_CYCLE_MODEL, unsigned DIM>
bool CellsGenerator<CELL_CYCLE_MODEL, DIM>::CellsCanDifferentiate()
{
	using namespace boost::mpl;
	using namespace boost;
	typedef typename if_<is_same<CELL_CYCLE_MODEL, FixedDurationGenerationBasedCellCycleModel>,
					     integral_c<unsigned, 1>,
					     typename if_<is_same<CELL_CYCLE_MODEL, StochasticDurationGenerationBasedCellCycleModel>,
									  integral_c<unsigned, 1>,
									  integral_c<unsigned, 0>
									  >::type
     					  >::type selector_t;

    unsigned selector = selector_t();

    if (selector==1)
    {
    	// With FixedDuration or Stochastic cell cycle models, cells can differentiate
    	return true;
    }
    else
    {
    	return false;
    }
}

#endif /* CELLSGENERATOR_HPP_ */
