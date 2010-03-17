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
 * It is templated for different types of cell model.
 */
template<class CELL_CYCLE_MODEL, unsigned DIM>
class CellsGenerator
{
public:

    /**
     * Return the typical cell cycle duration for a transit cell, in hours.
     * Used when giving cells random ages - the ages will follow a uniform
     * distribution with this value as the upper limit.
     */
    double GetTypicalTransitCellCycleTime();

    /**
     * Return the typical cell cycle duration for a stem cell, in hours.
     * Used when giving cells random ages - the ages will follow a uniform
     * distribution with this value as the upper limit.
     */
    double GetTypicalStemCellCycleTime();

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
     *
     */
    virtual void GenerateForCrypt(std::vector<TissueCell>& rCells,
                                  TetrahedralMesh<2,2>& rMesh,
                                  const std::vector<unsigned> locationIndices,
                                  bool randomBirthTimes,
                                  double y0 = 0.3,
                                  double y1 = 2.0,
                                  double y2 = 3.0,
                                  double y3 = 4.0,
                                  bool initialiseCells = false);

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
void CellsGenerator<CELL_CYCLE_MODEL,DIM>::GenerateForCrypt(std::vector<TissueCell>& rCells,
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

        double typical_transit_cycle_time = GetTypicalTransitCellCycleTime();
        double typical_stem_cycle_time = GetTypicalStemCellCycleTime();

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
            cell_type = CellsCanDifferentiate() ? DIFFERENTIATED : TRANSIT;
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

///////////////////////////////
//
//template<unsigned DIM>
//double CellsGenerator<IngeWntSwatCellCycleModel,DIM>::GetTypicalTransitCellCycleTime()
//{
//    return 16.0;
//}
//
//template<unsigned DIM>
//double CellsGenerator<IngeWntSwatCellCycleModel,DIM>::GetTypicalStemCellCycleTime()
//{
//    return 16.0;
//}
//

//template<unsigned DIM>
//double CellsGenerator<SingleOdeWntCellCycleModel,DIM>::GetTypicalTransitCellCycleTime()
//{
//    return 16.0;
//}
//
//
//template<unsigned DIM>
//double CellsGenerator<SingleOdeWntCellCycleModel,DIM>::GetTypicalStemCellCycleTime()
//{
//    return 16.0;
//}
//
//
//template<unsigned DIM>
//double CellsGenerator<StochasticWntCellCycleModel,DIM>::GetTypicalTransitCellCycleTime()
//{
//    return 16.0;
//}
//
//
//template<unsigned DIM>
//double CellsGenerator<StochasticWntCellCycleModel,DIM>::GetTypicalStemCellCycleTime()
//{
//    return 16.0;
//}
//
//
//


template<class CELL_CYCLE_MODEL, unsigned DIM>
double CellsGenerator<CELL_CYCLE_MODEL, DIM>::GetTypicalTransitCellCycleTime()
{
	// Here the returned value will depend on what type of model CELL_CYCLE_MODEL
	// is. The following madness just basically does
	//
	// if(CELL_CYCLE_MODEL==TysonNovakCellCycleModel)
	//   return ...
	// else if (CELL_CYCLE_MODEL==WntCellCycleModel)
	//   return ..
	// else if ...
	//
	using namespace boost::mpl;
	using namespace boost;
    typedef typename if_<is_same<CELL_CYCLE_MODEL, TysonNovakCellCycleModel>,
					     integral_c<unsigned, 1>,
					     typename if_<is_same<CELL_CYCLE_MODEL, WntCellCycleModel>,
									  integral_c<unsigned, 2>,
									  typename if_<is_same<CELL_CYCLE_MODEL, StochasticWntCellCycleModel>,
											       integral_c<unsigned, 2>,
											       typename if_<is_same<CELL_CYCLE_MODEL, SingleOdeWntCellCycleModel>,
															    integral_c<unsigned, 2>,
															    typename if_<is_same<CELL_CYCLE_MODEL, IngeWntSwatCellCycleModel>,
																     		 integral_c<unsigned, 2>,
																     		 integral_c<unsigned, 0>
																     		>::type
															   >::type
											      >::type
									  >::type
					    >::type selector_t;
    unsigned selector = selector_t();

    switch (selector)
    {
		case 1:
			return 1.25;
		case 2:
			return 16.0;
		default:
			return TissueConfig::Instance()->GetTransitCellG1Duration()
				+ TissueConfig::Instance()->GetSG2MDuration();
    }
}


template<class CELL_CYCLE_MODEL, unsigned DIM>
double CellsGenerator<CELL_CYCLE_MODEL, DIM>::GetTypicalStemCellCycleTime()
{
	// Here the returned value will depend on what type of model CELL_CYCLE_MODEL
	// is. The following madness just basically does
	//
	// if(CELL_CYCLE_MODEL==TysonNovakCellCycleModel)
	//   return ...
	// else if (CELL_CYCLE_MODEL==WntCellCycleModel)
	//   return ..
	// else if ...
	//
	using namespace boost::mpl;
	using namespace boost;
    typedef typename if_<is_same<CELL_CYCLE_MODEL, TysonNovakCellCycleModel>,
					     integral_c<unsigned, 1>,
					     typename if_<is_same<CELL_CYCLE_MODEL, WntCellCycleModel>,
									  integral_c<unsigned, 2>,
									  typename if_<is_same<CELL_CYCLE_MODEL, StochasticWntCellCycleModel>,
											       integral_c<unsigned, 2>,
											       typename if_<is_same<CELL_CYCLE_MODEL, SingleOdeWntCellCycleModel>,
															    integral_c<unsigned, 2>,
															    typename if_<is_same<CELL_CYCLE_MODEL, IngeWntSwatCellCycleModel>,
																     		 integral_c<unsigned, 2>,
																     		 integral_c<unsigned, 0>
																     		>::type
															   >::type
											      >::type
									  >::type
					    >::type selector_t;
    unsigned selector = selector_t();
    switch (selector)
    {
		case 1:
			return 1.25;
		case 2:
			return 16.0;
		default:
			return TissueConfig::Instance()->GetStemCellG1Duration()
				+ TissueConfig::Instance()->GetSG2MDuration();
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

    if(selector==1)
    {
    	// With FixedDuration or Stochastic cell cycle models, cells can differentiate
    	return true;
    }
    else
    {
    	return false;
    }
}


//typedef typename if_<is_same<CELL_CYCLE_MODEL, TysonNovakCellCycleModel>,integral_c<unsigned, 1>,typename if_<is_same<CELL_CYCLE_MODEL, WntCellCycleModel>,integral_c<unsigned, 2>,typename if_<is_same<CELL_CYCLE_MODEL, StochasticWntCellCycleModel>,integral_c<unsigned, 2>,typename if_<is_same<CELL_CYCLE_MODEL, SingleOdeWntCellCycleModel>,integral_c<unsigned, 2>,typename if_<is_same<CELL_CYCLE_MODEL, IngeWntSwatCellCycleModel>,integral_c<unsigned, 2>,integral_c<unsigned, 0>	>::type>::type>::type>::type>::type selector_t;


#endif /* CELLSGENERATOR_HPP_ */
