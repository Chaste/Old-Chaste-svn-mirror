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
#ifndef ABSTRACTCELLSGENERATOR_HPP_
#define ABSTRACTCELLSGENERATOR_HPP_

#include "TetrahedralMesh.hpp"
#include "TissueCell.hpp"

/**
 * A helper class for generating a vector of cells for a given mesh.
 * 
 * It is subclassed for different types of cell model.
 */
template<unsigned DIM>
class AbstractCellsGenerator
{
public:

    /**
     * Default constructor.
     */
    AbstractCellsGenerator()
    {}

    /**
     * @return a pointer to a new cell cycle model (should be implemented in subclasses)
     */
    virtual AbstractCellCycleModel* CreateCellCycleModel()=0;
    
    /**
     * Return the typical cell cycle duration for a transit cell, in hours.
     * Used when giving cells random ages - the ages will follow a uniform
     * distribution with this value as the upper limit.
     */
    virtual double GetTypicalTransitCellCycleTime()=0;
    
    /**
     * Return the typical cell cycle duration for a stem cell, in hours.
     * Used when giving cells random ages - the ages will follow a uniform
     * distribution with this value as the upper limit.
     */
    virtual double GetTypicalStemCellCycleTime()=0;
    
    /**
     * Whether cells are able to fully differentiate.
     * Defaults to false unless overridden.
     */
    virtual bool CellsCanDifferentiate();

    /**
     * Destructor.
     */    
    virtual ~AbstractCellsGenerator()
    {}

    /**
     * Generates cells of a specified cell cycle type under the correct
     * crypt conditions and gives random ages if required,
     * or gives them an age of 0.0 - creates least work for solver startup.
     *
     * @param rCells  An empty cells vector for this function to fill up
     * @param rMesh  The crypt mesh (can be cylindrical)
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
                                 bool randomBirthTimes,
                                 double y0 = 0.3,
                                 double y1 = 2.0,
                                 double y2 = 3.0,
                                 double y3 = 4.0,
                                 bool initialiseCells = false);
};


template<unsigned DIM>
bool AbstractCellsGenerator<DIM>::CellsCanDifferentiate()
{
    return false;
}

template<unsigned DIM>
void AbstractCellsGenerator<DIM>::GenerateForCrypt(std::vector<TissueCell>& rCells,
                                 TetrahedralMesh<2,2>& rMesh,
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
    
    RandomNumberGenerator *p_random_num_gen = RandomNumberGenerator::Instance();
    unsigned num_cells = rMesh.GetNumNodes();

    AbstractCellCycleModel* p_cell_cycle_model = NULL;
    double typical_transit_cycle_time;
    double typical_stem_cycle_time;

    rCells.clear();
    rCells.reserve(num_cells);

    for (unsigned i=0; i<num_cells; i++)
    {
        CellType cell_type;
        unsigned generation;

        double y = rMesh.GetNode(i)->GetPoint().rGetLocation()[1];

        p_cell_cycle_model = CreateCellCycleModel();
        typical_transit_cycle_time = this->GetTypicalTransitCellCycleTime();
        typical_stem_cycle_time = GetTypicalStemCellCycleTime();
        
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
        /**
         * \todo
         * Not all cell cycle models use the generation, only some of them are generation based.
         * The following line does not do any harm, but does give redundant information to some of the models.
         * We could create an AbstractGenerationBasedCellCycleModel, but not a priority at the moment (see #839).
         */
        p_cell_cycle_model->SetGeneration(generation);
        
        TissueCell cell(cell_type, HEALTHY, p_cell_cycle_model);
        if (initialiseCells)
        {
            cell.InitialiseCellCycleModel();
        }

        cell.SetBirthTime(birth_time);
        rCells.push_back(cell);
    }
}

#endif /*ABSTRACTCELLSGENERATOR_HPP_*/
