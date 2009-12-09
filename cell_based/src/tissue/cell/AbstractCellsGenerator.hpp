/*

Copyright (C) University of Oxford, 2005-2009

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
#include "AbstractSimpleGenerationBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"


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
     * Fills a vector of cells with a specified cell cycle model, to match
     * a given number of cells. Gives them birth times of 0 for node 0,
     * -1 for node 1, -2 for node 2 etc...
     *
     * @param rCells  An empty vector of cells to fill up.
     * @param numCells  The number of cells to generate.
     */
    void GenerateBasic(std::vector<TissueCell>& rCells,
                       const unsigned numCells);

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
};



#endif /*ABSTRACTCELLSGENERATOR_HPP_*/
