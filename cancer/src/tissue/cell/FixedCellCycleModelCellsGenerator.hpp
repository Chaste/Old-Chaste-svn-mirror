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
#ifndef FIXEDCELLCYCLEMODELCELLSGENERATOR_HPP_
#define FIXEDCELLCYCLEMODELCELLSGENERATOR_HPP_

#include "AbstractCellsGenerator.hpp"
#include "FixedCellCycleModel.hpp"

/**
 * A helper class for generating a vector of cells with 
 * FixedCellCycleModels for a given mesh.
 */
template<unsigned DIM>
class FixedCellCycleModelCellsGenerator : public AbstractCellsGenerator<DIM>
{
public:

    /**
     * @return a pointer to a new FixedCellCycleModel.
     */
    AbstractCellCycleModel* CreateCellCycleModel();
    
    /**
     * Fills a vector of cells with a specified cell cycle model, to match
     * a given mesh. Gives them birth times of 0 for node 0,
     * -1 for node 1, -2 for node 2 etc...
     *
     * @param rCells  An empty vector of cells to fill up.
     * @param rMesh  The mesh the cells should be associated with.
     */
    void GenerateBasic(std::vector<TissueCell>& rCells,
                               TetrahedralMesh<DIM,DIM>& rMesh);

    /**
     * @return default cell cycle time for a transit cell.
     */
    double GetTypicalTransitCellCycleTime();
    
    /**
     * @return default cell cycle time for a transit cell.
     */
    double GetTypicalStemCellCycleTime();

    /**
     * @return true (cells can always differentiate).
     */
    bool CellsCanDifferentiate();
};


#endif /*FIXEDCELLCYCLEMODELCELLSGENERATOR_HPP_*/
