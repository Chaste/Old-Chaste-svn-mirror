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
 * A helper class for generating a vector of cells for a given mesh
 */
template<unsigned DIM>
class FixedCellCycleModelCellsGenerator : public AbstractCellsGenerator<DIM>
{
public:

    virtual AbstractCellCycleModel* CreateCellCycleModel();
    
    /**
     * Fills a vector of cells with a specified cell cycle model, to match
     * a given mesh. Gives them birth times of 0 for node 0,
     * -1 for node 1, -2 for node 2 etc...
     *
     * @param rCells  An empty vector of cells to fill up.
     * @param rMesh  The mesh the cells should be associated with.
     */
    virtual void GenerateBasic(std::vector<TissueCell>& rCells,
                               TetrahedralMesh<DIM,DIM>& rMesh);
    
    virtual double GetTypicalTransitCellCycleTime();
    
    virtual double GetTypicalStemCellCycleTime();

    virtual bool CellsCanDifferentiate();
};


template<unsigned DIM>
AbstractCellCycleModel* FixedCellCycleModelCellsGenerator<DIM>::CreateCellCycleModel()
{
    return new FixedCellCycleModel();
}

template<unsigned DIM>
double FixedCellCycleModelCellsGenerator<DIM>::GetTypicalTransitCellCycleTime()
{
    return CancerParameters::Instance()->GetTransitCellG1Duration()
            + CancerParameters::Instance()->GetSG2MDuration();
}

template<unsigned DIM>
double FixedCellCycleModelCellsGenerator<DIM>::GetTypicalStemCellCycleTime()
{
    return CancerParameters::Instance()->GetStemCellG1Duration()
            + CancerParameters::Instance()->GetSG2MDuration();
}

template<unsigned DIM>
bool FixedCellCycleModelCellsGenerator<DIM>::CellsCanDifferentiate()
{
    return true;
}

template<unsigned DIM>
void FixedCellCycleModelCellsGenerator<DIM>::GenerateBasic(
    std::vector<TissueCell>& rCells,
    TetrahedralMesh<DIM,DIM>& rMesh)
{
    rCells.clear();
    rCells.reserve(rMesh.GetNumNodes());
    for (unsigned i=0; i<rMesh.GetNumNodes(); i++)
    {
        AbstractCellCycleModel* p_cell_cycle_model = CreateCellCycleModel();
        TissueCell cell(STEM, HEALTHY, p_cell_cycle_model);        
        double birth_time = 0.0 - i;
        cell.SetNodeIndex(i);
        cell.SetBirthTime(birth_time);
        rCells.push_back(cell);
    }
}


#endif /*FIXEDCELLCYCLEMODELCELLSGENERATOR_HPP_*/
