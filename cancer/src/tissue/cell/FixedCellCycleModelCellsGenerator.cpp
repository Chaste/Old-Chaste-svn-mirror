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
#include "FixedCellCycleModelCellsGenerator.hpp"

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
        cell.SetLocationIndex(i);
        cell.SetBirthTime(birth_time);
        rCells.push_back(cell);
    }
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class FixedCellCycleModelCellsGenerator<1>;
template class FixedCellCycleModelCellsGenerator<2>;
template class FixedCellCycleModelCellsGenerator<3>;
