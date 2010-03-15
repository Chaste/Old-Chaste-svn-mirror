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
#include "IngeWntSwatCellCycleModelCellsGenerator.hpp"


template<unsigned DIM>
IngeWntSwatCellCycleModelCellsGenerator<DIM>::IngeWntSwatCellCycleModelCellsGenerator(unsigned hypothesis)
    : mHypothesis(hypothesis)
{
}


template<unsigned DIM>
AbstractCellCycleModel* IngeWntSwatCellCycleModelCellsGenerator<DIM>::CreateCellCycleModel()
{
	IngeWntSwatCellCycleModel* p_cell_cycle_model = new IngeWntSwatCellCycleModel();
	p_cell_cycle_model->SetDimension(DIM);
	p_cell_cycle_model->SetHypothesis(mHypothesis);
    return p_cell_cycle_model;
}


template<unsigned DIM>
double IngeWntSwatCellCycleModelCellsGenerator<DIM>::GetTypicalTransitCellCycleTime()
{
    return 16.0;
}

template<unsigned DIM>
double IngeWntSwatCellCycleModelCellsGenerator<DIM>::GetTypicalStemCellCycleTime()
{
    return 16.0;
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//template class IngeWntSwatCellCycleModelCellsGenerator<1>;
template class IngeWntSwatCellCycleModelCellsGenerator<2>;
//template class IngeWntSwatCellCycleModelCellsGenerator<3>;
