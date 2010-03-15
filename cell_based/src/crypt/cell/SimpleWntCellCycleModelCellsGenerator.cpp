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
#include "SimpleWntCellCycleModelCellsGenerator.hpp"


template<unsigned DIM>
AbstractCellCycleModel* SimpleWntCellCycleModelCellsGenerator<DIM>::CreateCellCycleModel()
{
//	SimpleWntCellCycleModel p_cell_cycle_model = new SimpleWntCellCycleModel();
//	p_cell_cycle_model->SetDimension(DIM);
	return new SimpleWntCellCycleModel(DIM);
}


template<unsigned DIM>
double SimpleWntCellCycleModelCellsGenerator<DIM>::GetTypicalTransitCellCycleTime()
{
    return TissueConfig::Instance()->GetTransitCellG1Duration()
            + TissueConfig::Instance()->GetSG2MDuration();
}


template<unsigned DIM>
double SimpleWntCellCycleModelCellsGenerator<DIM>::GetTypicalStemCellCycleTime()
{
    return TissueConfig::Instance()->GetStemCellG1Duration()
            + TissueConfig::Instance()->GetSG2MDuration();
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//template class SimpleWntCellCycleModelCellsGenerator<1>;
template class SimpleWntCellCycleModelCellsGenerator<2>;
//template class SimpleWntCellCycleModelCellsGenerator<3>;
