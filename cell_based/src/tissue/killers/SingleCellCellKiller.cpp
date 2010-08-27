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
#include "SingleCellCellKiller.hpp"

template<unsigned DIM>
SingleCellCellKiller<DIM>::SingleCellCellKiller(AbstractTissue<DIM>* pTissue, unsigned number)
: AbstractCellKiller<DIM>(pTissue),
  mNumber(number)
{
}

template<unsigned DIM>
unsigned SingleCellCellKiller<DIM>::GetNumber() const
{
    return mNumber;
}

template<unsigned DIM>
void SingleCellCellKiller<DIM>::TestAndLabelCellsForApoptosisOrDeath()
{
	if (this->mpTissue->GetNumRealCells()==0)
	{
		return;
	}

	typename AbstractTissue<DIM>::Iterator cell_iter = this->mpTissue->Begin();

	for (unsigned i=0; ( (i<mNumber) && (cell_iter!=this->mpTissue->End()) ); i++)
	{
		++cell_iter;
	}

	cell_iter->Kill();
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class SingleCellCellKiller<1>;
template class SingleCellCellKiller<2>;
template class SingleCellCellKiller<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SingleCellCellKiller)
