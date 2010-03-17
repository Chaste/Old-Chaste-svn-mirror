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
#include "AbstractCryptStatistics.hpp"
#include "WildTypeCellMutationState.hpp"
#include "LabelledCellMutationState.hpp"

void AbstractCryptStatistics::LabelSPhaseCells()
{
    for (AbstractTissue<2>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        if (cell_iter->GetCellCycleModel()->GetCurrentCellCyclePhase()== S_PHASE)
        {
            // This should only be done for healthy or labelled populations, not mutants (at the moment anyway)
            assert( cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>()
            		|| cell_iter->GetMutationState()->IsType<LabelledCellMutationState>() );
            boost::shared_ptr<AbstractCellMutationState> p_labelled(new LabelledCellMutationState);
            cell_iter->SetMutationState(p_labelled);
        }
    }
}

void AbstractCryptStatistics::LabelAllCellsAsHealthy()
{
	boost::shared_ptr<AbstractCellMutationState> p_wildtype(new WildTypeCellMutationState);
    for (AbstractTissue<2>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        cell_iter->SetMutationState(p_wildtype);
    }
}

std::vector<bool> AbstractCryptStatistics::GetWhetherCryptSectionCellsAreLabelled(std::vector<TissueCell*> cryptSection)
{
    std::vector<bool> crypt_section_labelled(cryptSection.size());

    for (unsigned vector_index=0; vector_index<cryptSection.size(); vector_index++)
    {
        if (cryptSection[vector_index]->GetMutationState()->IsType<LabelledCellMutationState>())
        {
            crypt_section_labelled[vector_index] = true;
        }
        else
        {
            crypt_section_labelled[vector_index] = false;
        }
    }

    return crypt_section_labelled;
}
