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
#include "CellPropertyRegistry.hpp"

AbstractCryptStatistics::AbstractCryptStatistics(MeshBasedTissue<2>& rCrypt)
    : mrCrypt(rCrypt)
{}

AbstractCryptStatistics::~AbstractCryptStatistics()
{
}

void AbstractCryptStatistics::LabelSPhaseCells()
{
    for (AbstractTissue<2>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
    	boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        if (cell_iter->GetCellCycleModel()->GetCurrentCellCyclePhase()== S_PHASE)
        {
            // This should only be done for wild type cells (at the moment anyway)
            assert(cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>());

            // Label this cell
            if (!cell_iter->rGetCellPropertyCollection().HasProperty<CellLabel>());
            {
            	cell_iter->AddCellProperty(p_label);
            }
        }
    }
}

void AbstractCryptStatistics::LabelAllCellsAsHealthy()
{
    boost::shared_ptr<AbstractCellMutationState> p_wildtype(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());
    for (AbstractTissue<2>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        cell_iter->SetMutationState(p_wildtype);
        cell_iter->RemoveCellProperty<CellLabel>();
    }
}

std::vector<bool> AbstractCryptStatistics::AreCryptSectionCellsLabelled(std::vector<TissueCellPtr>& rCryptSection)
{
    std::vector<bool> crypt_section_labelled(rCryptSection.size());

    for (unsigned vector_index=0; vector_index<rCryptSection.size(); vector_index++)
    {
        if (rCryptSection[vector_index]->rGetCellPropertyCollection().HasProperty<CellLabel>())
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
