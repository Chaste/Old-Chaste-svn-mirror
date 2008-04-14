#include "AbstractCryptStatistics.hpp"

void AbstractCryptStatistics::LabelSPhaseCells()
{
    for (AbstractTissue<2>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        if ((*cell_iter).GetCellCycleModel()->GetCurrentCellCyclePhase()== S_PHASE)
        {   // This should only be done for healthy or labelled populations, not mutants (at the moment anyway)
            assert((*cell_iter).GetMutationState() == HEALTHY || (*cell_iter).GetMutationState() == LABELLED);
            (*cell_iter).SetMutationState(LABELLED);
        }
    } 
}

void AbstractCryptStatistics::LabelAllCellsAsHealthy()
{
    for (AbstractTissue<2>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        (*cell_iter).SetMutationState(HEALTHY);
    }    
}

std::vector<bool> AbstractCryptStatistics::GetWhetherCryptSectionCellsAreLabelled(std::vector<TissueCell*> cryptSection)
{
    std::vector<bool> crypt_section_labelled(cryptSection.size()) ;
    
    for (unsigned vector_index=0; vector_index<cryptSection.size(); vector_index++)
    {
        if (cryptSection[vector_index]->GetMutationState() == LABELLED)
        {
            crypt_section_labelled[vector_index]=true;
        }
        else
        {   
            crypt_section_labelled[vector_index]=false;
        }
    }
            
    return crypt_section_labelled;
}
