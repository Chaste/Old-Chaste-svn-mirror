#ifndef CRYPTSTATISTICS_CPP_
#define CRYPTSTATISTICS_CPP_

#include "CryptStatistics.hpp"

/*
 * PRIVATE FUNCTIONS -----------------------------------------------------------------
 */

bool CryptStatistics::CellIsInSection(double xBottom, double xTop, double yTop, const c_vector<double,2>& cellPosition, double widthOfSection)
{
    c_vector<double,2> intercept;

    if(xBottom==xTop)
    {
        intercept[0] = xTop;
        intercept[1] = cellPosition[1];
    }
    else
    {
        double m = (yTop)/(xTop-xBottom); // gradient of line
    
        intercept[0] = (m*m*xBottom + cellPosition[0] + m*cellPosition[1])/(1+m*m);
        intercept[1] = m*(intercept[0] - xBottom);
    }
        
    c_vector<double,2> vec_from_A_to_B = mrCrypt.rGetMesh().GetVectorFromAtoB(intercept, cellPosition);
    double dist = norm_2(vec_from_A_to_B);
   
    return (dist <= widthOfSection);
}

bool CryptStatistics::CellIsInSectionPeriodic(double xBottom, double xTop, double yTop, const c_vector<double,2>& cellPosition, double widthOfSection)
{
    bool is_in_section=false;
    
    c_vector<double,2> intercept;
    double crypt_width = CancerParameters::Instance()->GetCryptWidth();

    double m; // gradient of line
    double offset;
    
    if(xBottom<xTop)
    {
        offset = -crypt_width;    
    }
    else
    {
        offset = crypt_width;
    }
    
    m = (yTop)/(xTop-xBottom+offset); // gradient of line
    
    // 1st Line        
    intercept[0] = (m*m*xBottom + cellPosition[0] + m*cellPosition[1])/(1+m*m);
    intercept[1] = m*(intercept[0] - xBottom);
    
    c_vector<double,2> vec_from_A_to_B = mrCrypt.rGetMesh().GetVectorFromAtoB(intercept, cellPosition);
    double dist = norm_2(vec_from_A_to_B);

    if(dist < widthOfSection)
    {
        is_in_section=true;
    }

    // 2nd Line        
    intercept[0] = (m*m*(xBottom-offset) + cellPosition[0] + m*cellPosition[1])/(1+m*m);
    intercept[1] = m*(intercept[0] - (xBottom-offset));
    
    vec_from_A_to_B = mrCrypt.rGetMesh().GetVectorFromAtoB(intercept, cellPosition);
    dist = norm_2(vec_from_A_to_B);

    if(dist < widthOfSection)
    {
        is_in_section=true;
    }

    return is_in_section;
}
    
/*
 * PUBLIC FUNCTIONS -----------------------------------------------------------------
 */
std::vector<TissueCell*> CryptStatistics::GetCryptSection(double xBottom, double xTop, double yTop, bool periodic)
{
    assert(yTop>0.0);
    std::list<std::pair<TissueCell*, double> > cells_list; // the second entry is the y value (needed for sorting)
    
    if (fabs(xTop-xBottom)<0.5*CancerParameters::Instance()->GetCryptWidth())
    {
        // the periodic version isn't needed, ignore even if periodic was set to true
        periodic = false;
    }    
    
    // loop over cells and add to the store if they are within a cell's radius of the
    // specified line  
    for (Tissue<2>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        if(periodic)
        {
            if(CellIsInSectionPeriodic(xBottom, xTop, yTop, cell_iter.rGetLocation()))
            {
                // set up a pair, equal to (cell,y_val) and insert
                std::pair<TissueCell*, double> pair(&(*cell_iter), cell_iter.rGetLocation()[1]);
                cells_list.push_back(pair);
            }
        }
        else
        {
            if(CellIsInSection(xBottom, xTop, yTop, cell_iter.rGetLocation()))
            {
                // set up a pair, equal to (cell,y_val) and insert
                std::pair<TissueCell*, double> pair(&(*cell_iter), cell_iter.rGetLocation()[1]);
                cells_list.push_back(pair);
            }
        }
    }

    // sort the list
    cells_list.sort(CellsHeightComparison);

    // copy to a vector
    std::vector<TissueCell*> ordered_cells;
    for(std::list<std::pair<TissueCell*, double> >::iterator iter = cells_list.begin();
        iter!=cells_list.end();
        iter++)
    {
        ordered_cells.push_back(iter->first);
    }

    return ordered_cells;
}

std::vector<TissueCell*> CryptStatistics::GetCryptSectionPeriodic(double xBottom, double xTop, double yTop, bool periodic) 
{
   return GetCryptSection(xBottom,xTop,yTop,true);
}   

    
void CryptStatistics::LabelSPhaseCells()
{
    for (Tissue<2>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        if ((*cell_iter).GetCellCycleModel()->GetCurrentCellCyclePhase()== S)
        {
            assert((*cell_iter).GetMutationState() == HEALTHY);
            (*cell_iter).SetMutationState(LABELLED);
        }
    } 
 
}
    
void CryptStatistics::LabelAllCellsAsHealthy()
{
    for (Tissue<2>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        (*cell_iter).SetMutationState(HEALTHY);
    }    
} 

std::vector<bool> CryptStatistics::GetWhetherCryptSectionCellsAreLabelled(double xBottom, 
                                                         double xTop, 
                                                         double yTop, 
                                                         bool periodic)
{
    std::vector<TissueCell*> crypt_section = GetCryptSectionPeriodic(xBottom,xTop,yTop,periodic);
    std::vector<bool> crypt_section_labelled(crypt_section.size()) ;
    
    for (unsigned vector_index=0; vector_index<crypt_section.size(); vector_index++)
    {
        if (crypt_section[vector_index]->GetMutationState() == LABELLED)
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


#endif /*CRYPTSTATISTICS_CPP_*/
