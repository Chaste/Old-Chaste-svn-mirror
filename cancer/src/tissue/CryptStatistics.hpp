#ifndef CRYPTSTATISTICS_HPP_
#define CRYPTSTATISTICS_HPP_

#include "Tissue.cpp"
#include "RandomNumberGenerator.hpp"

/** This global function is to allow the list of cells in to be compared in
 *  terms of their y-value and std::list.sort() to be called
 */
bool CellsHeightComparison(const std::pair<TissueCell*, double> lhs, const std::pair<TissueCell*, double> rhs)
{
    return lhs.second < rhs.second;
}

class CryptStatistics 
{
private:
    Tissue<2>& mrCrypt;

    /**
     *  Method computing the perpendicular distance from the cell to the line from (xBottom,0) to (xTop,yTop), 
     *  and returning if the distance is within the specified width to the section (defaults to 1.0)
     */  
    bool CellIsInSection(double xBottom, double xTop, double yTop, const c_vector<double,2>& cellPosition, double widthOfSection=0.5)
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
    

    /**
     *  Method computing the perpendicular distance from the cell to the line from (xBottom,0) to 
     *  (xTop,yTop), taking into account periodicity, and returning if the distance is within the 
     *  specified width to the section (defaults to 1.0). Done by considering the two possible lines
     *  and checking if cells are within range.
     */  
    bool CellIsInSectionPeriodic(double xBottom, double xTop, double yTop, const c_vector<double,2>& cellPosition, double widthOfSection=1.0)
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
    
    
public :            

    /** 
     *  Constructor
     * 
     *  @param rCrypt The crypt
     */
    CryptStatistics(Tissue<2>& rCrypt)
        : mrCrypt(rCrypt)
    {
    }
    
    
    /**
     *  Get all cells within a cell width of the section defined as the line between points (xBottom,0)
     *  and (xTop,yTop)
     * 
     *  Periodicity can be taken into account (if xTop and xBottom are more than half a crypt 
     *  width apart then a more realistic section will be across the periodic boundary), using the 
     *  final parameter. This obviously requires the mesh to be cylindrical.
     * 
     * @param xBottom    (defaults to a random number U[0,crypt_width])
     * @param xTop  (defaults to a random number U[0,crypt_width])
     * @param yTop  (defaults to crypt_length +2, to get the cells near the top)
     * @param periodic  (defaults to false)
     * 
     * @return  an ordered list of pointes to TissueCells from the bottom to the top of the crypt.
     */
     std::vector<TissueCell*> GetCryptSection(double xBottom = RandomNumberGenerator::Instance()->ranf()*CancerParameters::Instance()->GetCryptWidth(), 
                                             double xTop = RandomNumberGenerator::Instance()->ranf()*CancerParameters::Instance()->GetCryptWidth(), 
                                             double yTop = CancerParameters::Instance()->GetCryptLength() + 2.0, 
                                             bool periodic = false)
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
    
    /** 
     *  Get all cells with a cell width of the line defined by the points (xBottom,0)
     *  and (xTop,yTop), taking into account periodicity
     * 
     *  If xTop and xBottom are more than half a crypt width apart then a more realistic section
     *  will be across the periodic boundary.
     */
    std::vector<TissueCell*> GetCryptSectionPeriodic(double xBottom = RandomNumberGenerator::Instance()->ranf()*CancerParameters::Instance()->GetCryptWidth(), 
                                             double xTop = RandomNumberGenerator::Instance()->ranf()*CancerParameters::Instance()->GetCryptWidth(), 
                                             double yTop = CancerParameters::Instance()->GetCryptLength() + 2.0, 
                                             bool periodic = false)
     {
        return GetCryptSection(xBottom,xTop,yTop,true);
     }
     
    void LabelSPhaseCells()
    {
     
        for (Tissue<2>::Iterator cell_iter = mrCrypt.Begin();
             cell_iter != mrCrypt.End();
             ++cell_iter)
        {
            if ((*cell_iter).GetCellCycleModel()->GetCurrentCellCyclePhase()== S)
            {
                (*cell_iter).SetMutationState(LABELLED);
            }
        } 
     
    }
    
    std::vector<bool> GetWhetherCryptSectionCellsAreLabelled(double xBottom = RandomNumberGenerator::Instance()->ranf()*CancerParameters::Instance()->GetCryptWidth(), 
                                             double xTop = RandomNumberGenerator::Instance()->ranf()*CancerParameters::Instance()->GetCryptWidth(), 
                                             double yTop = CancerParameters::Instance()->GetCryptLength() + 2.0, 
                                             bool periodic = false)
    {
    
        std::vector<bool> crypt_section_labelled ;
        
        //std::vector<TissueCell*> GetCryptSectionPeriodic(xBottom,xTop,yTop,periodic);
        
        return crypt_section_labelled;
    
    }
    
    
};


#endif /*CRYPTSTATISTICS_HPP_*/

