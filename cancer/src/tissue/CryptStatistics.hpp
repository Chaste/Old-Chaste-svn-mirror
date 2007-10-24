#ifndef CRYPTSTATISTICS_HPP_
#define CRYPTSTATISTICS_HPP_

#include "Tissue.cpp"

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
    bool CellIsInSection(double xBottom, double xTop, double yTop, const c_vector<double,2>& cellPosition, double widthOfSection=1.0)
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
   
        return (dist < widthOfSection);
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
    
    
    std::vector<TissueCell*> GetCryptSection(double xBottom, double xTop, double yTop)
    {
        assert(yTop>0.0);
        std::list<std::pair<TissueCell*, double> > cells_list; // the second entry is the y value (needed for sorting)
        
        // loop over cells and add to the store if they are within a cell's radius of the
        // specified line  
        for (Tissue<2>::Iterator cell_iter = mrCrypt.Begin();
             cell_iter != mrCrypt.End();
             ++cell_iter)
        {
            if(CellIsInSection(xBottom, xTop, yTop, cell_iter.rGetLocation()))
            {
                // set up a pair, equal to (cell,y_val) and insert
                std::pair<TissueCell*, double> pair(&(*cell_iter), cell_iter.rGetLocation()[1]);
                cells_list.push_back(pair);
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
};


#endif /*CRYPTSTATISTICS_HPP_*/

