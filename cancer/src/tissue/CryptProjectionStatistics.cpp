#include "CryptProjectionStatistics.hpp"

/** This global function is to allow the list of cells in to be compared in
 *  terms of their y-value and std::list.sort() to be called
 */
bool CellsRadiusComparison(const std::pair<TissueCell*, double> lhs, const std::pair<TissueCell*, double> rhs)
{
    return lhs.second < rhs.second;
}

/**
 * Overridden CellIsInSection method.
 * 
 * @param angle  The angle between the crypt section and the x axis in the projection 
 * @param cellPosition  The vector of a cell's position
 * @param widthOfSection The width of the section
 */ 
bool CryptProjectionStatistics::CellIsInSection(double angle, const c_vector<double,2>& cellPosition, double widthOfSection)
{
    // Get corresponding 3D position of closest point on line
    c_vector<double,2> line_position;
    line_position[0] = norm_2(cellPosition)*cos(angle);
    line_position[1] = norm_2(cellPosition)*sin(angle);
    
    double distance_between_cell_and_line = norm_2(cellPosition - line_position);
    
    return ( distance_between_cell_and_line <= widthOfSection);
}

    
/**
 * Overridden GetCryptSection method. Takes in an angle from the 
 * interval (-pi, pi].
 * 
 * @param angle  The angle between the crypt section and the x axis in the projection 
 * 
 */ 
std::vector<TissueCell*> CryptProjectionStatistics::GetCryptSection(double angle)
{   
    if (angle == DBL_MAX)
    {
        angle = M_PI - 2*M_PI*RandomNumberGenerator::Instance()->ranf();
    }
    
    assert(angle>=-M_PI && angle<=M_PI);
    
    std::list<std::pair<TissueCell*, double> > cells_list; // the second entry is the radius (needed for sorting)
    
    
    // Loop over cells and add to the store if they are within a cell's radius of the
    // specified line  
    for (AbstractTissue<2>::Iterator cell_iter = mrCrypt.Begin();
         cell_iter != mrCrypt.End();
         ++cell_iter)
    {
        
        if( CellIsInSection(angle, cell_iter.rGetLocation()) )
        {
            // Set up a pair, equal to (cell,r) and insert
            std::pair<TissueCell*, double> pair(&(*cell_iter), norm_2(cell_iter.rGetLocation()));
            cells_list.push_back(pair);
        }            
    }

    // Sort the list
    cells_list.sort(CellsRadiusComparison);

    // Copy to a vector
    std::vector<TissueCell*> ordered_cells;
    for(std::list<std::pair<TissueCell*, double> >::iterator iter = cells_list.begin();
        iter!=cells_list.end();
        iter++)
    {
        ordered_cells.push_back(iter->first);
    }

    return ordered_cells;                     
}

