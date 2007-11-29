#ifndef CRYPTSTATISTICS_HPP_
#define CRYPTSTATISTICS_HPP_

#include "Tissue.cpp"
#include "RandomNumberGenerator.hpp"

class CryptStatistics 
{
private:
    Tissue<2>& mrCrypt;

    /**
     *  Method computing the perpendicular distance from the cell to the line from (xBottom,0) to (xTop,yTop), 
     *  and returning if the distance is within the specified width to the section (defaults to 0.5)
     */  
    bool CellIsInSection(double xBottom, double xTop, double yTop, const c_vector<double,2>& cellPosition, double widthOfSection=0.5);

    

    /**
     *  Method computing the perpendicular distance from the cell to the line from (xBottom,0) to 
     *  (xTop,yTop), taking into account periodicity, and returning if the distance is within the 
     *  specified width to the section (defaults to 1.0). Done by considering the two possible lines
     *  and checking if cells are within range.
     */  
    bool CellIsInSectionPeriodic(double xBottom, double xTop, double yTop, const c_vector<double,2>& cellPosition, double widthOfSection=1.0);

    
public :            

    /** 
     *  Constructor
     * 
     *  @param rCrypt The crypt
     */
    CryptStatistics(Tissue<2>& rCrypt)
        : mrCrypt(rCrypt) {};
    
    
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
     * 
     * Note that placing calls to functions with side-effects (eg. changing the random seed)
     * in the default arguments is DANGEROUS.  There is no guarantee that the compiler will 
     * execute these in a sensible order.
     * It appears that Intel goes left-to-right and Gcc goes right-to-left.
     */
     std::vector<TissueCell*> GetCryptSection(double xBottom = DBL_MAX, //RandomNumberGenerator::Instance()->ranf()*CancerParameters::Instance()->GetCryptWidth(), 
                                             double xTop = DBL_MAX, //RandomNumberGenerator::Instance()->ranf()*CancerParameters::Instance()->GetCryptWidth(), 
                                             double yTop = CancerParameters::Instance()->GetCryptLength() + 2.0, 
                                             bool periodic = false);
    
   
    
    /** 
     *  Get all cells with a cell width of the line defined by the points (xBottom,0)
     *  and (xTop,yTop), taking into account periodicity
     * 
     *  If xTop and xBottom are more than half a crypt width apart then a more realistic section
     *  will be across the periodic boundary.
     * Note that placing calls to functions with side-effects (eg. changing the random seed)
     * in the default arguments is DANGEROUS.  There is no guarantee that the compiler will 
     * execute these in a sensible order.
     * It appears that Intel goes left-to-right and Gcc goes right-to-left.
     */
    std::vector<TissueCell*> GetCryptSectionPeriodic(double xBottom = DBL_MAX, //RandomNumberGenerator::Instance()->ranf()*CancerParameters::Instance()->GetCryptWidth(), 
                                             double xTop = DBL_MAX, //RandomNumberGenerator::Instance()->ranf()*CancerParameters::Instance()->GetCryptWidth(), 
                                             double yTop = CancerParameters::Instance()->GetCryptLength() + 2.0);
    
    /**
     * To recreate the Meineke labelling experiments
     * 
     * Cells which are in S phase have their mutation state changed 
     * from 'HEALTHY' to 'LABELLED'.
     */ 
    void LabelSPhaseCells();
    
    /**
     * Sets all the cells in the crypt to have a mutation
     * state of 'HEALTHY'
     */
    void LabelAllCellsAsHealthy();
    
    
    /**
     *  Get all cells within a cell width of the section defined as the line between points (xBottom,0)
     *  and (xTop,yTop). If a patricular cell is labelled then the boolean true is returned.
     * 
     *  Periodicity can be taken into account (if xTop and xBottom are more than half a crypt 
     *  width apart then a more realistic section will be across the periodic boundary), using the 
     *  final parameter. This obviously requires the mesh to be cylindrical.
     * 
     * @param xBottom    (defaults to a random number U[0,crypt_width])
     * @param xTop  (defaults to a random number U[0,crypt_width])
     * @param yTop  (defaults to crypt_length +2, to get the cells near the top)
     * @param periodic  (defaults to true)
     * 
     * @return  a standard vector of booleans which states whether a labelled cell is present at a corresponing position.
     */
    std::vector<bool> GetWhetherCryptSectionCellsAreLabelled(double xBottom = RandomNumberGenerator::Instance()->ranf()*CancerParameters::Instance()->GetCryptWidth(), 
                                             double xTop = RandomNumberGenerator::Instance()->ranf()*CancerParameters::Instance()->GetCryptWidth(), 
                                             double yTop = CancerParameters::Instance()->GetCryptLength() + 2.0, 
                                             bool periodic = true);
    
};


#endif /*CRYPTSTATISTICS_HPP_*/

