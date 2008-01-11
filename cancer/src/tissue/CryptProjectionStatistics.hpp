#ifndef CRYPTPROJECTIONSTATISTICS_HPP_
#define CRYPTPROJECTIONSTATISTICS_HPP_

#include "AbstractCryptStatistics.hpp"
#include "CryptProjectionSpringSystem.hpp"

class CryptProjectionStatistics : public AbstractCryptStatistics
{
protected:

    /**
     * CellIsInSection method.
     * 
     * @param angle  The angle between the crypt section and the x axis in the projection 
     * @param cellPosition  The vector of a cell's position
     * @param widthOfSection The width of the section
     */ 
    bool CellIsInSection(double angle, const c_vector<double,2>& cellPosition, double widthOfSection=0.6);
    
    
public:

    /** 
     *  Constructor
     * 
     *  @param rCrypt  The crypt
     */
    CryptProjectionStatistics(Tissue<2>& rCrypt)
        : AbstractCryptStatistics(rCrypt) 
    {}
    
    /**
     * GetCryptSection method. Takes in an angle from the 
     * interval (-pi, pi].
     * 
     * @param angle  The angle between the crypt section and the x axis in the projection 
     * 
     */ 
    std::vector<TissueCell*> GetCryptSection(double angle = DBL_MAX);
    
     
};

#endif /*CRYPTPROJECTIONSTATISTICS_HPP_*/
