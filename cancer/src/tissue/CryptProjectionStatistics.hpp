#ifndef CRYPTPROJECTIONSTATISTICS_HPP_
#define CRYPTPROJECTIONSTATISTICS_HPP_

#include "CryptStatistics.hpp"
#include "CryptProjectionSpringSystem.hpp"

class CryptProjectionStatistics : public CryptStatistics
{
protected:
    /**
     * Overridden CellIsInSection method.
     * 
     * @param angle  The angle between the crypt section and the x axis in the projection 
     * @param cellPosition  The vector of a cell's position
     * @param widthOfSection The width of the section
     */ 
    bool CellIsInSection(double angle, const c_vector<double,2>& cellPosition, double widthOfSection=0.5);
    
    /**
     * Because the above method overloads the method name of the base class as well as 
     * overriding it we need a straight overriding method (for the Intel compiler)  
     */  
    bool CellIsInSection(double xBottom, double xTop, double yTop, const c_vector<double,2>& cellPosition, double widthOfSection=0.5)
    {
        assert(0);
        return false;
    }
    
public:

    /** 
     *  Constructor
     * 
     *  @param rCrypt  The crypt
     */
    CryptProjectionStatistics(Tissue<2>& rCrypt)
        : CryptStatistics(rCrypt) 
    {}
    
    /**
     * Overridden GetCryptSection method. Takes in an angle from the 
     * interval (-pi, pi].
     * 
     * @param angle  The angle between the crypt section and the x axis in the projection 
     * 
     */ 
    std::vector<TissueCell*> GetCryptSection(double angle = DBL_MAX,
                                             double unused1 = 0.0,
                                             double unused2 = 0.0,
                                             bool unused3 = false);
    
     
};

#endif /*CRYPTPROJECTIONSTATISTICS_HPP_*/
