#ifndef CHASTECUBOID_HPP_
#define CHASTECUBOID_HPP_

#include "ChastePoint.hpp"

/*
 * This class defines a 3D cuboid and provides a method to check 
 * if a given point is contained in the volume.
 */
class ChasteCuboid
{
private:
    ChastePoint<3> mrPointA;
    ChastePoint<3> mrPointB;
    
public:
    /**
     * The cuboid is defined by any two space-diagonal opposite corners.
     * 
     * @param rPointA Any vertex of the cuboid
     * @param rPointB The space-diagonal opposite corner of pointA
     */   
    ChasteCuboid(ChastePoint<3>& rPointA, ChastePoint<3>& rPointB);
        
    /*
     * Checks if a given point is contained in the cuboid.
     * 
     * @param pointToCheck Point to be checked to be contained in.
     */
    bool DoesContain(ChastePoint<3>& rPointToCheck);
};


#endif /*CHASTECUBOID_HPP_*/
