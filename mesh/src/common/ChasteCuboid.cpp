#include "ChasteCuboid.hpp"

ChasteCuboid::ChasteCuboid(ChastePoint<3>& rPointA, ChastePoint<3>& rPointB) : mrPointA(rPointA), mrPointB(rPointB)
{
    for (unsigned dim=0; dim<3; dim++)
    {
        if (mrPointA[dim] > mrPointB[dim])
        {
            EXCEPTION("Attempt to create a cuboid with MinCorner greater than MaxCorner in some dimension");
        }
    }
}

bool ChasteCuboid::DoesContain(const ChastePoint<3>& rPointToCheck)
{
    bool inside=true;
    for (unsigned dim=0; dim<3; dim++)
    {
        if (rPointToCheck[dim] < mrPointA[dim] - 100*DBL_EPSILON 
            || mrPointB[dim] + 100* DBL_EPSILON < rPointToCheck[dim])
        {
            inside=false;
            break;
        }
    }
    return inside;
}

