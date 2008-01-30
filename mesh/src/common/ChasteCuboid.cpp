#include "ChasteCuboid.hpp"

ChasteCuboid::ChasteCuboid(ChastePoint<3>& rPointA, ChastePoint<3>& rPointB): mrPointA(rPointA), mrPointB(rPointB)
{
}

bool ChasteCuboid::DoesContain(const ChastePoint<3>& rPointToCheck)
{
    for (unsigned dim=0; dim<3; dim++)
    {
        if (rPointToCheck[dim] != mrPointA[dim])
        {
            if (rPointToCheck[dim] > mrPointA[dim])
            {
                if (rPointToCheck[dim] > mrPointB[dim])
                {
                    return false;
                }
            }
            else
            {
                if (rPointToCheck[dim] < mrPointB[dim])
                {
                    return false;  
                }
            }
        }
    }
                    
    return true;
}
