#ifndef NODEMAP_HPP_
#define NODEMAP_HPP_

#include <vector>
#include "Exception.hpp"

class NodeMap
{
private:
    std::vector<unsigned > mMap;
    
public:
    NodeMap(unsigned size)
    {
        mMap.reserve(size);
    }
    
    void Reserve(unsigned size)
    {
        mMap.reserve(size);
    }
    
    void ResetToIdentity()
    {
        for (unsigned oldIndex=0; oldIndex<mMap.size(); oldIndex++)
        {
            mMap[oldIndex] = oldIndex;
        }
    }
    void SetNewIndex(unsigned oldIndex, unsigned newIndex)
    {
        mMap[oldIndex] = newIndex;
    }
    
    void SetDeleted(unsigned index)
    {
        mMap[index] = UINT_MAX;
    }
    
    unsigned GetNewIndex(unsigned oldIndex) const
    {
        if (mMap[oldIndex] == UINT_MAX)
        {
            EXCEPTION("Node has been deleted");
        }
        return (unsigned) mMap[oldIndex];
    }
};


#endif /*NODEMAP_HPP_*/
