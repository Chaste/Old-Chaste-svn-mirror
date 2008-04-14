/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef NODEMAP_HPP_
#define NODEMAP_HPP_

#include <vector>
#include "Exception.hpp"

class NodeMap
{
private:
    std::vector<unsigned> mMap;
    
public:
    NodeMap(unsigned size)
    {
        // this used to be reserve, but this acts oddly: 
        // eg: mMap.reserve(2); mMap[0]=1;
        // runs and mMap[0] returns 1, but mMap.size() returns 0
        mMap.resize(size); 
    }
    
    void Resize(unsigned size)
    {
        mMap.resize(size);
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
    
    bool IsDeleted(unsigned index)
    {
        return (mMap[index]==UINT_MAX);
    }
    
    unsigned GetNewIndex(unsigned oldIndex) const
    {
        if (mMap[oldIndex] == UINT_MAX)
        {
            EXCEPTION("Node has been deleted");
        }
        return (unsigned) mMap[oldIndex];
    }
    
    bool IsIdentityMap()
    {
        for(unsigned i=0; i<mMap.size(); i++)
        {
            if(mMap[i]!=i)
            {
                return false;
            }
        }
        return true;
    }
    
    unsigned Size()
    {
        return mMap.size();
    }
  
};


#endif /*NODEMAP_HPP_*/
