#ifndef NODEMAP_HPP_
#define NODEMAP_HPP_

#include <vector>
#include "Exception.hpp"

class NodeMap 
{
	private:
	std::vector<int > mMap;
	
	public:
	NodeMap(int size)
	{
		mMap.reserve(size);
	}
	
	void SetNewIndex(unsigned oldIndex, unsigned newIndex)
	{
		mMap[oldIndex] = (int) newIndex;
	}
	
	void SetDeleted(unsigned index)
	{
		mMap[index] = -1;
	}
	
	unsigned GetNewIndex(unsigned oldIndex)
	{
		if (mMap[oldIndex] < 0)
		{
			EXCEPTION("Node has been deleted");
		}
		return (unsigned) mMap[oldIndex];
	}
};


#endif /*NODEMAP_HPP_*/
