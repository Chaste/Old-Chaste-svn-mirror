#ifndef _NODE_HPP_
#define _NODE_HPP_

#include "Point.hpp"
#include <iostream>

template<int SPACE_DIM>
class Node
{
private:
    int mIndex;
    Point<SPACE_DIM> mPoint;
    
    bool mIsBoundaryNode;

	//TODO:
	//bool mIsDirichletNode;
    
public:

	Node() {}
	
	Node(int index, Point<SPACE_DIM> point, bool isBoundaryNode=false)
	{
		mPoint = point;
		mIndex = index;
		mIsBoundaryNode = isBoundaryNode;
	}
	
	Node(int index, bool isBoundaryNode=false, double v1=0, double v2=0, double v3=0)
	{
		mPoint = Point<SPACE_DIM>(v1, v2, v3);
		mIndex = index;
		mIsBoundaryNode = isBoundaryNode;
	}

	void SetPoint(Point<SPACE_DIM> point)
	{
		mPoint = point;
	}
	
	/**
	 * This method should only be called during mesh generation.
	 */
	void SetIndex(int index)
	{
		mIndex = index;
	}
	
	void SetAsBoundaryNode(bool value=true)
	{
		mIsBoundaryNode = value;
	}
			
	Point<SPACE_DIM> GetPoint() const
	{
		return mPoint;
	}
	
	int GetIndex() const
	{
		return mIndex;
	}

	bool IsBoundaryNode() const
	{
		return mIsBoundaryNode;
	}
};


#endif //_NODE_HPP_
