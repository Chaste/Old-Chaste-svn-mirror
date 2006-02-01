#ifndef _NODE_HPP_
#define _NODE_HPP_

#include "Point.hpp"
#include <iostream>
#include <set>

template<int SPACE_DIM>
class Node
{
private:
    int mIndex;
    Point<SPACE_DIM> mPoint;
    
    bool mIsBoundaryNode;

    // Set of pointers to elements containing this node as a vertex
    std::set<const void *> mElements;

	//TODO:
	//bool mIsDirichletNode;
	Node() {}
	    
public:

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
    
    /**
     * Add an element that contains this node.
     * 
     * @param pElement Pointer to the element to add.
     */
    void AddElement(const void *pElement)
    {
        mElements.insert(pElement);
    }

    /**
     * Return a set of pointers to elements containing this node as a vertex.
     */
    std::set<const void *> &rGetContainingElements()
    {
        return mElements;
    }

};


#endif //_NODE_HPP_
