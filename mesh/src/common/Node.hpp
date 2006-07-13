#ifndef _NODE_HPP_
#define _NODE_HPP_

#include "Point.hpp"
#include <set>

template<int SPACE_DIM>
class Node
{
private:
    int mIndex;
    Point<SPACE_DIM> mPoint;
    
    bool mIsBoundaryNode;

    // Set of indices of elements containing this node as a vertex
    std::set<unsigned> mElementIndices;
    std::set<unsigned>::const_iterator mIterator;
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

	//Note setting the point in space is dangerous
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

	const Point<SPACE_DIM> &rGetPoint() const
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
     * @param index of the element to add.
     */
  
    void AddElement(unsigned index)
    {
        mElementIndices.insert(index);
        mIterator=mElementIndices.begin();
    }

    /**
     * Return a set of pointers to elements containing this node as a vertex.
     */
    std::set<unsigned> &rGetContainingElementIndices()
    {
        return mElementIndices;
    }
    
    const int GetNumContainingElements()
    {
        return mElementIndices.size();
    }
    
    
    unsigned GetNextContainingElementIndex()
    {
        unsigned current_containing_element = *mIterator; 
        mIterator++;
        
        if (mIterator == mElementIndices.end())
        {
           mIterator = mElementIndices.begin();
        }
        return current_containing_element;
    }

};


#endif //_NODE_HPP_
