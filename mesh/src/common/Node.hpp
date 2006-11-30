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
    bool mIsDeleted;
    
    // Set of indices of elements containing this node as a vertex
    std::set<unsigned> mElementIndices;
    std::set<unsigned>::const_iterator mElementIterator;
    std::set<unsigned> mBoundaryElementIndices;
    std::set<unsigned>::const_iterator mBoundaryElementIterator;
    //TODO:
    //bool mIsDirichletNode;
    Node()
{}

public:

    Node(int index, Point<SPACE_DIM> point, bool isBoundaryNode=false)
    {
        mPoint = point;
        mIndex = index;
        mIsBoundaryNode = isBoundaryNode;
        mIsDeleted = false;
    }
    
    Node(int index, bool isBoundaryNode=false, double v1=0, double v2=0, double v3=0)
    {
        mPoint = Point<SPACE_DIM>(v1, v2, v3);
        mIndex = index;
        mIsBoundaryNode = isBoundaryNode;
        mIsDeleted = false;
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
        mElementIterator=mElementIndices.begin();
    }
    
    /**
     * Remove an element that contains this node.
     * 
     * @param index of the element to be removed.
     */
    
    void RemoveElement(unsigned index)
    {
        int count = mElementIndices.erase(index);
        if (count == 0)
        {
            EXCEPTION("Tried to remove an index which was not in the set");
        }
        mElementIterator = mElementIndices.begin();
    }
    
    /**
     * Remove an boundary element that contains this node.
     * 
     * @param index of the boundary element to be removed.
     */
    
    void RemoveBoundaryElement(unsigned index)
    {
        int count = mBoundaryElementIndices.erase(index);
        if (count == 0)
        {
            EXCEPTION("Tried to remove an index which was not in the set");
        }
        mBoundaryElementIterator = mBoundaryElementIndices.begin();
    }
    
    /**
     * Add an boundary element that contains this node.
     * 
     * @param index of the element to add.
     */
    
    void AddBoundaryElement(unsigned index)
    {
        mBoundaryElementIndices.insert(index);
        mBoundaryElementIterator=mBoundaryElementIndices.begin();
    }
    
    /**
     * Return a set of pointers to elements containing this node as a vertex.
     */
    std::set<unsigned> &rGetContainingElementIndices()
    {
        return mElementIndices;
    }
    
    /**
     * Return a set of pointers to boundary elements containing this node as a vertex.
     */
    std::set<unsigned> &rGetContainingBoundaryElementIndices()
    {
        return mBoundaryElementIndices;
    }
    
    unsigned GetNumContainingElements()
    {
        return mElementIndices.size();
    }
    
    unsigned GetNumBoundaryElements()
    {
        return mBoundaryElementIndices.size();
    }
    
    
    unsigned GetNextBoundaryElementIndex()
    {
        unsigned current_boundary_element = *mBoundaryElementIterator;
        mBoundaryElementIterator++;
        
        if (mBoundaryElementIterator == mBoundaryElementIndices.end())
        {
            mBoundaryElementIterator = mBoundaryElementIndices.begin();
        }
        return current_boundary_element;
    }
    
    unsigned GetNextContainingElementIndex()
    {
        unsigned current_containing_element = *mElementIterator;
        mElementIterator++;
        
        if (mElementIterator == mElementIndices.end())
        {
            mElementIterator = mElementIndices.begin();
        }
        return current_containing_element;
    }
    
    
    /**
     * Mark a node as having been removed from the mesh
     */
    void MarkAsDeleted()
    {
        mIsDeleted = true;
    }
    
    bool IsDeleted()
    {
        return mIsDeleted;
    }
};


#endif //_NODE_HPP_
