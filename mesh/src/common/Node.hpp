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
    c_vector<double, SPACE_DIM> mLocation;
    
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
        mLocation = point.rGetLocation();
        mPoint = point;
        mIndex = index;
        mIsBoundaryNode = isBoundaryNode;
        mIsDeleted = false;
    }
    
    Node(int index, bool isBoundaryNode=false, double v1=0, double v2=0, double v3=0)
    {
        mLocation[0] = v1;
        if (SPACE_DIM > 1)
        {
            mLocation[1] = v2;
            if (SPACE_DIM > 2)
            {
                mLocation[2] = v3;
            }
        }
        mPoint.rGetLocation() = mLocation;
        
        mIndex = index;
        mIsBoundaryNode = isBoundaryNode;
        mIsDeleted = false;
    }
    
    /**
     * Note setting the point in space is dangerous
     * Jacobian and JacobianDeterminant of element need to be updated
     */
    void SetPoint(Point<SPACE_DIM> point)
    {
        mLocation = point.rGetLocation();
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
        return Point<SPACE_DIM>(mLocation);
    }
    
    /**
     * DON'T UPDATE THE POINT YOU GET!
     */
    const Point<SPACE_DIM> &rGetPoint() const
    {
        return mPoint;
    }
    
    /**
     * The returned location may not be modified; if you want that functionality use
     * rGetModifiableLocation instead.
     */
    const c_vector<double, SPACE_DIM> &rGetLocation() const
    {
        return mLocation;
    }
    
    /**
     * If you modify the returned location,
     * Jacobian and JacobianDeterminant of element need to be updated.
     */
    c_vector<double, SPACE_DIM> &rGetModifiableLocation()
    {
        return mLocation;
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
