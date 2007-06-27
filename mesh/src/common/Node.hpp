#ifndef _NODE_HPP_
#define _NODE_HPP_

//#include "ConformingTetrahedralMesh.hpp"
#include "Point.hpp"
#include <set>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ConformingTetrahedralMesh;

template<unsigned SPACE_DIM>
class Node
{
private:
    unsigned mIndex;
    
    c_vector<double, SPACE_DIM> mLocation;
    
    bool mIsBoundaryNode;
    bool mIsDeleted;
    
    // Set of indices of elements containing this node as a vertex
    std::set<unsigned> mElementIndices;
    std::set<unsigned> mBoundaryElementIndices;

    /**
     * Extraction of commonality between the constructors
     */
    void CommonConstructor(unsigned index, bool isBoundaryNode)
    {
        mIndex = index;
        mIsBoundaryNode = isBoundaryNode;
        mIsDeleted = false;
        mElementIterator = ContainingElementsBegin();
        mBoundaryElementIterator = ContainingBoundaryElementsBegin();
    }

public:
    ~Node()
    {
    }

    Node(unsigned index, Point<SPACE_DIM> point, bool isBoundaryNode=false)
    {
        mLocation = point.rGetLocation();
        CommonConstructor(index, isBoundaryNode);
    }
    
    Node(unsigned index, std::vector<double> coords, bool isBoundaryNode=false)
    {
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            mLocation(i) = coords.at(i);
        }
        CommonConstructor(index, isBoundaryNode);
    }
    
    Node(unsigned index, c_vector<double, SPACE_DIM> location, bool isBoundaryNode=false)
    {
        mLocation = location;
        CommonConstructor(index, isBoundaryNode);
    }
    
    Node(unsigned index, bool isBoundaryNode=false, double v1=0, double v2=0, double v3=0)
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
        CommonConstructor(index, isBoundaryNode);
    }
    
    /**
     * Note setting the point in space is dangerous
     * Jacobian and JacobianDeterminant of element need to be updated
     */
    void SetPoint(Point<SPACE_DIM> point)
    {
        mLocation = point.rGetLocation();
    }
    
    /**
     * This method should only be called during mesh generation.
     */
    void SetIndex(unsigned index)
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
     * The returned location may not be modified; if you want that functionality use
     * rGetModifiableLocation instead.
     */
    const c_vector<double, SPACE_DIM>& rGetLocation() const
    {
        assert(!mIsDeleted);
        return mLocation;
    }
    
    /**
     * If you modify the returned location,
     * Jacobian and JacobianDeterminant of elements need to be updated.
     *
     * Don't forget to assign the result of this call to a reference!
     */
    c_vector<double, SPACE_DIM> &rGetModifiableLocation()
    {
        assert(!mIsDeleted);
        return mLocation;
    }
    
    unsigned GetIndex() const
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
        mElementIterator = mElementIndices.begin();
    }
    
    /**
     * Remove an element that contains this node.
     * 
     * @param index of the element to be removed.
     */
    
    void RemoveElement(unsigned index)
    {
        unsigned count = mElementIndices.erase(index);
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
        unsigned count = mBoundaryElementIndices.erase(index);
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
    
    /**
     * Deprecated iteration method; use the ContainingBoundaryElementIterator instead.
     */
    unsigned GetNextBoundaryElementIndex()
    {
        unsigned current_boundary_element = *mBoundaryElementIterator;
        ++mBoundaryElementIterator;
        
        if (!(mBoundaryElementIterator != ContainingBoundaryElementsEnd()))
        {
            mBoundaryElementIterator = ContainingBoundaryElementsBegin();
        }
        return current_boundary_element;
    }
    
    /**
     * Deprecated iteration method; use the ContainingElementIterator instead.
     */
    unsigned GetNextContainingElementIndex()
    {
        unsigned current_containing_element = *mElementIterator;
        ++mElementIterator;
        
        if (!(mElementIterator != ContainingElementsEnd()))
        {
            mElementIterator = ContainingElementsBegin();
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
    
    /**
     * Determine if a node lives within a flagged element.
     */
    template <unsigned ELEMENT_DIM>
    bool IsFlagged(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
    {
        bool in_flagged_element = false;
        for (ContainingElementIterator it = ContainingElementsBegin();
             it != ContainingElementsEnd();
             ++it)
        {
            if (rMesh.GetElement(*it)->IsFlagged())
            {
                in_flagged_element = true;
                break;
            }
        }
        return in_flagged_element;
    }
    
    /**
     * An iterator over the indices of elements which contain this node.
     */
    class ContainingElementIterator
    {
    public:
        ContainingElementIterator(std::set<unsigned>::const_iterator indexIterator)
            : mIndexIterator(indexIterator)
        {}
        
        /**
         * A default constructor allows users to declare an iterator without assigning to it
         */
        ContainingElementIterator()
        {}
        
        const unsigned& operator*() const
        {
            return *mIndexIterator;
        }
        
        bool operator!=(const ContainingElementIterator& other) const
        {
            return mIndexIterator != other.mIndexIterator;
        }
        
        ContainingElementIterator& operator++()
        {
            ++mIndexIterator;
            return *this;
        }
    private:
        std::set<unsigned>::const_iterator mIndexIterator;
    };
    
    ContainingElementIterator ContainingElementsBegin() const
    {
        return ContainingElementIterator(mElementIndices.begin());
    }
    
    ContainingElementIterator ContainingElementsEnd() const
    {
        return ContainingElementIterator(mElementIndices.end());
    }
    
    /**
     * An iterator over the indices of boundary elements which contain this node.
     */
    class ContainingBoundaryElementIterator
    {
    public:
        ContainingBoundaryElementIterator(std::set<unsigned>::const_iterator indexIterator)
            : mIndexIterator(indexIterator)
        {}
        
        /**
         * A default constructor allows users to declare an iterator without assigning to it
         */
        ContainingBoundaryElementIterator()
        {}
        
        const unsigned& operator*() const
        {
            return *mIndexIterator;
        }
        
        bool operator!=(const ContainingBoundaryElementIterator& other) const
        {
            return mIndexIterator != other.mIndexIterator;
        }
        
        ContainingBoundaryElementIterator& operator++()
        {
            ++mIndexIterator;
            return *this;
        }
    private:
        std::set<unsigned>::const_iterator mIndexIterator;
    };
    
    ContainingBoundaryElementIterator ContainingBoundaryElementsBegin() const
    {
        return ContainingBoundaryElementIterator(mBoundaryElementIndices.begin());
    }
    
    ContainingBoundaryElementIterator ContainingBoundaryElementsEnd() const
    {
        return ContainingBoundaryElementIterator(mBoundaryElementIndices.end());
    }

private:
    typename Node<SPACE_DIM>::ContainingElementIterator mElementIterator;
    typename Node<SPACE_DIM>::ContainingBoundaryElementIterator mBoundaryElementIterator;

};


#endif //_NODE_HPP_
