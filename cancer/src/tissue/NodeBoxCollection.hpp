#ifndef NODEBOXCOLLECTION_HPP_
#define NODEBOXCOLLECTION_HPP_

#include "Node.hpp"

//#include <boost/serialization/set.hpp>
//#include <boost/serialization/vector.hpp>

/** A small class for a nD 'box' defined by its min/max x/y/z values which
 *  contains a list of nodes located in that box
 */ 
template<unsigned DIM>
class NodeBox
{
private:
    /** Coordinates of the box, in the form (for 2D) (xmin, xmax, ymin, ymax) (etc) */
    c_vector<double, 2*DIM> mMinAndMaxValues;
    /** Nodes contained in this box */
    std::set< Node<DIM>* > mNodesContained;
        
public:
    /** Constructor just takes in the extremal values of the box */
    NodeBox(c_vector<double, 2*DIM> minAndMaxValues);
    
    /** Get the coordinates of the box, in the form (for 2D) (xmin, xmax, ymin, ymax) (etc) */
    c_vector<double, 2*DIM>& rGetMinAndMaxValues();
    
    /** Add a node to this box */
    void AddNode(Node<DIM>* p_node);
    /** Remove a node from this box */
    void RemoveNode(Node<DIM>* p_node);
    
    /** Get all the nodes in this box */
    std::set< Node<DIM>* >& rGetNodesContained();
};

/**
 *  A collection of 'boxes' partitioning the domain with information on which nodes are located in which box
 */
template<unsigned DIM>
class NodeBoxCollection 
{
private:
    /** A vector of boxes to store rough node positions */
    std::vector< NodeBox<DIM> > mBoxes;
    
    /** The domain being partitioned */
    c_vector<double,2*DIM> mDomainSize;

    /** The width of each box */
    double mCutOffLength;
    
    /** Number of boxes in each direction */
    c_vector<unsigned,DIM> mNumBoxesEachDirection;
    
    ////////////////////////
    // helper methods
    ////////////////////////
    /** Whether a box is on the bottom row of all the boxes */
    bool IsBottomRow(unsigned boxIndex)
    {
        return boxIndex % mNumBoxesEachDirection(1)==0;
    }

    /** Whether a box is on the top row of all the boxes */
    bool IsTopRow(unsigned boxIndex)
    {
        return boxIndex % mNumBoxesEachDirection(1)==mNumBoxesEachDirection(1)-1;
    }

    /** Whether a box is on the left side of the boxes */
    bool IsLeftColumn(unsigned boxIndex)
    {
        return boxIndex < mNumBoxesEachDirection(1);
    }
    
    /** Whether a box is on the right side of the boxes */
    bool IsRightColumn(unsigned boxIndex)
    {
        return boxIndex >= mBoxes.size() - mNumBoxesEachDirection(1);
    }

    
public:
    /** Constructor takes in the width of each box (cutOffLength) and the size of the domain, in the form
     *  (xmin, xmax, ymin, ymax) (etc)
     */ 
    NodeBoxCollection(double cutOffLength, c_vector<double, 2*DIM> domainSize);

    /** Calculate which box this node is contained in */
    unsigned CalculateContainingBox(Node<DIM>* pNode);

    /** Get a box */
    NodeBox<DIM>& rGetBox(unsigned boxIndex);

    /** Get the number of boxes */
    unsigned GetNumBoxes();
    
    /** Returns a set of all the local boxes, ie itself and its nearest-neighbours */
    std::set<unsigned> GetLocalBoxes(unsigned boxIndex);
};
    

#endif /*NODEBOXCOLLECTION_HPP_*/
