/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
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
