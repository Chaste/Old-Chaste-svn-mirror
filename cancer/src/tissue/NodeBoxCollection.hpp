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


/**
 * A small class for a nD 'box' defined by its min/max x/y/z values which
 * contains a list of nodes located in that box
 */
template<unsigned DIM>
class NodeBox
{
private:

    /** Coordinates of the box, in the form (for 2D) (xmin, xmax, ymin, ymax) (etc). */
    c_vector<double, 2*DIM> mMinAndMaxValues;

    /** Nodes contained in this box. */
    std::set< Node<DIM>* > mNodesContained;

public:

    /**
     * Constructor just takes in the extremal values of the box.
     *
     * @param minAndMaxValues the extremal values
     */
    NodeBox(c_vector<double, 2*DIM> minAndMaxValues);

    /** Get the coordinates of the box, in the form (for 2D) (xmin, xmax, ymin, ymax) (etc). */
    c_vector<double, 2*DIM>& rGetMinAndMaxValues();

    /**
     * Add a node to this box.
     *
     * @param pNode address of the node to be added
     */
    void AddNode(Node<DIM>* pNode);

    /**
     * Remove a node from this box.
     *
     * @param pNode address of the node to be removed
     */
    void RemoveNode(Node<DIM>* pNode);

    /** Get all the nodes in this box. */
    std::set< Node<DIM>* >& rGetNodesContained();
};


/**
 * A collection of 'boxes' partitioning the domain with information on which nodes are located in which box
 *
 * not archived - constructed in NodeBasedTissue constructor.
 */
template<unsigned DIM>
class NodeBoxCollection
{
private:

    /** A vector of boxes to store rough node positions. */
    std::vector< NodeBox<DIM> > mBoxes;

    /** The domain being partitioned. */
    c_vector<double, 2*DIM> mDomainSize;

    /** The width of each box. */
    double mCutOffLength;

    /** Number of boxes in each direction. */
    c_vector<unsigned, DIM> mNumBoxesEachDirection;

    /** The boxes local (itself and nearest neighbour) to a given box. */
    std::vector< std::set<unsigned> > mLocalBoxes;

public:

    /**
     * Constructor
     *
     * @param cutOffLength the width of each box (cutOffLength)
     * @param domainSize the size of the domain, in the form (xmin, xmax, ymin, ymax) (etc)
     */
    NodeBoxCollection(double cutOffLength, c_vector<double, 2*DIM> domainSize);

    /**
     * Calculate which box this node is contained in.
     *
     * @param pNode address of the node
     */
    unsigned CalculateContainingBox(Node<DIM>* pNode);

    /**
     * Get a box.
     *
     * @param boxIndex the index of the box to return
     * @return a NodeBox
     */
    NodeBox<DIM>& rGetBox(unsigned boxIndex);

    /** Get the number of boxes. */
    unsigned GetNumBoxes();

    /** Set up the local boxes (ie itself and its nearest-neighbours) for each of the boxes. */
    void CalculateLocalBoxes();

    /**
     * Get the set of all the local boxes, i.e. itself and its nearest-neighbours.
     *
     * @param boxIndex the index of the box
     */
    std::set<unsigned> GetLocalBoxes(unsigned boxIndex);

    /**
     *  Compute all the pairs of (potentially) connected nodes, ie nodes which are in a local box
     *  to the box containing the first node. **Note that the user still has to check that the node
     *  pairs are less than the cut-off distance apart.** The pairs are checked so that index1 < index2,
     *  so each connected pair of nodes is only in the set once.
     *
     *  @param rNodes all the nodes to be consider
     *  @param rNodePairs the return value, a set of pairs of nodes
     */
    void CalculateNodePairs(std::vector<Node<DIM>*>& rNodes, std::set<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs);
};


#endif /*NODEBOXCOLLECTION_HPP_*/
