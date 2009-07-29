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
#include "NodeBoxCollection.hpp"

/////////////////////////////////////////////////////////////////////////////
// NodeBox methods
/////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
NodeBox<DIM>::NodeBox(c_vector<double, 2*DIM> minAndMaxValues)
{
    mMinAndMaxValues = minAndMaxValues;
}

template<unsigned DIM>
c_vector<double, 2*DIM>& NodeBox<DIM>::rGetMinAndMaxValues()
{
    return mMinAndMaxValues;
}

template<unsigned DIM>
void NodeBox<DIM>::AddNode(Node<DIM>* pNode)
{
    mNodesContained.insert(pNode);
}

template<unsigned DIM>
void NodeBox<DIM>::RemoveNode(Node<DIM>* pNode)
{
    mNodesContained.erase(pNode);
}

template<unsigned DIM>
std::set< Node<DIM>* >& NodeBox<DIM>::rGetNodesContained()
{
    return mNodesContained;
}


/////////////////////////////////////////////////////////////////////////////
// NodeBoxCollection methods
/////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
NodeBoxCollection<DIM>::NodeBoxCollection(double cutOffLength, c_vector<double, 2*DIM> domainSize)
    : mDomainSize(domainSize),
      mCutOffLength(cutOffLength)
{
    assert(DIM==2); /// \todo 3d node box collection
    switch (DIM)
    {
/// commented out as the 1d case is not tested or covered - do we really care about 1d? 
/// if so \todo 1d node box collection
//
//        case 1:
//        {
//            mNumBoxesEachDirection(0) = 0;
//            double box_min_x = domainSize(0);
//            while (box_min_x <=  domainSize(1))
//            {
//                c_vector<double, 2*DIM> box_coords;
//                box_coords(0) = box_min_x;
//                box_coords(1) = box_min_x + cutOffLength;
//
//                NodeBox<DIM> new_box(box_coords);
//                mBoxes.push_back(new_box);
//                mNumBoxesEachDirection(0)++;
//
//                box_min_x += cutOffLength;
//            }
//
//            break;
//        }
        case 2:
        {
            mNumBoxesEachDirection(0) = 0;
            double box_min_x = domainSize(0);
            while (box_min_x <=  domainSize(1))
            {
                double box_min_y = domainSize(2);
                mNumBoxesEachDirection(1) = 0;
                while (box_min_y <= domainSize(3))
                {
                    c_vector<double, 2*DIM> box_coords;
                    box_coords(0) = box_min_x;
                    box_coords(1) = box_min_x + cutOffLength;
                    box_coords(2) = box_min_y;
                    box_coords(3) = box_min_y + cutOffLength;

                    NodeBox<DIM> new_box(box_coords);
                    mBoxes.push_back(new_box);
                    mNumBoxesEachDirection(1)++;

                    box_min_y += cutOffLength;
                }
                mNumBoxesEachDirection(0)++;
                box_min_x += cutOffLength;
            }
            assert(mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1)==mBoxes.size());
            break;
        }
    }
    CalculateLocalBoxes();
}

template<unsigned DIM>
unsigned NodeBoxCollection<DIM>::CalculateContainingBox(Node<DIM>* pNode)
{
    assert(DIM==2);

    double x = pNode->rGetLocation()[0];
    double y = pNode->rGetLocation()[1];

    for (unsigned i=0; i<DIM; i++)
    {
        assert(pNode->rGetLocation()[i] >= mDomainSize(2*i));
        assert(pNode->rGetLocation()[i] <= mDomainSize(2*i+1));
    }
    unsigned box_x_index = (unsigned) floor((x-mDomainSize(0))/mCutOffLength);
    unsigned box_y_index = (unsigned) floor((y-mDomainSize(2))/mCutOffLength);
    assert(mNumBoxesEachDirection(1)*box_x_index + box_y_index < mBoxes.size());
    return mNumBoxesEachDirection(1)*box_x_index + box_y_index;
}

template<unsigned DIM>
NodeBox<DIM>& NodeBoxCollection<DIM>::rGetBox(unsigned boxIndex)
{
    assert(boxIndex < mBoxes.size());
    return mBoxes[boxIndex];
}

template<unsigned DIM>
unsigned NodeBoxCollection<DIM>::GetNumBoxes()
{
    return mBoxes.size();
}

template<unsigned DIM>
void NodeBoxCollection<DIM>::CalculateLocalBoxes()
{
    assert(DIM==2);

    mLocalBoxes.clear();
    for (unsigned box_index=0; box_index<mBoxes.size(); box_index++)
    {
        std::set<unsigned> local_boxes;
        local_boxes.insert(box_index);

        if (!IsBottomRow(box_index))
        {
            local_boxes.insert(box_index-1);
        }

        if (!IsTopRow(box_index))
        {
            local_boxes.insert(box_index+1);
        }

        if (!IsLeftColumn(box_index))
        {
            local_boxes.insert(box_index-mNumBoxesEachDirection(1));
        }

        if (!IsRightColumn(box_index))
        {
            local_boxes.insert(box_index+mNumBoxesEachDirection(1));
        }

        if ( (!IsBottomRow(box_index)) && (!IsLeftColumn(box_index)) )
        {
            local_boxes.insert(box_index-mNumBoxesEachDirection(1)-1);
        }

        if ( (!IsBottomRow(box_index)) && (!IsRightColumn(box_index)) )
        {
            local_boxes.insert(box_index+mNumBoxesEachDirection(1)-1);
        }

        if ( (!IsTopRow(box_index)) && (!IsRightColumn(box_index)) )
        {
            local_boxes.insert(box_index+mNumBoxesEachDirection(1)+1);
        }

        if ( (!IsTopRow(box_index)) && (!IsLeftColumn(box_index)) )
        {
            local_boxes.insert(box_index-mNumBoxesEachDirection(1)+1);
        }

        mLocalBoxes.push_back(local_boxes);
    }
}

template<unsigned DIM>
std::set<unsigned> NodeBoxCollection<DIM>::GetLocalBoxes(unsigned boxIndex)
{
    assert(boxIndex < mLocalBoxes.size());
    return mLocalBoxes[boxIndex];
}

template<unsigned DIM>
void NodeBoxCollection<DIM>::CalculateNodePairs(std::vector<Node<DIM>*>& rNodes, std::set<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs)
{
    rNodePairs.clear();
    for (unsigned node_index=0; node_index<rNodes.size(); node_index++)
    {
        // Get the box containing this node
        unsigned box_index = CalculateContainingBox(rNodes[node_index]);

        // Get the local boxes to this node
        std::set<unsigned> local_boxes_indices = GetLocalBoxes(box_index);

        // Loop over all the local boxes
        for (std::set<unsigned>::iterator iter = local_boxes_indices.begin();
             iter != local_boxes_indices.end();
             iter++)
        {
            NodeBox<DIM>& r_box = mBoxes[*iter];
            std::set< Node<DIM>* >& r_contained_nodes = r_box.rGetNodesContained();

            // Get all the nodes in the local boxes in the original node
            for (typename std::set<Node<DIM>*>::iterator node_iter = r_contained_nodes.begin();
                node_iter != r_contained_nodes.end();
                ++node_iter)
            {
                unsigned index2 = (*node_iter)->GetIndex();

                // If node_index1 < node_index2 add the pair to the set.
                if (node_index < index2)
                {
                    rNodePairs.insert(std::pair<Node<DIM>*,Node<DIM>*>(rNodes[node_index],rNodes[index2]));
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class NodeBox<1>;
template class NodeBox<2>;
template class NodeBox<3>;
template class NodeBoxCollection<1>;
template class NodeBoxCollection<2>;
template class NodeBoxCollection<3>;
