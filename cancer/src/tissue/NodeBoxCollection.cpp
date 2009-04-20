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
void NodeBox<DIM>::AddNode(Node<DIM>* p_node)
{
    mNodesContained.insert(p_node);
}
   
   
template<unsigned DIM>
void NodeBox<DIM>::RemoveNode(Node<DIM>* p_node)
{
    mNodesContained.erase(p_node);
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
    assert(DIM<3); //todo: 3d
    switch (DIM)
    {
        case 1:
        {
            double box_min_x = domainSize(0);
            while (box_min_x <  domainSize(1))
            {
                c_vector<double, 2*DIM> box_coords;
                box_coords(0) = box_min_x;
                box_coords(1) = box_min_x + cutOffLength;

                NodeBox<DIM> new_box(box_coords);
                mBoxes.push_back(new_box);
                mNumBoxesEachDirection(0)++;

                box_min_x += cutOffLength;
            }
            
            break;
        }
        case 2:
        {
            double box_min_x = domainSize(0);
            while (box_min_x <  domainSize(1))
            {
                double box_min_y = domainSize(2);
                mNumBoxesEachDirection(1) = 0;
                while (box_min_y < domainSize(3))
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
//    
//    for(unsigned i=0; i<mBoxes.size(); i++)
//    {
//        std::cout << i << ":  ";
//        
//        std::set<unsigned> local_boxes = GetLocalBoxes(i);
//        for(std::set<unsigned>::iterator iter = local_boxes.begin();
//            iter!= local_boxes.end();
//            iter++)
//        {
//            std::cout << *iter << " ";
//        }
//        std::cout << "\n";
//    }
}

template<unsigned DIM>
unsigned NodeBoxCollection<DIM>::CalculateContainingBox(Node<DIM>* pNode)
{
    assert(DIM==2);

    double x,y;
    x = pNode->rGetLocation()[0];
    y = pNode->rGetLocation()[1];
        
    unsigned box_x_index = (unsigned) floor((x-mDomainSize(0))/mCutOffLength);
    unsigned box_y_index = (unsigned) floor((y-mDomainSize(2))/mCutOffLength);
        
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

///\todo: maybe store these rather than calculate them repeatedly
template<unsigned DIM>
std::set<unsigned> NodeBoxCollection<DIM>::GetLocalBoxes(unsigned boxIndex)
{
    assert(boxIndex < mBoxes.size());
    
    std::set<unsigned> local_boxes;
    local_boxes.insert(boxIndex);

    assert(DIM==2);
    local_boxes.insert(boxIndex);

    if(!IsBottomRow(boxIndex))
    {
        local_boxes.insert(boxIndex-1);
    }
    
    if(!IsTopRow(boxIndex))
    {
        local_boxes.insert(boxIndex+1);
    }

    if(!IsLeftColumn(boxIndex))
    {
        local_boxes.insert(boxIndex-mNumBoxesEachDirection(1));
    }

    if(!IsRightColumn(boxIndex))
    {
        local_boxes.insert(boxIndex+mNumBoxesEachDirection(1));
    }
    
    if( (!IsBottomRow(boxIndex)) && (!IsLeftColumn(boxIndex)) )
    {
        local_boxes.insert(boxIndex-mNumBoxesEachDirection(1)-1);
    }

    if( (!IsBottomRow(boxIndex)) && (!IsRightColumn(boxIndex)) )
    {
        local_boxes.insert(boxIndex+mNumBoxesEachDirection(1)-1);
    }
    
    if( (!IsTopRow(boxIndex)) && (!IsRightColumn(boxIndex)) )
    {
        local_boxes.insert(boxIndex+mNumBoxesEachDirection(1)+1);
    }

    if( (!IsTopRow(boxIndex)) && (!IsLeftColumn(boxIndex)) )
    {
        local_boxes.insert(boxIndex-mNumBoxesEachDirection(1)+1);
    }
    
    return local_boxes;
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
