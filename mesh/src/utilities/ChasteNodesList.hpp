/*

Copyright (C) University of Oxford, 2005-2010

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


#ifndef CHASTENODESLIST_HPP_
#define CHASTENODESLIST_HPP_

#include "AbstractChasteRegion.hpp"
#include "Node.hpp"
#include "ChastePoint.hpp"

#include <vector>
using namespace std;
/**
 * This class defines a 3D cuboid and provides a method to check
 */
template <unsigned SPACE_DIM>
class ChasteNodesList : public AbstractChasteRegion<SPACE_DIM>
{
private:

    /** A vector to store the list of nodes*/
    std::vector< Node<SPACE_DIM>*> mListOfNodes;

public:

    /**
     * Constructor
     *
     * @param rNodesList a standard vector of (pointer to) nodes
     */
    ChasteNodesList(const std::vector<Node<SPACE_DIM>*> rNodesList) :
        mListOfNodes (rNodesList)
    {
    }


    /**
     * Checks if a given point is contained in the ndoe list.
     *
     * @param rPointToCheck Point to be checked whether it is a node in the list.
     */

    bool DoesContain(const ChastePoint<SPACE_DIM>& rPointToCheck) const
    {
        bool returned_value = false;
        for (unsigned index = 0; index < mListOfNodes.size(); index++)
        {
            if (mListOfNodes[index]->GetPoint().IsSamePoint(rPointToCheck))
            {
                returned_value = true;
                break;
            }
        }

        return returned_value;
    }


};

#endif /*CHASTENODESLIST_HPP_*/
