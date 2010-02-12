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

#include "DistanceMapCalculator.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::EuclideanDistanceTwoPoints(
        const c_vector<double, SPACE_DIM>& pointA,
        const c_vector<double, SPACE_DIM>& pointB) const
{
    return norm_2(pointA-pointB);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::CartToEucliDistance(
        c_vector<double, SPACE_DIM>& cartDistance) const
{
    return norm_2(cartDistance);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::DistanceMapCalculator(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
    : mrMesh(rMesh)
{
    mNumNodes = mrMesh.GetNumNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::ComputeDistanceMap(
        const std::vector<unsigned>& rOriginSurface,
        std::vector<double>& rNodeDistances)
{
    rNodeDistances.resize(mNumNodes);
     /*
     * Matrix of distances along each dimension (initialised to +inf)
     */
    std::vector< c_vector<double, SPACE_DIM> >  cart_distances(mNumNodes);
    for (unsigned index=0; index<mNumNodes; index++)
    {
        for (unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            cart_distances[index][dim] = DBL_MAX;
        }
    }

    /*
     * Queue of nodes to be processed (initialised with the nodes defining the surface)
     */
    std::queue<unsigned> active_nodes;
    for (unsigned index=0; index<rOriginSurface.size(); index++)
    {
        active_nodes.push(rOriginSurface[index]);

        for (unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            cart_distances[rOriginSurface[index]][dim] = 0.0;
        }
    }

    WorkOnLocalQueue(active_nodes, cart_distances);
    
    //Finish up
    for (unsigned index=0; index<mNumNodes; index++)
    {
        rNodeDistances[index] = CartToEucliDistance(cart_distances[index]);
    }

}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::WorkOnLocalQueue(std::queue<unsigned>& activeNodeIndexQueue, std::vector< c_vector<double, SPACE_DIM> >& cartDistances )
{
   while (!activeNodeIndexQueue.empty())
    {
        // Get a pointer to the next node in the queue
        unsigned current_node_index = activeNodeIndexQueue.front();
        Node<SPACE_DIM>* p_current_node = mrMesh.GetNode(current_node_index);

        // Loop over the elements containing the given node
        for(typename Node<SPACE_DIM>::ContainingElementIterator element_iterator = p_current_node->ContainingElementsBegin();
            element_iterator != p_current_node->ContainingElementsEnd();
            ++element_iterator)
        {
            // Get a pointer to the container element
            Element<ELEMENT_DIM, SPACE_DIM>* p_containing_element = mrMesh.GetElement(*element_iterator);

           // Loop over the nodes of the element
           for(unsigned node_local_index=0;
               node_local_index<p_containing_element->GetNumNodes();
               node_local_index++)
           {
                Node<SPACE_DIM>* p_neighbour_node = p_containing_element->GetNode(node_local_index);
                unsigned neighbour_node_index = p_neighbour_node->GetIndex();

                // Avoid revisiting the active node
                if(neighbour_node_index != current_node_index)
                {
                    // Test if we have found a shorter path from the origin surface to the current neighbour through current node
                    if ((CartToEucliDistance(cartDistances[current_node_index])
                          + EuclideanDistanceTwoPoints(p_current_node->rGetLocation(), p_neighbour_node->rGetLocation()))
                        < CartToEucliDistance(cartDistances[neighbour_node_index])*(1.0 - DBL_EPSILON) )
                    {
                        cartDistances[neighbour_node_index] = cartDistances[current_node_index]
                                                               + (p_current_node->rGetLocation() - p_neighbour_node->rGetLocation());
                        activeNodeIndexQueue.push(neighbour_node_index);
                    }
                }
           }
        }

        activeNodeIndexQueue.pop();
    }
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class DistanceMapCalculator<1, 1>;
template class DistanceMapCalculator<1, 2>;
template class DistanceMapCalculator<2, 2>;
template class DistanceMapCalculator<1, 3>;
//template class DistanceMapCalculator<2, 3>;
template class DistanceMapCalculator<3, 3>;
