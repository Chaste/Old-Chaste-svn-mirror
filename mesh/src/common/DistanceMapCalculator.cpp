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
#include "ParallelTetrahedralMesh.hpp" //For dynamic cast
//#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::DistanceMapCalculator(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
    : mrMesh(rMesh),
      mWorkOnEntireMesh(true),
      mNumHalosPerProcess(NULL)
{
    mNumNodes = mrMesh.GetNumNodes();

    ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* p_parallel_mesh = dynamic_cast<ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>*>(&mrMesh);
    if ( PetscTools::IsSequential() || p_parallel_mesh == NULL)
    {
        //It's a non-distributed mesh...
        mLo=0;
        mHi=mNumNodes;
    }
    else
    {
        //It's a parallel (distributed) mesh (p_parallel_mesh is non-null and we are running in parallel)
        mWorkOnEntireMesh=false;
        mLo = mrMesh.GetDistributedVectorFactory()->GetLow();
        mHi = mrMesh.GetDistributedVectorFactory()->GetHigh();
        //Get local halo information
        p_parallel_mesh->GetHaloNodeIndices(mHaloNodeIndices);
        //Share information on the number of halo nodes
        unsigned my_size=mHaloNodeIndices.size();
        mNumHalosPerProcess=new unsigned[PetscTools::GetNumProcs()];
        MPI_Alltoall(&my_size, 1, MPI_UNSIGNED, 
                     mNumHalosPerProcess, 1, MPI_UNSIGNED, PETSC_COMM_WORLD);
    }
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
        rNodeDistances[index] = DBL_MAX;
    }

    for (unsigned index=0; index<rOriginSurface.size(); index++)
    {
        PushLocal(rOriginSurface[index]);
        rNodeDistances[rOriginSurface[index]] = 0.0;

        for (unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            cart_distances[rOriginSurface[index]][dim] = 0.0;
        }
    }

    WorkOnLocalQueue(cart_distances, rNodeDistances);

    if (mWorkOnEntireMesh==false)
    {
        //Work out how to take non-local updates
        std::vector<double> best_distances(mNumNodes);
        MPI_Allreduce(&rNodeDistances[0], &best_distances[0], mNumNodes, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
        for (unsigned index=0;index<mNumNodes;index++)
        {
            if (best_distances[index] < rNodeDistances[index])
            {
                PushLocal(index);
                rNodeDistances[index]=best_distances[index];
            }
        }
        if (!mActiveNodeIndexQueue.empty())
        {
            std::cout<<mActiveNodeIndexQueue.size()<<"\n";
        }
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::WorkOnLocalQueue(std::vector< c_vector<double, SPACE_DIM> >& cartDistances,
                                                                     std::vector<double>& rNodeDistances)
{
   while (!mActiveNodeIndexQueue.empty())
    {
        // Get a pointer to the next node in the queue
        unsigned current_node_index = mActiveNodeIndexQueue.front();
        //PRINT_VARIABLES(activeNodeIndexQueue.size(), current_node_index);
        mActiveNodeIndexQueue.pop();
        
        try 
        {
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
                    ///\todo - This line is dangerous if it throws
                    Node<SPACE_DIM>* p_neighbour_node = p_containing_element->GetNode(node_local_index);
                    unsigned neighbour_node_index = p_neighbour_node->GetIndex();
    
                    // Avoid revisiting the active node
                    if(neighbour_node_index != current_node_index)
                    {
                        // Test if we have found a shorter path from the origin surface to the current neighbour through current node
                        if ((rNodeDistances[current_node_index]
                            + norm_2(p_current_node->rGetLocation() - p_neighbour_node->rGetLocation()))
                            < rNodeDistances[neighbour_node_index]*(1.0 - DBL_EPSILON) )
                        {
                            cartDistances[neighbour_node_index] = cartDistances[current_node_index]
                                                                   + (p_current_node->rGetLocation() - p_neighbour_node->rGetLocation());
                            //This will save some sqrts later...
                            rNodeDistances[neighbour_node_index] = norm_2(cartDistances[neighbour_node_index]);                                      
                            PushLocal(neighbour_node_index);
                        }
                    }
               }//Node
           }//Element
        }//Try
        catch (Exception &e)
        {
            //Node in the queue doesn't belong to process
            NEVER_REACHED;
        }

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
