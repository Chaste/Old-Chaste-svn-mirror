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
#include "Debug.hpp"

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
        MPI_Allgather(&my_size, 1, MPI_UNSIGNED, 
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

    bool non_empty_queue=true;
    unsigned round_counter=0;
    while (non_empty_queue)
    {
        WorkOnLocalQueue(cart_distances, rNodeDistances);
        non_empty_queue=UpdateQueueFromRemote(cart_distances, rNodeDistances);
        //Sanity - check that we aren't doing this very many times
        assert(round_counter++ <= PetscTools::GetNumProcs()+2);
    }


    if (mWorkOnEntireMesh==false)
    {
        //Update all processes with the best values from everywhere
        //Take a local copy
        std::vector<double> local_distances=rNodeDistances;
        //Share it back into the vector
        MPI_Allreduce( &local_distances[0], &rNodeDistances[0], mNumNodes, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::UpdateQueueFromRemote(std::vector< c_vector<double, SPACE_DIM> >& rCartDistances,
                                                                     std::vector<double>& rNodeDistances)
{
    if (mWorkOnEntireMesh)
    {
        //This update does nowt.
        return false;
    }
    for (unsigned bcast_process=0; bcast_process<PetscTools::GetNumProcs(); bcast_process++)
    {
        //Process packs cart0/cart1/...euclid/index into a 1-d array
        double* cart_exchange = new double[ SPACE_DIM * mNumHalosPerProcess[bcast_process] ];
        double* dist_exchange = new double[  mNumHalosPerProcess[bcast_process] ];
        unsigned* index_exchange = new unsigned[  mNumHalosPerProcess[bcast_process] ];
        if (PetscTools::GetMyRank() == bcast_process)
        {
            //Broadcaster fills the array
            for (unsigned index=0; index<mHaloNodeIndices.size();index++)
            {
                for (unsigned j=0; j<SPACE_DIM; j++)
                {
                    cart_exchange[ index*(SPACE_DIM) + j] =  rCartDistances[mHaloNodeIndices[index]][j];
                }
                dist_exchange[index] = rNodeDistances[mHaloNodeIndices[index]];
                index_exchange[index] = (double) mHaloNodeIndices[index]; 
            }
        }

        //Broadcast - this is can be done by casting indices to double and packing everything
        //into a single array.  That would be better for latency, but this is probably more readable.
        MPI_Bcast(cart_exchange, (SPACE_DIM) * mNumHalosPerProcess[bcast_process], MPI_DOUBLE, 
                  bcast_process, PETSC_COMM_WORLD);
        MPI_Bcast(dist_exchange, mNumHalosPerProcess[bcast_process], MPI_DOUBLE, 
                  bcast_process, PETSC_COMM_WORLD);
        MPI_Bcast(index_exchange, mNumHalosPerProcess[bcast_process], MPI_UNSIGNED, 
                  bcast_process, PETSC_COMM_WORLD);
        if (PetscTools::GetMyRank() != bcast_process)
        {
            //Receiving process take updates
            for (unsigned index=0; index<mNumHalosPerProcess[bcast_process];index++)
            {
                unsigned global_index=index_exchange[index];
                //Is it a better answer?
                if (dist_exchange[index] < rNodeDistances[global_index]*(1.0-DBL_EPSILON))
                {
                    //Copy across - this may be unnecessary when PushLocal isn't going to push because it's not local.
                    rNodeDistances[global_index] = dist_exchange[index];
                    for (unsigned j=0; j<SPACE_DIM; j++)
                    {
                         rCartDistances[global_index][j] = cart_exchange[ index*(SPACE_DIM) + j];
                    }
                    PushLocal(global_index);
                }
            }
        }
        delete [] cart_exchange;
        delete [] dist_exchange;
        delete [] index_exchange;
    }
    //Is any queue non-empty?
    bool non_empty_queue=PetscTools::ReplicateBool(!mActiveNodeIndexQueue.empty());
    return(non_empty_queue);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::WorkOnLocalQueue(std::vector< c_vector<double, SPACE_DIM> >& rCartDistances,
                                                                     std::vector<double>& rNodeDistances)
{
   while (!mActiveNodeIndexQueue.empty())
    {
        // Get the next index in the queue
        unsigned current_node_index = mActiveNodeIndexQueue.front();
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
                            rCartDistances[neighbour_node_index] = rCartDistances[current_node_index]
                                                                   + (p_current_node->rGetLocation() - p_neighbour_node->rGetLocation());
                            //This will save some sqrts later...
                            rNodeDistances[neighbour_node_index] = norm_2(rCartDistances[neighbour_node_index]);                                      
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
