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
#include "DistributedTetrahedralMesh.hpp" //For dynamic cast

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::DistanceMapCalculator(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
    : mrMesh(rMesh),
      mWorkOnEntireMesh(true),
      mNumHalosPerProcess(NULL),
      mRoundCounter(0u)
{
    mNumNodes = mrMesh.GetNumNodes();
    DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* p_distributed_mesh = dynamic_cast<DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>*>(&mrMesh);
    if ( PetscTools::IsSequential() || p_distributed_mesh == NULL)
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
        p_distributed_mesh->GetHaloNodeIndices(mHaloNodeIndices);
        //Share information on the number of halo nodes
        unsigned my_size=mHaloNodeIndices.size();
        mNumHalosPerProcess=new unsigned[PetscTools::GetNumProcs()];
        MPI_Allgather(&my_size, 1, MPI_UNSIGNED, 
                     mNumHalosPerProcess, 1, MPI_UNSIGNED, PETSC_COMM_WORLD);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::ComputeDistanceMap(
        const std::vector<unsigned>& rSourceNodeIndices,
        std::vector<double>& rNodeDistances)
{
    rNodeDistances.resize(mNumNodes);
    for (unsigned index=0; index<mNumNodes; index++)
    {
        rNodeDistances[index] = DBL_MAX;
    }
    /*
     * Vector of witness points (a closest point in original source set of points)
     */
    std::vector< c_vector<double, SPACE_DIM> >  witness_points(mNumNodes);

    for (unsigned source_index=0; source_index<rSourceNodeIndices.size(); source_index++)
    {
        unsigned node_index=rSourceNodeIndices[source_index];
        PushLocal(node_index);
        rNodeDistances[node_index] = 0.0;
        //If we have the correct information, then we can set the witness
        if (mLo<=node_index && node_index<mHi)
        {
            witness_points[node_index]=mrMesh.GetNode(node_index)->rGetLocation();
        }
    }

    bool non_empty_queue=true;
    mRoundCounter=0;
    while (non_empty_queue)
    {
        WorkOnLocalQueue(witness_points, rNodeDistances);
        non_empty_queue=UpdateQueueFromRemote(witness_points, rNodeDistances);
        //Sanity - check that we aren't doing this very many times
        if (mRoundCounter++ > 3 * PetscTools::GetNumProcs())
        {
            //This line will be hit if there's a parallel distributed mesh with a really bad partition
            NEVER_REACHED;
        }
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
bool DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::UpdateQueueFromRemote(std::vector< c_vector<double, SPACE_DIM> >& rWitnessPoints,
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
        double* witness_exchange = new double[ SPACE_DIM * mNumHalosPerProcess[bcast_process] ];
        double* dist_exchange = new double[  mNumHalosPerProcess[bcast_process] ];
        unsigned* index_exchange = new unsigned[  mNumHalosPerProcess[bcast_process] ];
        if (PetscTools::GetMyRank() == bcast_process)
        {
            //Broadcaster fills the array
            for (unsigned index=0; index<mHaloNodeIndices.size();index++)
            {
                for (unsigned j=0; j<SPACE_DIM; j++)
                {
                    witness_exchange[ index*(SPACE_DIM) + j] =  rWitnessPoints[mHaloNodeIndices[index]][j];
                }
                dist_exchange[index] = rNodeDistances[mHaloNodeIndices[index]];
                index_exchange[index] = mHaloNodeIndices[index]; 
            }
        }

        //Broadcast - this is can be done by casting indices to double and packing everything
        //into a single array.  That would be better for latency, but this is probably more readable.
        MPI_Bcast(witness_exchange, (SPACE_DIM) * mNumHalosPerProcess[bcast_process], MPI_DOUBLE, 
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
                if (dist_exchange[index] < rNodeDistances[global_index]*(1.0-DBL_EPSILON) )
                {
                    //Copy across - this may be unnecessary when PushLocal isn't going to push because it's not local.
                    rNodeDistances[global_index] = dist_exchange[index];
                    for (unsigned j=0; j<SPACE_DIM; j++)
                    {
                         rWitnessPoints[global_index][j] = witness_exchange[ index*(SPACE_DIM) + j];
                    }
                    PushLocal(global_index);
                }
            }
        }
        delete [] witness_exchange;
        delete [] dist_exchange;
        delete [] index_exchange;
    }
    //Is any queue non-empty?
    bool non_empty_queue=PetscTools::ReplicateBool(!mActiveNodeIndexQueue.empty());
    return(non_empty_queue);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::WorkOnLocalQueue(std::vector< c_vector<double, SPACE_DIM> >& rWitnessPoints,
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
                        
                        // Test if we have found a shorter path from the witness in the source to the current neighbour through current node  
                        //This will save some sqrts later...
                        double updated_distance = norm_2(p_neighbour_node->rGetLocation() - rWitnessPoints[current_node_index]);
                        if ( updated_distance < rNodeDistances[neighbour_node_index] * (1.0-DBL_EPSILON) )
                        {
                            rWitnessPoints[neighbour_node_index] = rWitnessPoints[current_node_index];
                            rNodeDistances[neighbour_node_index] = updated_distance;                                      
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
