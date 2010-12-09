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

#include "AbstractCardiacTissue.hpp"

#include "DistributedVector.hpp"
#include "AxisymmetricConductivityTensors.hpp"
#include "OrthotropicConductivityTensors.hpp"
#include "ChastePoint.hpp"
#include "AbstractChasteRegion.hpp"
#include "HeartEventHandler.hpp"
#include "PetscTools.hpp"
#include "PetscVecTools.hpp"

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::AbstractCardiacTissue(
            AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory,
            bool exchangeHalos)
    : mpMesh(pCellFactory->GetMesh()),
      mDoCacheReplication(true),
      mpDistributedVectorFactory(mpMesh->GetDistributedVectorFactory()),
      mMeshUnarchived(false),
      mpConductivityModifier(NULL),
      mExchangeHalos(exchangeHalos)
{
    //This constructor is called from the Initialise() method of the CardiacProblem class
    assert(pCellFactory != NULL);
    assert(pCellFactory->GetMesh() != NULL);

    if (mExchangeHalos)
    {
        mpMesh->CalculateNodeExchange(mNodesToSendPerProcess, mNodesToReceivePerProcess);
    }

    unsigned num_local_nodes = mpDistributedVectorFactory->GetLocalOwnership();
    unsigned ownership_range_low = mpDistributedVectorFactory->GetLow();
    mCellsDistributed.resize(num_local_nodes);

    try
    {
        for (unsigned local_index = 0; local_index < num_local_nodes; local_index++)
        {
            unsigned global_index = ownership_range_low + local_index;
            mCellsDistributed[local_index] = pCellFactory->CreateCardiacCellForNode(global_index);
            mCellsDistributed[local_index]->SetUsedInTissueSimulation();
        }

        pCellFactory->FinaliseCellCreation(&mCellsDistributed,
                                           mpDistributedVectorFactory->GetLow(),
                                           mpDistributedVectorFactory->GetHigh());
    }
    catch (const Exception& e)
    {
        // Errors thrown creating cells will often be process-specific
        PetscTools::ReplicateException(true);

        // Delete cells
        // Should really do this for other processes too, but this is all we need
        // to get memory testing to pass, and leaking when we're about to die isn't
        // that bad!
        for (std::vector<AbstractCardiacCell*>::iterator cell_iterator = mCellsDistributed.begin();
             cell_iterator != mCellsDistributed.end();
             ++cell_iterator)
        {
            delete (*cell_iterator);
        }

        throw e;
    }
    PetscTools::ReplicateException(false);

    // Halo nodes (if required)
    if (mExchangeHalos)
    {
        mpMesh->GetHaloNodeIndices( mHaloNodes );
        unsigned num_halo_nodes = mHaloNodes.size();
        mHaloCellsDistributed.resize( num_halo_nodes );

        try
        {
            for (unsigned local_index = 0; local_index < num_halo_nodes; local_index++)
            {
                unsigned global_index = mHaloNodes[local_index];
                mHaloCellsDistributed[local_index] = pCellFactory->CreateCardiacCellForNode(global_index);
                mHaloCellsDistributed[local_index]->SetUsedInTissueSimulation();
                mHaloGlobalToLocalIndexMap[global_index] = local_index;
            }

            // No need to call FinaliseCellCreation() as halo node cardiac cells will
            // never be stimulated (their values are communicated from the process that
            // owns them.
        }
        catch (const Exception& e)
        {
            // Errors thrown creating cells will often be process-specific
            PetscTools::ReplicateException(true);

            // Delete cells
            // Should really do this for other processes too, but this is all we need
            // to get memory testing to pass, and leaking when we're about to die isn't
            // that bad!
            for (std::vector<AbstractCardiacCell*>::iterator cell_iterator = mHaloCellsDistributed.begin();
                 cell_iterator != mHaloCellsDistributed.end();
                 ++cell_iterator)
            {
                delete (*cell_iterator);
            }

            throw e;
        }
        PetscTools::ReplicateException(false);
    }

    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    mIionicCacheReplicated.Resize( pCellFactory->GetNumberOfCells() );
    mIntracellularStimulusCacheReplicated.Resize( pCellFactory->GetNumberOfCells() );
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);

    if(HeartConfig::Instance()->IsMeshProvided() && HeartConfig::Instance()->GetLoadMesh())
    {
        mFibreFilePathNoExtension = "./" + HeartConfig::Instance()->GetMeshName();
    }
    else
    {
        // As of 10671 fibre orientation can only be defined when loading a mesh from disc.
        mFibreFilePathNoExtension = "";
    }
    CreateIntracellularConductivityTensor();
}

// Constructor used for archiving
template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::AbstractCardiacTissue(std::vector<AbstractCardiacCell*> & rCellsDistributed,
                                                                    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    : mpMesh(pMesh),
      mCellsDistributed(rCellsDistributed),
      mDoCacheReplication(true),
      mpDistributedVectorFactory(mpMesh->GetDistributedVectorFactory()),
      mMeshUnarchived(true)
{
    mIionicCacheReplicated.Resize(mpDistributedVectorFactory->GetProblemSize());
    mIntracellularStimulusCacheReplicated.Resize(mpDistributedVectorFactory->GetProblemSize());

    mFibreFilePathNoExtension = ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename();
    CreateIntracellularConductivityTensor();
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::~AbstractCardiacTissue()
{
    // Delete cells
    for (std::vector<AbstractCardiacCell*>::iterator cell_iterator = mCellsDistributed.begin();
         cell_iterator != mCellsDistributed.end();
         ++cell_iterator)
    {
        delete (*cell_iterator);
    }

    // Delete cells for halo nodes
    for (std::vector<AbstractCardiacCell*>::iterator cell_iterator = mHaloCellsDistributed.begin();
         cell_iterator != mHaloCellsDistributed.end();
         ++cell_iterator)
    {
        delete (*cell_iterator);
    }

    delete mpIntracellularConductivityTensors;

    // If the mesh was unarchived we need to free it explicitly.
    if (mMeshUnarchived)
    {
        delete mpMesh;
    }
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::MergeCells(const std::vector<AbstractCardiacCell*>& rOtherCells)
{
    assert(rOtherCells.size() == mCellsDistributed.size());
    for (unsigned i=0; i<rOtherCells.size(); i++)
    {
        if (rOtherCells[i] != NULL)
        {
            assert(mCellsDistributed[i] == NULL);
            mCellsDistributed[i] = rOtherCells[i];
        }
    }
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::CreateIntracellularConductivityTensor()
{
    HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
    mpConfig = HeartConfig::Instance();

    if (mpConfig->IsMeshProvided() && mpConfig->GetLoadMesh())
    {
        assert(mFibreFilePathNoExtension != "");

        switch (mpConfig->GetConductivityMedia())
        {
            case cp::media_type::Orthotropic:
            {
                mpIntracellularConductivityTensors = new OrthotropicConductivityTensors<ELEMENT_DIM,SPACE_DIM>;
                FileFinder ortho_file(mFibreFilePathNoExtension + ".ortho", RelativeTo::AbsoluteOrCwd);
                assert(ortho_file.Exists());
                mpIntracellularConductivityTensors->SetFibreOrientationFile(ortho_file);
                break;
            }

            case cp::media_type::Axisymmetric:
            {
                mpIntracellularConductivityTensors = new AxisymmetricConductivityTensors<ELEMENT_DIM,SPACE_DIM>;
                FileFinder axi_file(mFibreFilePathNoExtension + ".axi", RelativeTo::AbsoluteOrCwd);
                assert(axi_file.Exists());
                mpIntracellularConductivityTensors->SetFibreOrientationFile(axi_file);
                break;
            }

            case cp::media_type::NoFibreOrientation:
                /// \todo #1316 Create a class defining constant tensors to be used when no fibre orientation is provided.
                mpIntracellularConductivityTensors = new OrthotropicConductivityTensors<ELEMENT_DIM,SPACE_DIM>;
                break;

            default :
                NEVER_REACHED;
        }
    }
    else // Slab defined in config file or SetMesh() called; no fibre orientation assumed
    {
        /// \todo #1316 Create a class defining constant tensors to be used when no fibre orientation is provided.
        mpIntracellularConductivityTensors = new OrthotropicConductivityTensors<ELEMENT_DIM,SPACE_DIM>;
    }

    c_vector<double, SPACE_DIM> intra_conductivities;
    mpConfig->GetIntracellularConductivities(intra_conductivities);

    // this definition must be here (and not inside the if statement) because SetNonConstantConductivities() will keep
    // a pointer to it and we don't want it to go out of scope before Init() is called
    unsigned num_local_elements = mpMesh->GetNumLocalElements();
    std::vector<c_vector<double, SPACE_DIM> > hetero_intra_conductivities;

    if (mpConfig->GetConductivityHeterogeneitiesProvided())
    {
        try
        {
            assert(hetero_intra_conductivities.size()==0);
            hetero_intra_conductivities.resize(num_local_elements, intra_conductivities);
        }
        catch(std::bad_alloc &r_bad_alloc)
        {
#define COVERAGE_IGNORE
            std::cout << "Failed to allocate std::vector of size " << num_local_elements << std::endl;
            PetscTools::ReplicateException(true);
            throw r_bad_alloc;
#undef COVERAGE_IGNORE
        }
        PetscTools::ReplicateException(false);

        std::vector<AbstractChasteRegion<SPACE_DIM>* > conductivities_heterogeneity_areas;
        std::vector< c_vector<double,3> > intra_h_conductivities;
        std::vector< c_vector<double,3> > extra_h_conductivities;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                intra_h_conductivities,
                                                                extra_h_conductivities);

        unsigned local_element_index = 0;

        for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator it = mpMesh->GetElementIteratorBegin();
             it != mpMesh->GetElementIteratorEnd();
             ++it)
        {
//            unsigned element_index = it->GetIndex();
            // if element centroid is contained in the region
            ChastePoint<SPACE_DIM> element_centroid(it->CalculateCentroid());
            for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas.size(); region_index++)
            {
                if ( conductivities_heterogeneity_areas[region_index]->DoesContain(element_centroid) )
                {
                    hetero_intra_conductivities[local_element_index] = intra_h_conductivities[region_index];
                }
            }
            local_element_index++;
        }

        // freeing memory allocated by HeartConfig::Instance()->GetConductivityHeterogeneities
        for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas.size(); region_index++)
        {
            delete conductivities_heterogeneity_areas[region_index];
        }

        mpIntracellularConductivityTensors->SetNonConstantConductivities(&hetero_intra_conductivities);
    }
    else
    {
        mpIntracellularConductivityTensors->SetConstantConductivities(intra_conductivities);
    }

    mpIntracellularConductivityTensors->Init(this->mpMesh);
    HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);
}


template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::SetCacheReplication(bool doCacheReplication)
{
    mDoCacheReplication = doCacheReplication;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
bool AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::GetDoCacheReplication()
{
    return mDoCacheReplication;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
const c_matrix<double, SPACE_DIM, SPACE_DIM>& AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::rGetIntracellularConductivityTensor(unsigned elementIndex)
{
    assert( mpIntracellularConductivityTensors);
    if(mpConductivityModifier==NULL)
    {
        return (*mpIntracellularConductivityTensors)[elementIndex];
    }
    else
    {
        return mpConductivityModifier->rGetModifiedConductivityTensor(elementIndex, (*mpIntracellularConductivityTensors)[elementIndex]);
    }
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacCell* AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::GetCardiacCell( unsigned globalIndex )
{
    assert(mpDistributedVectorFactory->GetLow() <= globalIndex &&
           globalIndex < mpDistributedVectorFactory->GetHigh());
    return mCellsDistributed[globalIndex - mpDistributedVectorFactory->GetLow()];
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacCell* AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::GetCardiacCellOrHaloCell( unsigned globalIndex )
{
    std::map<unsigned, unsigned>::const_iterator node_position;
    //First search the halo
    if ((node_position=mHaloGlobalToLocalIndexMap.find(globalIndex)) != mHaloGlobalToLocalIndexMap.end())
    {
        return mHaloCellsDistributed[node_position->second];
    }
    //Then search the owned node
    if ( mpDistributedVectorFactory->GetLow() <= globalIndex &&
         globalIndex < mpDistributedVectorFactory->GetHigh()    )
    {
        //Found an owned node
        return this->mCellsDistributed[node_position->second];
    }
    //Not here
    std::stringstream message;
    message << "Requested node/halo " << globalIndex << " does not belong to processor " << PetscTools::GetMyRank();
    EXCEPTION(message.str().c_str());
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::SolveCellSystems(Vec existingSolution, double time, double nextTime, bool updateVoltage)
{
    HeartEventHandler::BeginEvent(HeartEventHandler::SOLVE_ODES);

    DistributedVector dist_solution = mpDistributedVectorFactory->CreateDistributedVector(existingSolution);
    DistributedVector::Stripe voltage(dist_solution, 0);

    try
    {
        for (DistributedVector::Iterator index = dist_solution.Begin();
             index != dist_solution.End();
             ++index)
        {
            // overwrite the voltage with the input value
            mCellsDistributed[index.Local]->SetVoltage( voltage[index] );

            if(!updateVoltage)
            {
                // solve
                // Note: Voltage is not be updated. The voltage is updated in the PDE solve.
                mCellsDistributed[index.Local]->ComputeExceptVoltage(time, nextTime);
            }
            else
            {
                // solve, including updating the voltage (for the operator-splitting implementation of the monodomain solver)
                mCellsDistributed[index.Local]->SolveAndUpdateState(time, nextTime);
                voltage[index] = mCellsDistributed[index.Local]->GetVoltage();
            }

            // update the Iionic and stimulus caches
            UpdateCaches(index.Global, index.Local, nextTime);
        }

        if(updateVoltage)
        {
        	dist_solution.Restore();
        }
    }
    catch (Exception &e)
    {
        PetscTools::ReplicateException(true);
        throw e;
    }

    PetscTools::ReplicateException(false);
    HeartEventHandler::EndEvent(HeartEventHandler::SOLVE_ODES);

    // Communicate new state variable values to halo nodes
    if (mExchangeHalos)
    {
        for ( unsigned rank_offset = 1; rank_offset < PetscTools::GetNumProcs(); rank_offset++ )
        {
            unsigned send_to      = (PetscTools::GetMyRank() + rank_offset) % (PetscTools::GetNumProcs());
            unsigned receive_from = (PetscTools::GetMyRank() + PetscTools::GetNumProcs()- rank_offset ) % (PetscTools::GetNumProcs());

            unsigned number_of_cells_to_send    = mNodesToSendPerProcess[send_to].size();
            unsigned number_of_cells_to_receive = mNodesToReceivePerProcess[receive_from].size();

            unsigned number_of_state_variables = mCellsDistributed[0]->GetNumberOfStateVariables();
            unsigned send_size = number_of_cells_to_send*number_of_state_variables;
            unsigned receive_size = number_of_cells_to_receive*number_of_state_variables;

            // Pack
            if ( number_of_cells_to_send > 0 )
            {
                double send_data[send_size];
//                send_data.resize(send_size);

                for ( unsigned cell = 0; cell < number_of_cells_to_send; cell++ )
                {
                    std::vector<double>& cell_data = mCellsDistributed[mHaloGlobalToLocalIndexMap[mNodesToSendPerProcess[send_to][cell]]]->rGetStateVariables();
                    for (unsigned state_variable = 0; state_variable < number_of_state_variables; state_variable++)
                    {
                        send_data[ cell*number_of_state_variables + state_variable ] = cell_data[state_variable];
                    }
                }

                // Send
                MPI_Send( &send_data[0],
                          send_size,
                          MPI_DOUBLE,
                          send_to,
                          0,
                          PETSC_COMM_WORLD );
            }

            if ( number_of_cells_to_receive > 0 )
            {
                // Receive
                double receive_data[receive_size];
                MPI_Status status;

                MPI_Recv( receive_data,
                          receive_size,
                          MPI_DOUBLE,
                          receive_from,
                          0,
                          PETSC_COMM_WORLD,
                          &status );


                // Unpack
                for ( unsigned cell = 0; cell < number_of_cells_to_receive; cell++ )
                {
                    std::vector<double> cell_data;
                    cell_data.resize(number_of_state_variables);
                    for (unsigned state_variable = 0; state_variable < number_of_state_variables; state_variable++)
                    {
                        cell_data[state_variable] = receive_data[ cell*number_of_state_variables + state_variable ];
                    }
                    mHaloCellsDistributed[mHaloGlobalToLocalIndexMap[mNodesToReceivePerProcess[receive_from][cell]]]->SetStateVariables(cell_data);
                }
            }
        }
    }

    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    if ( mDoCacheReplication )
    {
        ReplicateCaches();
    }
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
ReplicatableVector& AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::rGetIionicCacheReplicated()
{
    return mIionicCacheReplicated;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
ReplicatableVector& AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::rGetIntracellularStimulusCacheReplicated()
{
    return mIntracellularStimulusCacheReplicated;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::UpdateCaches(unsigned globalIndex, unsigned localIndex, double nextTime)
{
    mIionicCacheReplicated[globalIndex] = mCellsDistributed[localIndex]->GetIIonic();
    mIntracellularStimulusCacheReplicated[globalIndex] = mCellsDistributed[localIndex]->GetIntracellularStimulus(nextTime);
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::ReplicateCaches()
{
    mIionicCacheReplicated.Replicate(mpDistributedVectorFactory->GetLow(), mpDistributedVectorFactory->GetHigh());
    mIntracellularStimulusCacheReplicated.Replicate(mpDistributedVectorFactory->GetLow(), mpDistributedVectorFactory->GetHigh());
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
const std::vector<AbstractCardiacCell*>& AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::rGetCellsDistributed() const
{
    return mCellsDistributed;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
const AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::pGetMesh() const
{
    return mpMesh;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::SetConductivityModifier(AbstractConductivityModifier<ELEMENT_DIM,SPACE_DIM>* pModifier)
{
    assert(pModifier!=NULL);
    assert(mpConductivityModifier==NULL); // shouldn't be called twice for example, or with two different modifiers (remove this assert
                                          // if for whatever reason want to be able to overwrite modifiers
    mpConductivityModifier = pModifier;
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractCardiacTissue<1,1>;
template class AbstractCardiacTissue<1,2>;
template class AbstractCardiacTissue<1,3>;
template class AbstractCardiacTissue<2,2>;
template class AbstractCardiacTissue<3,3>;
