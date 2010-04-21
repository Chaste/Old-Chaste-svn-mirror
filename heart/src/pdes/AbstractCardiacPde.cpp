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

#include "AbstractCardiacPde.hpp"

#include "DistributedVector.hpp"
#include "AxisymmetricConductivityTensors.hpp"
#include "OrthotropicConductivityTensors.hpp"
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"
#include "HeartEventHandler.hpp"
#include "PetscTools.hpp"
#include "FakeBathCell.hpp"


template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::AbstractCardiacPde(
            AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory,
            const unsigned stride)
    : mpMesh(pCellFactory->GetMesh()),
      mStride(stride),
      mDoCacheReplication(true),
      mDoOneCacheReplication(true),
      mpDistributedVectorFactory(mpMesh->GetDistributedVectorFactory()),
      mMeshUnarchived(false)
{
    //This constructor is called from the Initialise() method of the CardiacProblem class
    assert(pCellFactory != NULL);
    assert(pCellFactory->GetMesh() != NULL);

    unsigned num_local_nodes = mpDistributedVectorFactory->GetLocalOwnership();
    unsigned ownership_range_low = mpDistributedVectorFactory->GetLow();
    mCellsDistributed.resize(num_local_nodes);

    for (unsigned local_index = 0; local_index < num_local_nodes; local_index++)
    {
        unsigned global_index = ownership_range_low + local_index;
        mCellsDistributed[local_index] = pCellFactory->CreateCardiacCellForNode(global_index);
    }

    pCellFactory->FinaliseCellCreation(&mCellsDistributed,
                                       mpDistributedVectorFactory->GetLow(),
                                       mpDistributedVectorFactory->GetHigh());

    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    mIionicCacheReplicated.Resize( pCellFactory->GetNumberOfCells() );
    mIntracellularStimulusCacheReplicated.Resize( pCellFactory->GetNumberOfCells() );
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);

    CreateIntracellularConductivityTensor();
}

// Constructor used for archiving
template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::AbstractCardiacPde(std::vector<AbstractCardiacCell*> & rCellsDistributed,
                                                              AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                                              const unsigned stride)
    : mpMesh(pMesh),
      mCellsDistributed(rCellsDistributed),
      mStride(stride),
      mDoCacheReplication(true),
      mDoOneCacheReplication(true),
      mpDistributedVectorFactory(mpMesh->GetDistributedVectorFactory()),
      mMeshUnarchived(true)
{
    mIionicCacheReplicated.Resize(mpDistributedVectorFactory->GetProblemSize());
    mIntracellularStimulusCacheReplicated.Resize(mpDistributedVectorFactory->GetProblemSize());

    CreateIntracellularConductivityTensor();
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::DeleteCells(bool deleteFakeCells)
{
    std::set<FakeBathCell*> fake_cells;
    for (std::vector<AbstractCardiacCell*>::iterator cell_iterator = mCellsDistributed.begin();
         cell_iterator != mCellsDistributed.end();
         ++cell_iterator)
    {
        // Only delete real cells, unless we were loaded from an archive
        FakeBathCell* p_fake = dynamic_cast<FakeBathCell*>(*cell_iterator);
        if (p_fake == NULL)
        {
            delete (*cell_iterator);
        }
        else
        {
            fake_cells.insert(p_fake);
        }
    }

    if (deleteFakeCells)
    {
        // Likewise for fake cells
        for (std::set<FakeBathCell*>::iterator it = fake_cells.begin();
             it != fake_cells.end();
             ++it)
        {
            delete (*it);
        }
    }
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::~AbstractCardiacPde()
{
    DeleteCells(mMeshUnarchived);

    delete mpIntracellularConductivityTensors;

    // If the mesh was unarchived we need to free it explicitly.
    if (mMeshUnarchived)
    {
        delete mpMesh;
    }
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::MergeCells(const std::vector<AbstractCardiacCell*>& rOtherCells)
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
void AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::CreateIntracellularConductivityTensor()
{
    HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
    mpConfig = HeartConfig::Instance();

    if (mpConfig->IsMeshProvided() && mpConfig->GetLoadMesh())
    {
        switch (mpConfig->GetConductivityMedia())
        {
            case cp::media_type::Orthotropic:
                mpIntracellularConductivityTensors = new OrthotropicConductivityTensors<SPACE_DIM>;
                mpIntracellularConductivityTensors->SetFibreOrientationFile(mpConfig->GetMeshName() + ".ortho");
                break;

            case cp::media_type::Axisymmetric:
                mpIntracellularConductivityTensors = new AxisymmetricConductivityTensors<SPACE_DIM>;
                mpIntracellularConductivityTensors->SetFibreOrientationFile(mpConfig->GetMeshName() + ".axi");
                break;

            case cp::media_type::NoFibreOrientation:
                /// \todo #1316 Create a class defining constant tensors to be used when no fibre orientation is provided.
                mpIntracellularConductivityTensors = new OrthotropicConductivityTensors<SPACE_DIM>;
                break;

            default :
                NEVER_REACHED;
        }
    }
    else // Slab defined in config file or SetMesh() called; no fibre orientation assumed
    {
        /// \todo #1316 Create a class defining constant tensors to be used when no fibre orientation is provided.
        mpIntracellularConductivityTensors = new OrthotropicConductivityTensors<SPACE_DIM>;
    }

    c_vector<double, SPACE_DIM> intra_conductivities;
    mpConfig->GetIntracellularConductivities(intra_conductivities);

    // this definition must be here (and not inside the if statement) because SetNonConstantConductivities() will keep
    // a pointer to it and we don't want it to go out of scope before Init() is called
    unsigned num_elements = mpMesh->GetNumElements();
    std::vector<c_vector<double, SPACE_DIM> > hetero_intra_conductivities;

    if (mpConfig->GetConductivityHeterogeneitiesProvided())
    {
        try
        {
            assert(hetero_intra_conductivities.size()==0);
            hetero_intra_conductivities.resize(num_elements, intra_conductivities);
        }
        catch(std::bad_alloc &badAlloc)
        {
#define COVERAGE_IGNORE
            std::cout << "Failed to allocate std::vector of size " << num_elements << std::endl;
            PetscTools::ReplicateException(true);
            throw badAlloc;
#undef COVERAGE_IGNORE
        }
        PetscTools::ReplicateException(false);

        std::vector<ChasteCuboid<SPACE_DIM> > conductivities_heterogeneity_areas;
        std::vector< c_vector<double,3> > intra_h_conductivities;
        std::vector< c_vector<double,3> > extra_h_conductivities;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                intra_h_conductivities,
                                                                extra_h_conductivities);

        for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator it = mpMesh->GetElementIteratorBegin();
             it != mpMesh->GetElementIteratorEnd();
             ++it)
        {
            unsigned element_index = it->GetIndex();            
            // if element centroid is contained in the region
            ChastePoint<SPACE_DIM> element_centroid(it->CalculateCentroid());
            for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas.size(); region_index++)
            {
                if ( conductivities_heterogeneity_areas[region_index].DoesContain(element_centroid) )
                {
                    hetero_intra_conductivities[element_index] = intra_h_conductivities[region_index];
                }
            }
        }

        mpIntracellularConductivityTensors->SetNonConstantConductivities(&hetero_intra_conductivities);
    }
    else
    {
        mpIntracellularConductivityTensors->SetConstantConductivities(intra_conductivities);
    }

    mpIntracellularConductivityTensors->Init();
    HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);
}


template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::SetCacheReplication(bool doCacheReplication)
{
    mDoCacheReplication = doCacheReplication;
    mDoOneCacheReplication = true; //If we have reset to false (even from false) then we may have created a new matrix-based assembler and need to run it once
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
bool AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::GetDoCacheReplication()
{
    return mDoCacheReplication;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
const c_matrix<double, SPACE_DIM, SPACE_DIM>& AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::rGetIntracellularConductivityTensor(unsigned elementIndex)
{
    assert( mpIntracellularConductivityTensors);
    return (*mpIntracellularConductivityTensors)[elementIndex];
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacCell* AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::GetCardiacCell( unsigned globalIndex )
{
    assert(mpDistributedVectorFactory->GetLow() <= globalIndex &&
           globalIndex < mpDistributedVectorFactory->GetHigh());
    return mCellsDistributed[globalIndex - mpDistributedVectorFactory->GetLow()];
}


template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::SolveCellSystems(Vec existingSolution, double time, double nextTime)
{
    HeartEventHandler::BeginEvent(HeartEventHandler::SOLVE_ODES);

    DistributedVector dist_solution = mpDistributedVectorFactory->CreateDistributedVector(existingSolution);
    DistributedVector::Stripe voltage(dist_solution, 0);
    for (DistributedVector::Iterator index = dist_solution.Begin();
         index != dist_solution.End();
         ++index)
    {
        // overwrite the voltage with the input value
        mCellsDistributed[index.Local]->SetVoltage( voltage[index] );
        try
        {
            // solve
            // Note: Voltage should not be updated. GetIIonic will be called later
            // and needs the old voltage. The voltage will be updated from the pde.
            mCellsDistributed[index.Local]->ComputeExceptVoltage(time, nextTime);
        }
        catch (Exception &e)
        {
            PetscTools::ReplicateException(true);
            throw e;
        }

        // update the Iionic and stimulus caches
        UpdateCaches(index.Global, index.Local, nextTime);
    }
    PetscTools::ReplicateException(false);
    HeartEventHandler::EndEvent(HeartEventHandler::SOLVE_ODES);

    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    if ( mDoCacheReplication || mDoOneCacheReplication )
    {
        ReplicateCaches();
        mDoOneCacheReplication=false;
    }
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
ReplicatableVector& AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::rGetIionicCacheReplicated()
{
    return mIionicCacheReplicated;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
ReplicatableVector& AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::rGetIntracellularStimulusCacheReplicated()
{
    return mIntracellularStimulusCacheReplicated;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::UpdateCaches(unsigned globalIndex, unsigned localIndex, double nextTime)
{
    mIionicCacheReplicated[globalIndex] = mCellsDistributed[localIndex]->GetIIonic();
    mIntracellularStimulusCacheReplicated[globalIndex] = mCellsDistributed[localIndex]->GetIntracellularStimulus(nextTime);
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::ReplicateCaches()
{
    mIionicCacheReplicated.Replicate(mpDistributedVectorFactory->GetLow(), mpDistributedVectorFactory->GetHigh());
    mIntracellularStimulusCacheReplicated.Replicate(mpDistributedVectorFactory->GetLow(), mpDistributedVectorFactory->GetHigh());
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
const std::vector<AbstractCardiacCell*>& AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::rGetCellsDistributed() const
{
    return mCellsDistributed;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
const AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>::pGetMesh() const
{
    return mpMesh;
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractCardiacPde<1,1>;
template class AbstractCardiacPde<1,2>;
template class AbstractCardiacPde<1,3>;
template class AbstractCardiacPde<2,2>;
template class AbstractCardiacPde<3,3>;
