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

#include "AbstractCardiacPde.hpp"

#include "DistributedVector.hpp"
#include "AxisymmetricConductivityTensors.hpp"
#include "OrthotropicConductivityTensors.hpp"
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"
#include "HeartEventHandler.hpp"
#include "PetscTools.hpp"
#include "FakeBathCell.hpp"


template <unsigned SPACE_DIM>
AbstractCardiacPde<SPACE_DIM>::AbstractCardiacPde(
            AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory,
            const unsigned stride)
    : mStride(stride),
      mDoCacheReplication(true),
      mDoOneCacheReplication(true)
{
    //This constructor is called from the Initialise() method of the CardiacProblem class
    assert(pCellFactory!=NULL);
    assert(pCellFactory->GetMesh()!=NULL);

    std::vector<unsigned>& r_nodes_per_processor = pCellFactory->GetMesh()->rGetNodesPerProcessor();

    // check number of processor agrees with definition in mesh
    if((r_nodes_per_processor.size() != 0) && (r_nodes_per_processor.size() != PetscTools::NumProcs()) )
    {
        EXCEPTION("Number of processors defined in mesh class not equal to number of processors used");
    }

    if(r_nodes_per_processor.size() != 0)
    {
        unsigned num_local_nodes = r_nodes_per_processor[ PetscTools::GetMyRank() ];
        DistributedVector::SetProblemSizePerProcessor(pCellFactory->GetMesh()->GetNumNodes(), num_local_nodes);
    }
    else
    {
        DistributedVector::SetProblemSize(pCellFactory->GetMesh()->GetNumNodes());
    }

    mCellsDistributed.resize(DistributedVector::End().Global-DistributedVector::Begin().Global);

    for (DistributedVector::Iterator index = DistributedVector::Begin();
         index != DistributedVector::End();
         ++index)
    {
        mCellsDistributed[index.Local] = pCellFactory->CreateCardiacCellForNode(index.Global);
    }
    pCellFactory->FinaliseCellCreation(&mCellsDistributed,
                                       DistributedVector::Begin().Global,
                                       DistributedVector::End().Global);


    mIionicCacheReplicated.resize( pCellFactory->GetNumberOfCells() );
    mIntracellularStimulusCacheReplicated.resize( pCellFactory->GetNumberOfCells() );

    mpConfig = HeartConfig::Instance();

    if (mpConfig->GetIsMeshProvided() && mpConfig->GetLoadMesh())
    {
        switch (mpConfig->GetConductivityMedia())
        {
            case media_type::Orthotropic:
                mpIntracellularConductivityTensors =  new OrthotropicConductivityTensors<SPACE_DIM>;
                mpIntracellularConductivityTensors->SetFibreOrientationFile(mpConfig->GetMeshName() + ".ortho");
                break;

            case media_type::Axisymmetric:
                mpIntracellularConductivityTensors =  new AxisymmetricConductivityTensors<SPACE_DIM>;
                mpIntracellularConductivityTensors->SetFibreOrientationFile(mpConfig->GetMeshName() + ".axi");
                break;

            case media_type::NoFibreOrientation:
                /// \todo: Create a class defining constant tensors to be used when no fibre orientation is provided
                mpIntracellularConductivityTensors =  new OrthotropicConductivityTensors<SPACE_DIM>;
                break;

            default :
                NEVER_REACHED;
        }
    }
    else // Slab defined in config file or SetMesh() called; no fibre orientation assumed
    {
        // See previous todo.
    mpIntracellularConductivityTensors =  new OrthotropicConductivityTensors<SPACE_DIM>;
    }

    c_vector<double, SPACE_DIM> intra_conductivities;
    mpConfig->GetIntracellularConductivities(intra_conductivities);

    // this definition must be here (and not inside the if statement) because SetNonConstantConductivities() will keep
    // a pointer to it and we don't want it to go out of scope before Init() is called
    unsigned num_elements = pCellFactory->GetMesh()->GetNumElements();
    std::vector<c_vector<double, SPACE_DIM> > hetero_intra_conductivities(num_elements);

    if (mpConfig->GetConductivityHeterogeneitiesProvided())
    {
        std::vector<ChasteCuboid> conductivities_heterogeneity_areas;
        std::vector< c_vector<double,3> > intra_h_conductivities;
        std::vector< c_vector<double,3> > extra_h_conductivities;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                intra_h_conductivities,
                                                                extra_h_conductivities);

        for (unsigned element_index=0; element_index<num_elements; element_index++)
        {
            for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas.size(); region_index++)
            {
                // if element centroid is contained in the region
                ChastePoint<SPACE_DIM> element_centroid(pCellFactory->GetMesh()->GetElement(element_index)->CalculateCentroid());
                if ( conductivities_heterogeneity_areas[region_index].DoesContain(element_centroid) )
                {
                    hetero_intra_conductivities[element_index] = intra_h_conductivities[region_index];
                }
                else
                {
                    hetero_intra_conductivities[element_index] = intra_conductivities;
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
}

template <unsigned SPACE_DIM>
AbstractCardiacPde<SPACE_DIM>::~AbstractCardiacPde()
{
    for (DistributedVector::Iterator index = DistributedVector::Begin();
         index != DistributedVector::End();
         ++index)
    {
        // Only delete real cells
        FakeBathCell* p_fake = dynamic_cast<FakeBathCell*>(mCellsDistributed[index.Local]);
        if (p_fake == NULL)
        {
            delete mCellsDistributed[index.Local];
        }
    }

    delete mpIntracellularConductivityTensors;
}

template <unsigned SPACE_DIM>
void AbstractCardiacPde<SPACE_DIM>::SetCacheReplication(bool doCacheReplication)
{
    mDoCacheReplication = doCacheReplication;
}

template <unsigned SPACE_DIM>
const c_matrix<double, SPACE_DIM, SPACE_DIM>& AbstractCardiacPde<SPACE_DIM>::rGetIntracellularConductivityTensor(unsigned elementIndex)
{
    assert( mpIntracellularConductivityTensors);
    return (*mpIntracellularConductivityTensors)[elementIndex];
}

template <unsigned SPACE_DIM>
AbstractCardiacCell* AbstractCardiacPde<SPACE_DIM>::GetCardiacCell( unsigned globalIndex )
{
    assert(DistributedVector::Begin().Global <= globalIndex &&
           globalIndex < DistributedVector::End().Global);
    return mCellsDistributed[globalIndex - DistributedVector::Begin().Global];
}


template <unsigned SPACE_DIM>
void AbstractCardiacPde<SPACE_DIM>::SolveCellSystems(Vec currentSolution, double currentTime, double nextTime)
{
    HeartEventHandler::BeginEvent(HeartEventHandler::SOLVE_ODES);

    DistributedVector dist_solution(currentSolution);
    DistributedVector::Stripe voltage(dist_solution, 0);
    for (DistributedVector::Iterator index = DistributedVector::Begin();
         index != DistributedVector::End();
         ++index)
    {
        // overwrite the voltage with the input value
        mCellsDistributed[index.Local]->SetVoltage( voltage[index] );
        try
        {
            // solve
            // Note: Voltage should not be updated. GetIIonic will be called later
            // and needs the old voltage. The voltage will be updated from the pde.
            mCellsDistributed[index.Local]->ComputeExceptVoltage(currentTime, nextTime);
        }
        catch (Exception &e)
        {
            PetscTools::ReplicateException(true);
            throw e;
        }

        // update the Iionic and stimulus caches
        UpdateCaches(index.Global, index.Local, nextTime);
    }
    HeartEventHandler::EndEvent(HeartEventHandler::SOLVE_ODES);

    PetscTools::ReplicateException(false);

    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    if ((mDoCacheReplication)||mDoOneCacheReplication)
    {
        ReplicateCaches();
        mDoOneCacheReplication=false;
    }
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
}

template <unsigned SPACE_DIM>
ReplicatableVector& AbstractCardiacPde<SPACE_DIM>::rGetIionicCacheReplicated()
{
    return mIionicCacheReplicated;
}

template <unsigned SPACE_DIM>
ReplicatableVector& AbstractCardiacPde<SPACE_DIM>::rGetIntracellularStimulusCacheReplicated()
{
    return mIntracellularStimulusCacheReplicated;
}

template <unsigned SPACE_DIM>
void AbstractCardiacPde<SPACE_DIM>::UpdateCaches(unsigned globalIndex, unsigned localIndex, double nextTime)
{
    mIionicCacheReplicated[globalIndex] = mCellsDistributed[localIndex]->GetIIonic();
    mIntracellularStimulusCacheReplicated[globalIndex] = mCellsDistributed[localIndex]->GetIntracellularStimulus(nextTime);
}

template <unsigned SPACE_DIM>
void AbstractCardiacPde<SPACE_DIM>::ReplicateCaches()
{
    mIionicCacheReplicated.Replicate(DistributedVector::Begin().Global, DistributedVector::End().Global);
    mIntracellularStimulusCacheReplicated.Replicate(DistributedVector::Begin().Global, DistributedVector::End().Global);
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractCardiacPde<1>;
template class AbstractCardiacPde<2>;
template class AbstractCardiacPde<3>;
