/*

Copyright (C) University of Oxford, 2005-2011

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

#include "ExtendedBidomainTissue.hpp"

#include "DistributedVector.hpp"
#include "OrthotropicConductivityTensors.hpp"
#include "AbstractStimulusFunction.hpp"
#include "ChastePoint.hpp"
#include "AbstractChasteRegion.hpp"
#include "HeartEventHandler.hpp"

template <unsigned SPACE_DIM>
ExtendedBidomainTissue<SPACE_DIM>::ExtendedBidomainTissue(
            AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, AbstractCardiacCellFactory<SPACE_DIM>* pCellFactorySecondCell, AbstractStimulusFactory<SPACE_DIM>* pExtracellularStimulusFactory)
    : AbstractCardiacTissue<SPACE_DIM>(pCellFactory),
      mpIntracellularConductivityTensorsSecondCell(NULL),
      mUserSuppliedExtracellularStimulus(false)
{
	//First, do the same that the abstract constructor does, but applied to the second cell

    assert(pCellFactorySecondCell != NULL);
    assert(pCellFactorySecondCell->GetMesh() != NULL);
    assert(pCellFactorySecondCell->GetNumberOfCells() == pCellFactory->GetNumberOfCells() );
    assert(pExtracellularStimulusFactory != NULL);
    assert(pExtracellularStimulusFactory->GetMesh() != NULL);
    assert(pExtracellularStimulusFactory->GetNumberOfCells() == pCellFactorySecondCell->GetNumberOfCells() );

    unsigned num_local_nodes = this->mpDistributedVectorFactory->GetLocalOwnership();
    unsigned ownership_range_low = this->mpDistributedVectorFactory->GetLow();
    mCellsDistributedSecondCell.resize(num_local_nodes);
    mGgapDistributed.resize(num_local_nodes);
    mExtracellularStimuliDistributed.resize(num_local_nodes);

    try
    {
        for (unsigned local_index = 0; local_index < num_local_nodes; local_index++)
        {
            unsigned global_index = ownership_range_low + local_index;
            mCellsDistributedSecondCell[local_index] = pCellFactorySecondCell->CreateCardiacCellForNode(global_index);
            mCellsDistributedSecondCell[local_index]->SetUsedInTissueSimulation();
            mGgapDistributed[local_index] = 0.0;//default. It will be changed by specific method later when user input will be obtained
            mExtracellularStimuliDistributed[local_index] = pExtracellularStimulusFactory->CreateStimulusForNode(global_index);
        }

        pCellFactorySecondCell->FinaliseCellCreation(&mCellsDistributedSecondCell,
                                           this->mpDistributedVectorFactory->GetLow(),
                                           this->mpDistributedVectorFactory->GetHigh());
    }
    catch (const Exception& e)
    {
#define COVERAGE_IGNORE //don't really know how to cover this...
        // Errors thrown creating cells will often be process-specific
        PetscTools::ReplicateException(true);
        // Should really do this for other processes too, but this is all we need
        // to get memory testing to pass, and leaking when we're about to die isn't
        // that bad! Delete second cells
        for (std::vector<AbstractCardiacCell*>::iterator cell_iterator = mCellsDistributedSecondCell.begin();
             cell_iterator != mCellsDistributedSecondCell.end();
             ++cell_iterator)
        {
            delete (*cell_iterator);
        }
        throw e;
#undef COVERAGE_IGNORE
    }
    PetscTools::ReplicateException(false);

    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    mIionicCacheReplicatedSecondCell.Resize( pCellFactory->GetNumberOfCells() );
    mIntracellularStimulusCacheReplicatedSecondCell.Resize( pCellFactorySecondCell->GetNumberOfCells() );
    mGgapCacheReplicated.Resize(pCellFactorySecondCell->GetNumberOfCells());//this is a bit of a hack...
    mExtracellularStimulusCacheReplicated.Resize(pExtracellularStimulusFactory->GetNumberOfCells());
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);

    //Creat the extracellular conductivity tensor
    CreateExtracellularConductivityTensors();
}

//archiving constructor
template <unsigned SPACE_DIM>
ExtendedBidomainTissue<SPACE_DIM>::ExtendedBidomainTissue(std::vector<AbstractCardiacCell*> & rCellsDistributed, std::vector<AbstractCardiacCell*> & rSecondCellsDistributed, std::vector<boost::shared_ptr<AbstractStimulusFunction> > & rExtraStimuliDistributed, std::vector<double> & rGgapsDistributed, AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh,c_vector<double, SPACE_DIM>  intracellularConductivitiesSecondCell)
        :   AbstractCardiacTissue<SPACE_DIM>(pMesh),
			mpIntracellularConductivityTensorsSecondCell(NULL),
			mIntracellularConductivitiesSecondCell(intracellularConductivitiesSecondCell),
			mCellsDistributedSecondCell(rSecondCellsDistributed),
			mExtracellularStimuliDistributed(rExtraStimuliDistributed),
			mGgapDistributed(rGgapsDistributed),
			mUserSuppliedExtracellularStimulus(false)
{
	//segfault guards in case we failed to load anything from the archive
	assert(mCellsDistributedSecondCell.size() > 0);
	assert(mExtracellularStimuliDistributed.size() > 0);
    assert(mGgapDistributed.size() > 0);
	//allocate memory for the caches
    mIionicCacheReplicatedSecondCell.Resize( this->mpDistributedVectorFactory->GetProblemSize());
    mIntracellularStimulusCacheReplicatedSecondCell.Resize( this->mpDistributedVectorFactory->GetProblemSize() );
    mGgapCacheReplicated.Resize(this->mpDistributedVectorFactory->GetProblemSize());
    mExtracellularStimulusCacheReplicated.Resize(this->mpDistributedVectorFactory->GetProblemSize());

    CreateIntracellularConductivityTensorSecondCell();
    CreateExtracellularConductivityTensors();
}


template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::SetGgapHeterogeneities ( std::vector<boost::shared_ptr<AbstractChasteRegion<SPACE_DIM> > >& rGgapHeterogeneityRegions, std::vector<double> rGgapValues)
{
    assert( rGgapHeterogeneityRegions.size() == rGgapValues.size() );//problem class (which calls this method should have thrown otherwise)
    mGgapHeterogeneityRegions = rGgapHeterogeneityRegions;
    mGgapValues =rGgapValues;
}

template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::CreateGGapConductivities()
{
    assert(mGgapHeterogeneityRegions.size() == mGgapValues.size());
    assert(this->mpMesh != NULL);
    
    unsigned num_local_nodes = this->mpDistributedVectorFactory->GetLocalOwnership();
    unsigned ownership_range_low = this->mpDistributedVectorFactory->GetLow();
    assert(mGgapDistributed.size() == num_local_nodes);//the constructor should have allocated memory.
    try
    {
        for (unsigned local_index = 0; local_index < num_local_nodes; local_index++)
        {
			mGgapDistributed[local_index] = mGGap;//assign default uniform value everywhere first

			//then change where and if necessary
			unsigned global_index = ownership_range_low + local_index;
			for (unsigned het_index = 0; het_index < mGgapHeterogeneityRegions.size(); het_index++)
			{
				if ( mGgapHeterogeneityRegions[het_index]->DoesContain ( this->mpMesh->GetNode(global_index)->GetPoint() ) )
				{
					mGgapDistributed[local_index] = mGgapValues[het_index];
				}
			}
        }
    }
    catch (const Exception& e)
    {
#define COVERAGE_IGNORE 
        PetscTools::ReplicateException(true);
        throw e;
#undef COVERAGE_IGNORE
    }
    PetscTools::ReplicateException(false);
}

template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::CreateIntracellularConductivityTensorSecondCell()
{
	HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
    this->mpConfig = HeartConfig::Instance();
    mpIntracellularConductivityTensorsSecondCell = new OrthotropicConductivityTensors<SPACE_DIM,SPACE_DIM>;

    // this definition must be here (and not inside the if statement) because SetNonConstantConductivities() will keep
    // a pointer to it and we don't want it to go out of scope before Init() is called
    unsigned num_elements = this->mpMesh->GetNumElements();
    std::vector<c_vector<double, SPACE_DIM> > hetero_intra_conductivities;

    c_vector<double, SPACE_DIM> intra_conductivities;
    this->mpConfig->GetIntracellularConductivities(intra_conductivities);//this one is used just for resizing

    if (this->mpConfig->GetConductivityHeterogeneitiesProvided())
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

        std::vector<boost::shared_ptr<AbstractChasteRegion<SPACE_DIM> > > conductivities_heterogeneity_areas;
        std::vector< c_vector<double,3> > intra_h_conductivities;
        std::vector< c_vector<double,3> > extra_h_conductivities;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                intra_h_conductivities,
                                                                extra_h_conductivities);
        unsigned local_element_index = 0;
        for (typename AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>::ElementIterator it = this->mpMesh->GetElementIteratorBegin();
             it != this->mpMesh->GetElementIteratorEnd();
             ++it)
        {
            //unsigned element_index = it->GetIndex();
            // if element centroid is contained in the region
            ChastePoint<SPACE_DIM> element_centroid(it->CalculateCentroid());
            for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas.size(); region_index++)
            {
                if ( conductivities_heterogeneity_areas[region_index]->DoesContain(element_centroid) )
                {
                	//We don't use ublas vector assignment here, because we might be getting a subvector of a 3-vector
                	for (unsigned i=0; i<SPACE_DIM; i++)
                	{
                		hetero_intra_conductivities[local_element_index][i] = intra_h_conductivities[region_index][i];
                	}
                }
            }
            local_element_index++;
        }

        mpIntracellularConductivityTensorsSecondCell->SetNonConstantConductivities(&hetero_intra_conductivities);
    }
    else
    {
    	mpIntracellularConductivityTensorsSecondCell->SetConstantConductivities(mIntracellularConductivitiesSecondCell);
    }

    mpIntracellularConductivityTensorsSecondCell->Init(this->mpMesh);
    HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);
}

template <unsigned SPACE_DIM>
bool ExtendedBidomainTissue<SPACE_DIM>::HasTheUserSuppliedExtracellularStimulus()
{
	return mUserSuppliedExtracellularStimulus;
}

template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::SetUserSuppliedExtracellularStimulus(bool flag)
{
	mUserSuppliedExtracellularStimulus = flag;
}

template <unsigned SPACE_DIM>
const std::vector<AbstractCardiacCell*>& ExtendedBidomainTissue<SPACE_DIM>::rGetSecondCellsDistributed() const
{
    return mCellsDistributedSecondCell;
}

template <unsigned SPACE_DIM>
const std::vector<double>& ExtendedBidomainTissue<SPACE_DIM>::rGetGapsDistributed() const
{
    return mGgapDistributed;
}

template <unsigned SPACE_DIM>
const std::vector<boost::shared_ptr<AbstractStimulusFunction> >& ExtendedBidomainTissue<SPACE_DIM>::rGetExtracellularStimulusDistributed() const
{
    return mExtracellularStimuliDistributed;
}


template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::CreateExtracellularConductivityTensors()
{

    mpExtracellularConductivityTensors =  new OrthotropicConductivityTensors<SPACE_DIM,SPACE_DIM>;
    c_vector<double, SPACE_DIM> extra_conductivities;
    this->mpConfig->GetExtracellularConductivities(extra_conductivities);

    // this definition must be here (and not inside the if statement) because SetNonConstantConductivities() will keep
    // a pointer to it and we don't want it to go out of scope before Init() is called
    unsigned num_elements = this->mpMesh->GetNumElements();
    std::vector<c_vector<double, SPACE_DIM> > hetero_extra_conductivities;

    if (this->mpConfig->GetConductivityHeterogeneitiesProvided())
    {
        try
        {
            assert(hetero_extra_conductivities.size()==0);
            //initialise with the values of teh default conductivity tensor
            hetero_extra_conductivities.resize(num_elements, extra_conductivities);
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

        std::vector<boost::shared_ptr<AbstractChasteRegion<SPACE_DIM> > > conductivities_heterogeneity_areas;
        std::vector< c_vector<double,3> > intra_h_conductivities;
        std::vector< c_vector<double,3> > extra_h_conductivities;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                intra_h_conductivities,
                                                                extra_h_conductivities);
        unsigned local_element_index = 0;
        for (typename AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>::ElementIterator iter = (this->mpMesh)->GetElementIteratorBegin();
             iter != (this->mpMesh)->GetElementIteratorEnd();
             ++iter)
        {
            //unsigned element_index = iter->GetIndex();
            // if element centroid is contained in the region
            ChastePoint<SPACE_DIM> element_centroid(iter->CalculateCentroid());
            for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas.size(); region_index++)
            {
                // if element centroid is contained in the region
                if ( conductivities_heterogeneity_areas[region_index]->DoesContain( element_centroid ) )
                {
                	//We don't use ublas vector assignment here, because we might be getting a subvector of a 3-vector
                	for (unsigned i=0; i<SPACE_DIM; i++)
                	{
                		hetero_extra_conductivities[local_element_index][i] = extra_h_conductivities[region_index][i];
                	}
                }
            }
            local_element_index++;
        }

        mpExtracellularConductivityTensors->SetNonConstantConductivities(&hetero_extra_conductivities);
    }
    else
    {
        mpExtracellularConductivityTensors->SetConstantConductivities(extra_conductivities);
    }
    mpExtracellularConductivityTensors->Init(this->mpMesh);
}

template <unsigned SPACE_DIM>
ExtendedBidomainTissue<SPACE_DIM>::~ExtendedBidomainTissue()
{
    // Delete (second) cells
    for (std::vector<AbstractCardiacCell*>::iterator cell_iterator = mCellsDistributedSecondCell.begin();
         cell_iterator != mCellsDistributedSecondCell.end();
         ++cell_iterator)
    {
        delete (*cell_iterator);
    }

    if (mpExtracellularConductivityTensors)
    {
        delete mpExtracellularConductivityTensors;
    }

    if (mpIntracellularConductivityTensorsSecondCell)
    {
    	delete mpIntracellularConductivityTensorsSecondCell;
    }
}

template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::SetIntracellularConductivitiesSecondCell(c_vector<double, SPACE_DIM> conductivities)
{
	for (unsigned i = 0; i < SPACE_DIM; i++)
	{
		mIntracellularConductivitiesSecondCell[i] = conductivities[i];
	}
}

template <unsigned SPACE_DIM>
c_vector<double, SPACE_DIM>  ExtendedBidomainTissue<SPACE_DIM>::GetIntracellularConductivitiesSecondCell() const
{
	return mIntracellularConductivitiesSecondCell;
}

template <unsigned SPACE_DIM>
const c_matrix<double, SPACE_DIM, SPACE_DIM>& ExtendedBidomainTissue<SPACE_DIM>::rGetExtracellularConductivityTensor(unsigned elementIndex)
{
    assert(mpExtracellularConductivityTensors);
    return (*mpExtracellularConductivityTensors)[elementIndex];
}

template <unsigned SPACE_DIM>
const c_matrix<double, SPACE_DIM, SPACE_DIM>& ExtendedBidomainTissue<SPACE_DIM>::rGetIntracellularConductivityTensorSecondCell(unsigned elementIndex)
{
	assert(mpIntracellularConductivityTensorsSecondCell);
	return (*mpIntracellularConductivityTensorsSecondCell)[elementIndex];
}

template <unsigned SPACE_DIM>
AbstractCardiacCell* ExtendedBidomainTissue<SPACE_DIM>::GetCardiacSecondCell( unsigned globalIndex )
{
    return mCellsDistributedSecondCell[globalIndex - this->mpDistributedVectorFactory->GetLow()];
}

template <unsigned SPACE_DIM>
boost::shared_ptr<AbstractStimulusFunction> ExtendedBidomainTissue<SPACE_DIM>::GetExtracellularStimulus( unsigned globalIndex )
{
    return mExtracellularStimuliDistributed[globalIndex - this->mpDistributedVectorFactory->GetLow()];
}

template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::SolveCellSystems(Vec existingSolution, double time, double nextTime, bool updateVoltage)
{
    HeartEventHandler::BeginEvent(HeartEventHandler::SOLVE_ODES);

    DistributedVector dist_solution = this->mpDistributedVectorFactory->CreateDistributedVector(existingSolution);
    DistributedVector::Stripe phi_i_first_cell(dist_solution, 0);
    DistributedVector::Stripe phi_i_second_cell(dist_solution, 1);
    DistributedVector::Stripe phi_e(dist_solution, 2);

    for (DistributedVector::Iterator index = dist_solution.Begin();
         index != dist_solution.End();
         ++index)
    {
    	double voltage_first_cell = phi_i_first_cell[index] - phi_e[index];
    	double voltage_second_cell = phi_i_second_cell[index] - phi_e[index];

        // overwrite the voltage with the input value
        this->mCellsDistributed[index.Local]->SetVoltage( voltage_first_cell );
        mCellsDistributedSecondCell[index.Local]->SetVoltage( voltage_second_cell );
        try
        {
            // solve
            // Note: Voltage should not be updated. GetIIonic will be called later
            // and needs the old voltage. The voltage will be updated from the pde.
            this->mCellsDistributed[index.Local]->ComputeExceptVoltage(time, nextTime);
            mCellsDistributedSecondCell[index.Local]->ComputeExceptVoltage(time, nextTime);
        }
        catch (Exception &e)
        {
#define COVERAGE_IGNORE
            PetscTools::ReplicateException(true);
            throw e;
#undef COVERAGE_IGNORE
        }

        // update the Iionic and stimulus caches
        this->UpdateCaches(index.Global, index.Local, nextTime);//in parent class
        UpdateAdditionalCaches(index.Global, index.Local, nextTime);//extended bidomain specific caches
    }
    PetscTools::ReplicateException(false);
    HeartEventHandler::EndEvent(HeartEventHandler::SOLVE_ODES);

    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    if ( this->mDoCacheReplication )
    {
        this->ReplicateCaches();
        ReplicateAdditionalCaches();//extended bidomain specific caches
    }
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
}

template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::UpdateAdditionalCaches(unsigned globalIndex, unsigned localIndex, double nextTime)
{
	mIionicCacheReplicatedSecondCell[globalIndex] = mCellsDistributedSecondCell[localIndex]->GetIIonic();
	mIntracellularStimulusCacheReplicatedSecondCell[globalIndex] = mCellsDistributedSecondCell[localIndex]->GetIntracellularStimulus(nextTime);
	mExtracellularStimulusCacheReplicated[globalIndex] = mExtracellularStimuliDistributed[localIndex]->GetStimulus(nextTime);
	mGgapCacheReplicated[globalIndex] = mGgapDistributed[localIndex];
}

template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::ReplicateAdditionalCaches()
{
	mIionicCacheReplicatedSecondCell.Replicate(this->mpDistributedVectorFactory->GetLow(), this->mpDistributedVectorFactory->GetHigh());
	mIntracellularStimulusCacheReplicatedSecondCell.Replicate(this->mpDistributedVectorFactory->GetLow(), this->mpDistributedVectorFactory->GetHigh());
	mExtracellularStimulusCacheReplicated.Replicate(this->mpDistributedVectorFactory->GetLow(), this->mpDistributedVectorFactory->GetHigh());
	mGgapCacheReplicated.Replicate(this->mpDistributedVectorFactory->GetLow(), this->mpDistributedVectorFactory->GetHigh());
}

template <unsigned SPACE_DIM>
ReplicatableVector& ExtendedBidomainTissue<SPACE_DIM>::rGetIionicCacheReplicatedSecondCell()
{
    return mIionicCacheReplicatedSecondCell;
}

template <unsigned SPACE_DIM>
ReplicatableVector& ExtendedBidomainTissue<SPACE_DIM>::rGetIntracellularStimulusCacheReplicatedSecondCell()
{
    return mIntracellularStimulusCacheReplicatedSecondCell;
}

template <unsigned SPACE_DIM>
ReplicatableVector& ExtendedBidomainTissue<SPACE_DIM>::rGetExtracellularStimulusCacheReplicated()
{
    return mExtracellularStimulusCacheReplicated;
}

template <unsigned SPACE_DIM>
ReplicatableVector& ExtendedBidomainTissue<SPACE_DIM>::rGetGgapCacheReplicated()
{
    return mGgapCacheReplicated;
}
     
template <unsigned SPACE_DIM>
double ExtendedBidomainTissue<SPACE_DIM>::GetAmFirstCell()
{
	return mAmFirstCell;
}

template <unsigned SPACE_DIM>
double ExtendedBidomainTissue<SPACE_DIM>::GetAmSecondCell()
{
	return mAmSecondCell;
}

template <unsigned SPACE_DIM>
double ExtendedBidomainTissue<SPACE_DIM>::GetAmGap()
{
	return mAmGap;
}

template <unsigned SPACE_DIM>
double ExtendedBidomainTissue<SPACE_DIM>::GetCmFirstCell()
{
	return mCmFirstCell;
}

template <unsigned SPACE_DIM>
double ExtendedBidomainTissue<SPACE_DIM>::GetCmSecondCell()
{
	return mCmSecondCell;
}

template <unsigned SPACE_DIM>
double ExtendedBidomainTissue<SPACE_DIM>::GetGGap()
{
	return mGGap;
}

template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::SetAmFirstCell(double value)
{
	mAmFirstCell = value;
}

template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::SetAmSecondCell(double value)
{
	mAmSecondCell = value;
}

template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::SetAmGap(double value)
{
	mAmGap = value;
}

template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::SetGGap(double value)
{
	mGGap = value;
}

template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::SetCmFirstCell(double value)
{
	mCmFirstCell = value;
}

template <unsigned SPACE_DIM>
void ExtendedBidomainTissue<SPACE_DIM>::SetCmSecondCell(double value)
{
	mCmSecondCell = value;
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class ExtendedBidomainTissue<1>;
template class ExtendedBidomainTissue<2>;
template class ExtendedBidomainTissue<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExtendedBidomainTissue)