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

#include "BidomainPde.hpp"

#include "DistributedVector.hpp"
#include "AxisymmetricConductivityTensors.hpp"
#include "OrthotropicConductivityTensors.hpp"
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"


template <unsigned SPACE_DIM>
BidomainPde<SPACE_DIM>::BidomainPde(
            AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
    : AbstractCardiacPde<SPACE_DIM>(pCellFactory, 2 /*mStride*/)
{
    mExtracellularStimulusCacheReplicated.resize( pCellFactory->GetNumberOfCells() );
    CreateExtracellularConductivityTensors();
}

template <unsigned SPACE_DIM>
BidomainPde<SPACE_DIM>::BidomainPde(std::vector<AbstractCardiacCell*> &rCellsDistributed,AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh)
        :  AbstractCardiacPde<SPACE_DIM>(rCellsDistributed, pMesh, 2u) // The 2 tells it this is a bidomain
{
    mExtracellularStimulusCacheReplicated.resize(this->mpDistributedVectorFactory->GetProblemSize());
    CreateExtracellularConductivityTensors();
}

template <unsigned SPACE_DIM>
void BidomainPde<SPACE_DIM>::CreateExtracellularConductivityTensors()
{
    if (this->mpConfig->IsMeshProvided() && this->mpConfig->GetLoadMesh())
    {
        switch (this->mpConfig->GetConductivityMedia())
        {
            case media_type::Orthotropic:
                mpExtracellularConductivityTensors =  new OrthotropicConductivityTensors<SPACE_DIM>;
                mpExtracellularConductivityTensors->SetFibreOrientationFile(this->mpConfig->GetMeshName() + ".ortho");
                break;

            case media_type::Axisymmetric:
                mpExtracellularConductivityTensors =  new AxisymmetricConductivityTensors<SPACE_DIM>;
                mpExtracellularConductivityTensors->SetFibreOrientationFile(this->mpConfig->GetMeshName() + ".axi");
                break;

            case media_type::NoFibreOrientation:
                mpExtracellularConductivityTensors =  new OrthotropicConductivityTensors<SPACE_DIM>;
                break;

            default :
                NEVER_REACHED;
        }
    }
    else // Slab defined in config file or SetMesh() called; no fibre orientation assumed
    {
        mpExtracellularConductivityTensors =  new OrthotropicConductivityTensors<SPACE_DIM>;
    }

    c_vector<double, SPACE_DIM> extra_conductivities;
    this->mpConfig->GetExtracellularConductivities(extra_conductivities);

    // this definition must be here (and not inside the if statement) because SetNonConstantConductivities() will keep
    // a pointer to it and we don't want it to go out of scope before Init() is called
    unsigned num_elements = this->mpMesh->GetNumElements();
    std::vector<c_vector<double, SPACE_DIM> > hetero_extra_conductivities(num_elements);

    if (this->mpConfig->GetConductivityHeterogeneitiesProvided())
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
                ChastePoint<SPACE_DIM> element_centroid(this->mpMesh->GetElement(element_index)->CalculateCentroid());
                if ( conductivities_heterogeneity_areas[region_index].DoesContain( element_centroid ) )
                {
                    hetero_extra_conductivities[element_index] = extra_h_conductivities[region_index];
                }
                else
                {
                    hetero_extra_conductivities[element_index] = extra_conductivities;
                }
            }
        }

        mpExtracellularConductivityTensors->SetNonConstantConductivities(&hetero_extra_conductivities);
    }
    else
    {
        mpExtracellularConductivityTensors->SetConstantConductivities(extra_conductivities);
    }

    mpExtracellularConductivityTensors->Init();
}

template <unsigned SPACE_DIM>
BidomainPde<SPACE_DIM>::~BidomainPde()
{
    /// \todo: #98 once the archiving constructor is creating the conductivity tensors properly we won't need this if statement
    if (mpExtracellularConductivityTensors)
    {
        delete mpExtracellularConductivityTensors;
    }
}


template <unsigned SPACE_DIM>
const c_matrix<double, SPACE_DIM, SPACE_DIM>& BidomainPde<SPACE_DIM>::rGetExtracellularConductivityTensor(unsigned elementIndex)
{
    assert(mpExtracellularConductivityTensors);
    return (*mpExtracellularConductivityTensors)[elementIndex];
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class BidomainPde<1>;
template class BidomainPde<2>;
template class BidomainPde<3>;
