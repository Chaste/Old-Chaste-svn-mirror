/*

Copyright (C) University of Oxford, 2008

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


#ifndef BIDOMAINPDE_HPP_
#define BIDOMAINPDE_HPP_

#include <vector>
#include "Node.hpp"
#include "AbstractStimulusFunction.hpp"
#include "AbstractCardiacPde.hpp"
#include <boost/numeric/ublas/matrix.hpp>


/**
 * BidomainPde class.
 *
 * The bidmain equation is of the form:
 *
 * A_m ( C_m d(V_m)/dt + I_ion ) = div ( sigma_i grad( V_m + phi_e ) ) + I_intra_stim
 *   and
 * div ( (sigma_i + sigma_e) grad phi_e    +   sigma_i (grad V_m) )  + I_extra_stim
 *
 * where V_m is the trans-membrane potential = phi_i - phi_e            (mV),
 *       phi_i is the intracellular potential                           (mV),
 *       phi_e is the intracellular potential                           (mV),
 * and   A_m is the surface area to volume ratio of the cell membrane   (1/cm),
 *       C_m is the membrane capacitance                                (uF/cm^2),
 *       sigma_i is the intracellular conductivity tensor               (mS/cm),
 *       sigma_e is the intracellular conductivity tensor               (mS/cm),
 * and   I_ion is the ionic current                                     (uA/cm^2),
 *       I_intra_stim is the internal stimulus                          (uA/cm^3),
 *       I_extra_stim is the external stimulus (a shock)                (uA/cm^3).
 */
template <unsigned SPACE_DIM>
class BidomainPde : public virtual AbstractCardiacPde<SPACE_DIM>
{

private:
    AbstractConductivityTensors<SPACE_DIM> *mpExtracellularConductivityTensors;

    ReplicatableVector mExtracellularStimulusCacheReplicated;

public:
    //Constructor
    BidomainPde(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            :  AbstractCardiacPde<SPACE_DIM>(pCellFactory, 2 /*mStride*/)
    {
        mExtracellularStimulusCacheReplicated.resize( pCellFactory->GetNumberOfCells() );

        if (this->mpConfig->GetIsMeshProvided() && this->mpConfig->GetLoadMesh())
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
		std::vector<c_vector<double, SPACE_DIM> > hetero_extra_conductivities; 

		if (this->mpConfig->GetConductivityHeterogeneitiesProvided())
		{
			unsigned num_elements = pCellFactory->GetMesh()->GetNumElements();
			hetero_extra_conductivities.reserve(num_elements);						

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
    
    ~BidomainPde()
    {
        delete mpExtracellularConductivityTensors;
    }

//    void SetExtracellularConductivityTensors(AbstractConductivityTensors<SPACE_DIM>* pExtracellularTensors)
//    {
//        mpExtracellularConductivityTensors = pExtracellularTensors;
//    }


    const c_matrix<double, SPACE_DIM, SPACE_DIM>& rGetExtracellularConductivityTensor(unsigned elementIndex)
    {
        assert(mpExtracellularConductivityTensors);
        return (*mpExtracellularConductivityTensors)[elementIndex];
    }

    /**
     * The bidomain pde also updates the extracellular stimulus cache
     */
    void UpdateCaches(unsigned globalIndex, unsigned localIndex, double nextTime)
    {
        AbstractCardiacPde<SPACE_DIM>::UpdateCaches(globalIndex, localIndex, nextTime);
        mExtracellularStimulusCacheReplicated[globalIndex] = this->mCellsDistributed[localIndex]->GetExtracellularStimulus(nextTime);
    }

    /**
     * The bidomain Pde also replicates the extracellular stimulus cache
     */
    void ReplicateCaches()
    {
        AbstractCardiacPde<SPACE_DIM>::ReplicateCaches();
        unsigned lo=DistributedVector::Begin().Global;
        unsigned hi=DistributedVector::End().Global;

        mExtracellularStimulusCacheReplicated.Replicate(lo, hi);
    }

    ReplicatableVector& rGetExtracellularStimulusCacheReplicated()
    {
        return mExtracellularStimulusCacheReplicated;
    }
};



#endif /*BIDOMAINPDE_HPP_*/
