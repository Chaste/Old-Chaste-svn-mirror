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

#include "HeartConfigRelatedCellFactory.hpp"


template<unsigned SPACE_DIM>
HeartConfigRelatedCellFactory<SPACE_DIM>::HeartConfigRelatedCellFactory()
    : AbstractCardiacCellFactory<SPACE_DIM>(),
      mDefaultIonicModel(HeartConfig::Instance()->GetDefaultIonicModel())
{
    // Read and store possible region definitions
    HeartConfig::Instance()->GetIonicModelRegions(mIonicModelRegions,
                                                  mIonicModelsDefined);

    // Read and store Stimuli
    try
    {
        HeartConfig::Instance()->GetStimuli(mStimuliApplied, mStimulatedAreas);
    }
    catch (Exception& e)
    {
        // No stimuli provided in XML so we should have hit a parsing exception
        NEVER_REACHED;
    }

    // Read and store Cell Heterogeneities
    try
    {
        HeartConfig::Instance()->GetCellHeterogeneities(mCellHeterogeneityAreas,
                                                        mScaleFactorGks,
                                                        mScaleFactorIto,
                                                        mScaleFactorGkr);
    }
    catch (Exception& e)
    {
        // No cell heterogeneities provided
    }
}


template<unsigned SPACE_DIM>
AbstractCardiacCell* HeartConfigRelatedCellFactory<SPACE_DIM>::CreateCellWithIntracellularStimulus(
        boost::shared_ptr<AbstractStimulusFunction> intracellularStimulus,
        unsigned nodeIndex)
{
    cp::ionic_models_available_type ionic_model = mDefaultIonicModel;

    for (unsigned ionic_model_region_index = 0;
         ionic_model_region_index < mIonicModelRegions.size();
         ++ionic_model_region_index)
    {
        if ( mIonicModelRegions[ionic_model_region_index].DoesContain(this->GetMesh()->GetNode(nodeIndex)->GetPoint()) )
        {
            ionic_model = mIonicModelsDefined[ionic_model_region_index];
            break;
        }
    }

    switch(ionic_model)
    {
        case(cp::ionic_models_available_type::LuoRudyI):
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver, intracellularStimulus);
            break;
        }

        case(cp::ionic_models_available_type::LuoRudyIBackwardEuler):
        {
            return new BackwardEulerLuoRudyIModel1991(intracellularStimulus);
            break;
        }

        case(cp::ionic_models_available_type::Fox2002BackwardEuler):
        {
            return new BackwardEulerFoxModel2002Modified(intracellularStimulus);
            break;
        }

        case(cp::ionic_models_available_type::DifrancescoNoble):
        {
            return new DiFrancescoNoble1985OdeSystem(this->mpSolver, intracellularStimulus);
            break;
        }

        case(cp::ionic_models_available_type::MahajanShiferaw):
        {
            return new Mahajan2008OdeSystem(this->mpSolver, intracellularStimulus);
            break;
        }

        case(cp::ionic_models_available_type::tenTusscher2006):
        {
            TenTusscher2006OdeSystem*  p_tt06_instance = new TenTusscher2006OdeSystem(this->mpSolver, intracellularStimulus);

            for (unsigned ht_index = 0;
                 ht_index < mCellHeterogeneityAreas.size();
                 ++ht_index)
            {
                if ( mCellHeterogeneityAreas[ht_index].DoesContain(this->GetMesh()->GetNode(nodeIndex)->GetPoint()) )
                {
                    p_tt06_instance->SetScaleFactorGks(mScaleFactorGks[ht_index]);
                    p_tt06_instance->SetScaleFactorIto(mScaleFactorIto[ht_index]);
                    p_tt06_instance->SetScaleFactorGkr(mScaleFactorGkr[ht_index]);
                }
            }

            return p_tt06_instance;
            break;
        }
        
        case(cp::ionic_models_available_type::Maleckar):
        {
             Maleckar2009OdeSystem*  p_maleckar_instance = new Maleckar2009OdeSystem(this->mpSolver, intracellularStimulus);

            for (unsigned ht_index = 0;
                 ht_index < mCellHeterogeneityAreas.size();
                 ++ht_index)
            {
                if ( mCellHeterogeneityAreas[ht_index].DoesContain(this->GetMesh()->GetNode(nodeIndex)->GetPoint()) )
                {
                    p_maleckar_instance->SetScaleFactorGks(mScaleFactorGks[ht_index]);
                    p_maleckar_instance->SetScaleFactorIto(mScaleFactorIto[ht_index]);
                    p_maleckar_instance->SetScaleFactorGkr(mScaleFactorGkr[ht_index]);
                }
            }

            return p_maleckar_instance;
            break;
        }
        
        case(cp::ionic_models_available_type::HodgkinHuxley):
        {
            return new HodgkinHuxleySquidAxon1952OriginalOdeSystem(this->mpSolver, intracellularStimulus);
            break;
        }

        case(cp::ionic_models_available_type::FaberRudy2000):
        {
            FaberRudy2000Version3*  p_faber_rudy_instance = new FaberRudy2000Version3(this->mpSolver, intracellularStimulus);

            for (unsigned ht_index = 0;
                 ht_index < mCellHeterogeneityAreas.size();
                 ++ht_index)
            {
                if ( mCellHeterogeneityAreas[ht_index].DoesContain(this->GetMesh()->GetNode(nodeIndex)->GetPoint()) )
                {
                    p_faber_rudy_instance->SetScaleFactorGks(mScaleFactorGks[ht_index]);
                    p_faber_rudy_instance->SetScaleFactorIto(mScaleFactorIto[ht_index]);
                    p_faber_rudy_instance->SetScaleFactorGkr(mScaleFactorGkr[ht_index]);
                }
            }

            return p_faber_rudy_instance;
            break;
        }

        case(cp::ionic_models_available_type::FaberRudy2000Optimised):
        {
            return new FaberRudy2000Version3Optimised(this->mpSolver, intracellularStimulus);
            break;
        }

        default:
        {
           //If the ionic model is not in the current enumeration then the XML parser will have picked it up before now!
           NEVER_REACHED;
        }
    }

    return NULL;
}


template<unsigned SPACE_DIM>
AbstractCardiacCell* HeartConfigRelatedCellFactory<SPACE_DIM>::CreateCardiacCellForTissueNode(unsigned nodeIndex)
{
    boost::shared_ptr<MultiStimulus> node_specific_stimulus(new MultiStimulus());

    // Check which of the defined stimuli contain the current node
    for (unsigned stimulus_index = 0;
         stimulus_index < mStimuliApplied.size();
         ++stimulus_index)
    {
        if ( mStimulatedAreas[stimulus_index].DoesContain(this->GetMesh()->GetNode(nodeIndex)->GetPoint()) )
        {
            node_specific_stimulus->AddStimulus(mStimuliApplied[stimulus_index]);
        }
    }

    return CreateCellWithIntracellularStimulus(node_specific_stimulus, nodeIndex);
}


// Explicit instantiation
template class HeartConfigRelatedCellFactory<1u>;
template class HeartConfigRelatedCellFactory<2u>;
template class HeartConfigRelatedCellFactory<3u>;
