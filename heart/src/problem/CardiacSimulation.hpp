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

#ifndef CARDIACSIMULATION_HPP_
#define CARDIACSIMULATION_HPP_

// Must go first
#include "CardiacSimulationArchiver.hpp"

#include <vector>
#include <ctime>
#include <memory>

#include "UblasIncludes.hpp"

#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include "PetscTools.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "TimeStepper.hpp"

#include "HeartConfig.hpp"

#include "AbstractStimulusFunction.hpp"
#include "MultiStimulus.hpp"
#include "SimpleStimulus.hpp"

#include "TetrahedralMesh.hpp"
#include "NonCachedTetrahedralMesh.hpp"
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "TrianglesMeshWriter.hpp"

#include "BackwardEulerFoxModel2002Modified.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BackwardEulerLuoRudyIModel1991.hpp"
#include "FoxModel2002Modified.hpp"
#include "FaberRudy2000Version3.hpp"
#include "FaberRudy2000Version3Optimised.hpp"
#include "DiFrancescoNoble1985OdeSystem.hpp"
#include "Mahajan2008OdeSystem.hpp"
#include "TenTusscher2006OdeSystem.hpp"
#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"

#include "OrthotropicConductivityTensors.hpp"

#include "Hdf5ToMeshalyzerConverter.hpp"
#include "PostProcessingWriter.hpp"
#include "Version.hpp"

/**
 *  \todo: dox, coverage, maybe own tests
 */
template<unsigned SPACE_DIM>
class HeartConfigRelatedCellFactory : public AbstractCardiacCellFactory<SPACE_DIM>
{
private:
    /** Default cardiac cell model to be used in all tissue (except heterogeneous regions)*/
    cp::ionic_models_available_type mDefaultIonicModel;
    /** List of axis-aligned box regions which contain heterogeneous cardiac ionic model types*/
    std::vector<ChasteCuboid> mIonicModelRegions;
    /** List of ionic model (size matches that of mIonicModelRegions)*/
    std::vector<cp::ionic_models_available_type> mIonicModelsDefined;
    
    /** List of axis-aligned box regions which represent areas to stimulate*/
    std::vector<ChasteCuboid> mStimulatedAreas;
    /** List of intracellular current stimuli to apply (size matches that of mStimulatedAreas)*/
    std::vector<boost::shared_ptr<SimpleStimulus> > mStimuliApplied;
    
    /** List of axis-aligned box regions which represent areas in which to give parametric heterogeneity (scaling gating parameters)*/
    std::vector<ChasteCuboid> mCellHeterogeneityAreas;
    /** List of scale factors for Gks scaling in each region (size of list matches that of mCellHeterogeneityAreas)*/
    std::vector<double> mScaleFactorGks;
    /** List of scale factors for Ito scaling in each region (size of list matches that of mCellHeterogeneityAreas)*/
    std::vector<double> mScaleFactorIto;
    /** List of scale factors for Gkr scaling in each region (size of list matches that of mCellHeterogeneityAreas)*/
    std::vector<double> mScaleFactorGkr;

public:
    HeartConfigRelatedCellFactory()
        : AbstractCardiacCellFactory<SPACE_DIM>(),
          mDefaultIonicModel(cp::ionic_models_available_type::LuoRudyI)
    {        
        // Read and store default ionic model and possible region definitions
        mDefaultIonicModel = HeartConfig::Instance()->GetDefaultIonicModel();
        HeartConfig::Instance()->GetIonicModelRegions(mIonicModelRegions,
                                                      mIonicModelsDefined);
    
        // Read and store Stimuli
        try
        {
            HeartConfig::Instance()->GetStimuli(mStimuliApplied, mStimulatedAreas);
        }
        catch(Exception& e)
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
        catch(Exception& e)
        {
            // No cell heterogeneities provided
        }
    }


    /**
     * Create the correct tissue cell for a given region in the mesh
     * @param intracellularStimulus is computed in CreateCardiacCellForTissueNode determined by the list of stimulation regions
     * @param nodeIndex is the global index within the mesh
     */
    AbstractCardiacCell* CreateCellWithIntracellularStimulus(boost::shared_ptr<AbstractStimulusFunction> intracellularStimulus, unsigned nodeIndex)
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


    /**
     * Create the correct stimulated tissue cell for a given region in the mesh
     * The stimulus is determined in this method (using the list of stimulation regions).  
     * The cardiac cell type (and parameters) are 
     * determined in the CreateCellWithIntracellularStimulus method
     * @param nodeIndex is the global index within the mesh
     */
    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
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

    ~HeartConfigRelatedCellFactory(void)
    {
    }
};

/**
 * A class which encapsulates the executable functionality.
 * 
 * Takes in a chaste parameters xml file and runs the relevant simulation.
 */
class CardiacSimulation
{
private:
    
    /**
     * Read parameters from the HeartConfig XML file.
     * 
     * @param parameterFileName a string containing the chaste simulation parameters XML file name.
     */
    void ReadParametersFromFile(std::string parameterFileName)
    {
        try
        {
            // Try the hardcoded schema location first
            HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
            HeartConfig::Instance()->SetParametersFile(parameterFileName);
        }
        catch (Exception& e)
        {
            if (e.CheckShortMessageContains("Missing file parsing configuration") == "")
            {
            // Try using the schema location given in the XML
                HeartConfig::Instance()->Reset();
                HeartConfig::Instance()->SetUseFixedSchemaLocation(false);
                HeartConfig::Instance()->SetParametersFile(parameterFileName);
            }
            else
            {
                throw e;
            }
        }    
    }
    
    template<class Problem, unsigned SPACE_DIM>
    void CreateAndRun()
    {
        Problem* p_problem;
        
        if (HeartConfig::Instance()->IsSimulationDefined())
        {
            HeartConfigRelatedCellFactory<SPACE_DIM> cell_factory;
            p_problem = new Problem(&cell_factory);

            p_problem->Initialise();
        }
        else // (HeartConfig::Instance()->IsSimulationResumed())
        {
            p_problem = CardiacSimulationArchiver<Problem>::Load(HeartConfig::Instance()->GetArchivedSimulationDir());                            
        }

        p_problem->Solve();

        /// \todo: #1158 ignoring periodic timestep now 
        if (HeartConfig::Instance()->GetCheckpointSimulation())
        {
            std::stringstream directory;
            directory << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms"; 
            CardiacSimulationArchiver<Problem>::Save(*p_problem, directory.str(), false);
        }
        
        delete p_problem;
    }
    
    /**
     * Run the simulation
     */
    void Run()
    {
        switch (HeartConfig::Instance()->GetDomain())
        {
            case cp::domain_type::Mono :
            {
                switch (HeartConfig::Instance()->GetSpaceDimension())
                {
                    case 3:
                    {
                        CreateAndRun<MonodomainProblem<3>,3>();
                        break;
                    }
    
                    case 2:
                    {
                        CreateAndRun<MonodomainProblem<2>,2>();
                        break;
                    }
    
                    case 1:
                    {
                        CreateAndRun<MonodomainProblem<1>,1>();
                        delete p_mono_problem;
    
                        break;
                    }
                    default :
                        EXCEPTION("Monodomain space dimension not supported: should be 1, 2 or 3");
                }
                break;
            }
    
            case cp::domain_type::Bi :
            {
                switch (HeartConfig::Instance()->GetSpaceDimension())
                {
                    case 3:
                    {
                        CreateAndRun<BidomainProblem<3>,3>();
                        break;
                    }
                    case 2:
                    {
                        CreateAndRun<BidomainProblem<2>,2>();
                        break;
                    }
                    case 1:
                    {
                        CreateAndRun<BidomainProblem<1>,1>();
                        break;
                    }
                    default :
                    {
                        EXCEPTION("Bidomain space dimension not supported: should be 1, 2 or 3");
                    }
                }
                break;
            }
            default :
            {
                //If the domain is not set correctly then the XML parser will have picked it up before now!
                NEVER_REACHED;
            }
        }
    }

public:

    /**
     * Constructor
     * 
     * This also runs the simulation immediately.
     * 
     * @param parameterFileName  The name of the chaste parameters xml file to use to run a simulation (not mandatory since HeartConfig may be set by hand)
     */
    CardiacSimulation(std::string parameterFileName)
    {
        //If we have been passed an XML file then parse the XML file, otherwise throw
        if (parameterFileName == "")
        {
            EXCEPTION("No XML file name given");
        }
        ReadParametersFromFile(parameterFileName);
        Run();
        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }    
  
};

#endif /*CARDIACSIMULATION_HPP_*/
