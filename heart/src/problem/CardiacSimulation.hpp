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
 *  todo: dox, coverage, maybe own tests
 */
template<unsigned SPACE_DIM>
class HeartConfigRelatedCellFactory : public AbstractCardiacCellFactory<SPACE_DIM>
{
private:
    cp::ionic_models_available_type mDefaultIonicModel;
    std::vector<ChasteCuboid> mIonicModelRegions;
    std::vector<cp::ionic_models_available_type> mIonicModelsDefined;
    
    std::vector<boost::shared_ptr<SimpleStimulus> > mStimuliApplied;
    std::vector<ChasteCuboid> mStimulatedAreas;
    
    std::vector<double> mScaleFactorGks;
    std::vector<double> mScaleFactorIto;
    std::vector<double> mScaleFactorGkr;
    std::vector<ChasteCuboid> mCellHeterogeneityAreas;

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
            // No stimuli provided
            std::cout << "Warning: No stimuli provided. Simulation will be run anyway." << std::endl;
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


    AbstractCardiacCell* CreateCellWithIntracellularStimulus(boost::shared_ptr<AbstractStimulusFunction> intracellularStimulus, unsigned node)
    {
        cp::ionic_models_available_type ionic_model = mDefaultIonicModel;

        for (unsigned ionic_model_region_index = 0;
             ionic_model_region_index < mIonicModelRegions.size();
             ++ionic_model_region_index)
        {
            if ( mIonicModelRegions[ionic_model_region_index].DoesContain(this->GetMesh()->GetNode(node)->GetPoint()) )
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
                    if ( mCellHeterogeneityAreas[ht_index].DoesContain(this->GetMesh()->GetNode(node)->GetPoint()) )
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
                    if ( mCellHeterogeneityAreas[ht_index].DoesContain(this->GetMesh()->GetNode(node)->GetPoint()) )
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
                EXCEPTION("Unknown ionic model!!!");
            }
        }

        return NULL;
    }


    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        boost::shared_ptr<MultiStimulus> node_specific_stimulus(new MultiStimulus());

        // Check which of the defined stimuli contain the current node
        for (unsigned stimulus_index = 0;
             stimulus_index < mStimuliApplied.size();
             ++stimulus_index)
        {
            if ( mStimulatedAreas[stimulus_index].DoesContain(this->GetMesh()->GetNode(node)->GetPoint()) )
            {
                node_specific_stimulus->AddStimulus(mStimuliApplied[stimulus_index]);
            }
        }

        return CreateCellWithIntracellularStimulus(node_specific_stimulus, node);
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
    /*
     * 
     * User-modifiable parameters.  Real values will be read from the config file.
     * Defaults are required in the declaration due to the lack of default constructor.
     * Values will be ignored, though.
     */
     
    /** Domain type, 'Mono' or 'Bi' in the xml */
    cp::domain_type mDomain;
    /** Space dimension (1,2 or 3) */
    unsigned mSpaceDimension;
    
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
                throw;
            }
        }
    
    
        if (HeartConfig::Instance()->IsSimulationDefined())
        {
            mDomain = HeartConfig::Instance()->GetDomain();
            mSpaceDimension = HeartConfig::Instance()->GetSpaceDimension();

        }
        else
        {
            mDomain = HeartConfig::Instance()->GetDomain();
            mSpaceDimension = HeartConfig::Instance()->GetSpaceDimension();
        }
    }
    
    /**
     * Run the simulation
     */
    void Run()
    {
//// TODO: remember to #include PetscArguments in executable, and #include PetscSetupAndFinalize in tests (as usual)
        
        switch(mDomain)
        {
            case cp::domain_type::Mono :
            {
                switch (mSpaceDimension)
                {
                    case 3:
                    {
                        MonodomainProblem<3>* p_mono_problem;
                        
                        if (HeartConfig::Instance()->IsSimulationDefined())
                        {
                            HeartConfigRelatedCellFactory<3> cell_factory;
                            p_mono_problem = new MonodomainProblem<3>(&cell_factory);
    
                            p_mono_problem->Initialise();
    
                            p_mono_problem->ConvertOutputToMeshalyzerFormat(true);
                        }
                        else // (HeartConfig::Instance()->IsSimulationResumed())
                        {
                            p_mono_problem = CardiacSimulationArchiver<MonodomainProblem<3> >::Load(HeartConfig::Instance()->GetArchivedSimulationDir());                            
                        }

                        p_mono_problem->Solve();

                        if (HeartConfig::Instance()->GetSaveSimulation())
                        {
                            std::stringstream directory;
                            directory << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms"; 
                            CardiacSimulationArchiver<MonodomainProblem<3> >::Save(*p_mono_problem, directory.str(), false);
                        }

                        delete p_mono_problem;

                        break;
                    }

                    case 2:
                    {
                        MonodomainProblem<2>* p_mono_problem;
                        
                        if (HeartConfig::Instance()->IsSimulationDefined())
                        {
                            HeartConfigRelatedCellFactory<2> cell_factory;
                            p_mono_problem = new MonodomainProblem<2>(&cell_factory);
    
                            p_mono_problem->Initialise();
    
                            p_mono_problem->ConvertOutputToMeshalyzerFormat(true);
                        }
                        else // (HeartConfig::Instance()->IsSimulationResumed())
                        {
                            p_mono_problem = CardiacSimulationArchiver<MonodomainProblem<2> >::Load(HeartConfig::Instance()->GetArchivedSimulationDir());                            
                        }

                        p_mono_problem->Solve();

                        if (HeartConfig::Instance()->GetSaveSimulation())
                        {
                            std::stringstream directory;
                            directory << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms"; 
                            CardiacSimulationArchiver<MonodomainProblem<2> >::Save(*p_mono_problem, directory.str(), false);
                        }

                        delete p_mono_problem;

                        break;
                    }

                    case 1:
                    {
                        MonodomainProblem<1>* p_mono_problem;
                        
                        if (HeartConfig::Instance()->IsSimulationDefined())
                        {
                            HeartConfigRelatedCellFactory<1> cell_factory;
                            p_mono_problem = new MonodomainProblem<1>(&cell_factory);
    
                            p_mono_problem->Initialise();
    
                            p_mono_problem->ConvertOutputToMeshalyzerFormat(true);
                        }
                        else // (HeartConfig::Instance()->IsSimulationResumed())
                        {
                            p_mono_problem = CardiacSimulationArchiver<MonodomainProblem<1> >::Load(HeartConfig::Instance()->GetArchivedSimulationDir());                            
                        }

                        p_mono_problem->Solve();

                        if (HeartConfig::Instance()->GetSaveSimulation())
                        {
                            std::stringstream directory;
                            directory << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms"; 
                            CardiacSimulationArchiver<MonodomainProblem<1> >::Save(*p_mono_problem, directory.str(), false);
                        }

                        delete p_mono_problem;

                        break;
                    }
                    default :
                        EXCEPTION("Space dimension not supported: should be 1, 2 or 3");
                }
                break;
            }

            case cp::domain_type::Bi :
            {
                switch (mSpaceDimension)
                {
                    case 3:
                    {
                        BidomainProblem<3>* p_bi_problem;
                        
                        if (HeartConfig::Instance()->IsSimulationDefined())
                        {
                            HeartConfigRelatedCellFactory<3> cell_factory;
                            p_bi_problem = new BidomainProblem<3>(&cell_factory);
    
                            p_bi_problem->Initialise();
    
                            p_bi_problem->ConvertOutputToMeshalyzerFormat(true);
                        }
                        else // (HeartConfig::Instance()->IsSimulationResumed())
                        {
                            p_bi_problem = CardiacSimulationArchiver<BidomainProblem<3> >::Load(HeartConfig::Instance()->GetArchivedSimulationDir());                            
                        }

                        p_bi_problem->Solve();

                        if (HeartConfig::Instance()->GetSaveSimulation())
                        {
                            std::stringstream directory;
                            directory << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms"; 
                            CardiacSimulationArchiver<BidomainProblem<3> >::Save(*p_bi_problem, directory.str(), false);
                        }

                        delete p_bi_problem;

                        break;
                    }
                    case 2:
                    {
                        BidomainProblem<2>* p_bi_problem;
                        
                        if (HeartConfig::Instance()->IsSimulationDefined())
                        {
                            HeartConfigRelatedCellFactory<2> cell_factory;
                            p_bi_problem = new BidomainProblem<2>(&cell_factory);
    
                            p_bi_problem->Initialise();
    
                            p_bi_problem->ConvertOutputToMeshalyzerFormat(true);
                        }
                        else // (HeartConfig::Instance()->IsSimulationResumed())
                        {
                            p_bi_problem = CardiacSimulationArchiver<BidomainProblem<2> >::Load(HeartConfig::Instance()->GetArchivedSimulationDir());                            
                        }

                        p_bi_problem->Solve();

                        if (HeartConfig::Instance()->GetSaveSimulation())
                        {
                            std::stringstream directory;
                            directory << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms"; 
                            CardiacSimulationArchiver<BidomainProblem<2> >::Save(*p_bi_problem, directory.str(), false);
                        }

                        delete p_bi_problem;

                        break;
                    }
                    case 1:
                    {
                        BidomainProblem<1>* p_bi_problem;
                        
                        if (HeartConfig::Instance()->IsSimulationDefined())
                        {
                            HeartConfigRelatedCellFactory<1> cell_factory;
                            p_bi_problem = new BidomainProblem<1>(&cell_factory);
    
                            p_bi_problem->Initialise();
    
                            p_bi_problem->ConvertOutputToMeshalyzerFormat(true);
                        }
                        else // (HeartConfig::Instance()->IsSimulationResumed())
                        {
                            p_bi_problem = CardiacSimulationArchiver<BidomainProblem<1> >::Load(HeartConfig::Instance()->GetArchivedSimulationDir());                            
                        }

                        p_bi_problem->Solve();

                        if (HeartConfig::Instance()->GetSaveSimulation())
                        {
                            std::stringstream directory;
                            directory << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms"; 
                            CardiacSimulationArchiver<BidomainProblem<1> >::Save(*p_bi_problem, directory.str(), false);
                        }

                        delete p_bi_problem;
                        break;
                    }
                    default :
                    {
                        EXCEPTION("Space dimension not supported: should be 1, 2 or 3");
                    }
                }
                break;
            }
            default :
            {
                EXCEPTION("Unknown domain type: should be 'Mono' or 'Bi'");
            }
        }

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
    
public:

    /**
     * Constructor
     * 
     * This also runs the simulation immediately.
     * 
     * @param parameterFileName  The name of the chaste parameters xml file to run a simulation with
     */
    CardiacSimulation(std::string parameterFileName)
        : mDomain(cp::domain_type::Mono),
          mSpaceDimension(3u)
    {
        ReadParametersFromFile(parameterFileName);
        Run();
    }    
  
};

#endif /*CARDIACSIMULATION_HPP_*/
