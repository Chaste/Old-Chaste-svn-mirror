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
#include "ChasteParameters.hpp"

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

// Path to the parameter file
std::string parameter_file;

/*
 * User-modifiable parameters.  Real values will be read from the config file.
 * Defaults are required in the declaration due to the lack of default constructor.
 * Values will be ignored, though.
 */
std::string  output_directory = "/";
domain_type domain = domain_type::Mono;

ionic_models_available_type default_ionic_model = ionic_models_available_type::LuoRudyI;
std::vector<ChasteCuboid> ionic_model_regions;
std::vector<ionic_models_available_type> ionic_models_defined;

std::vector<boost::shared_ptr<SimpleStimulus> > stimuli_applied;
std::vector<ChasteCuboid> stimulated_areas;

std::vector<double> scale_factor_gks;
std::vector<double> scale_factor_ito;
std::vector<double> scale_factor_gkr;
std::vector<ChasteCuboid> cell_heterogeneity_areas;

template<unsigned SPACE_DIM>
class ChasteSlabCellFactory : public AbstractCardiacCellFactory<SPACE_DIM>
{
public:
    ChasteSlabCellFactory() : AbstractCardiacCellFactory<SPACE_DIM>()
    {
    }


    AbstractCardiacCell* CreateCellWithIntracellularStimulus(boost::shared_ptr<AbstractStimulusFunction> intracellularStimulus, unsigned node)
    {
        ionic_models_available_type ionic_model = default_ionic_model;

        for (unsigned ionic_model_region_index = 0;
             ionic_model_region_index < ionic_model_regions.size();
             ++ionic_model_region_index)
        {
            if ( ionic_model_regions[ionic_model_region_index].DoesContain(this->GetMesh()->GetNode(node)->GetPoint()) )
            {
                ionic_model = ionic_models_defined[ionic_model_region_index];
                break;
            }
        }

        switch(ionic_model)
        {
            case(ionic_models_available_type::LuoRudyI):
                return new LuoRudyIModel1991OdeSystem(this->mpSolver, intracellularStimulus);
                break;

            case(ionic_models_available_type::LuoRudyIBackwardEuler):
                return new BackwardEulerLuoRudyIModel1991(intracellularStimulus);
                break;

            case(ionic_models_available_type::Fox2002BackwardEuler):
                return new BackwardEulerFoxModel2002Modified(intracellularStimulus);
                break;

            case(ionic_models_available_type::DifrancescoNoble):
                return new DiFrancescoNoble1985OdeSystem(this->mpSolver, intracellularStimulus);
                break;

            case(ionic_models_available_type::MahajanShiferaw):
                return new Mahajan2008OdeSystem(this->mpSolver, intracellularStimulus);
                break;

            case(ionic_models_available_type::tenTusscher2006):
                {
                    TenTusscher2006OdeSystem*  tt06_instance = new TenTusscher2006OdeSystem(this->mpSolver, intracellularStimulus);

                    for (unsigned ht_index = 0;
                         ht_index < cell_heterogeneity_areas.size();
                         ++ht_index)
                    {
                        if ( cell_heterogeneity_areas[ht_index].DoesContain(this->GetMesh()->GetNode(node)->GetPoint()) )
                        {
                            tt06_instance->SetScaleFactorGks(scale_factor_gks[ht_index]);
                            tt06_instance->SetScaleFactorIto(scale_factor_ito[ht_index]);
                            tt06_instance->SetScaleFactorGkr(scale_factor_gkr[ht_index]);
                        }
                    }

                    return tt06_instance;
                    break;
                }


            case(ionic_models_available_type::HodgkinHuxley):
                return new HodgkinHuxleySquidAxon1952OriginalOdeSystem(this->mpSolver, intracellularStimulus);
                break;

            case(ionic_models_available_type::FaberRudy2000):
                {
                    FaberRudy2000Version3*  faber_rudy_instance = new FaberRudy2000Version3(this->mpSolver, intracellularStimulus);

                    for (unsigned ht_index = 0;
                         ht_index < cell_heterogeneity_areas.size();
                         ++ht_index)
                    {
                        if ( cell_heterogeneity_areas[ht_index].DoesContain(this->GetMesh()->GetNode(node)->GetPoint()) )
                        {
                            faber_rudy_instance->SetScaleFactorGks(scale_factor_gks[ht_index]);
                            faber_rudy_instance->SetScaleFactorIto(scale_factor_ito[ht_index]);
                            faber_rudy_instance->SetScaleFactorGkr(scale_factor_gkr[ht_index]);
                        }
                    }

                    return faber_rudy_instance;
                    break;
                }

            case(ionic_models_available_type::FaberRudy2000Optimised):
                return new FaberRudy2000Version3Optimised(this->mpSolver, intracellularStimulus);
                break;

            default:
                EXCEPTION("Unknown ionic model!!!");
        }

        return NULL;
    }


    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        boost::shared_ptr<MultiStimulus> node_specific_stimulus(new MultiStimulus());

        // Check which of the defined stimuli contain the current node
        for (unsigned stimulus_index = 0;
             stimulus_index < stimuli_applied.size();
             ++stimulus_index)
        {
            if ( stimulated_areas[stimulus_index].DoesContain(this->GetMesh()->GetNode(node)->GetPoint()) )
            {
                node_specific_stimulus->AddStimulus(stimuli_applied[stimulus_index]);
            }
        }

        return CreateCellWithIntracellularStimulus(node_specific_stimulus, node);
    }

    ~ChasteSlabCellFactory(void)
    {
    }
};

void ReadParametersFromFile()
{
    try
    {
        // Try using the schema location given in the XML first
        HeartConfig::Instance()->SetUseFixedSchemaLocation(false);
        HeartConfig::Instance()->SetParametersFile(parameter_file);
    }
    catch (Exception& e)
    {
        // Use the fixed location, but warn the user
        std::cerr << "Failed to load parameters file using schema specified in file"
                  << " (error was: " << e.GetMessage() << ");"
                  << " using built-in default schema location." << std::endl << std::flush;
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        HeartConfig::Instance()->SetParametersFile(parameter_file);
    }

    output_directory = HeartConfig::Instance()->GetOutputDirectory();

    domain = HeartConfig::Instance()->GetDomain();

    // Read and store default ionic model and possible region definitions
    default_ionic_model = HeartConfig::Instance()->GetDefaultIonicModel();
    HeartConfig::Instance()->GetIonicModelRegions(ionic_model_regions,
                                                  ionic_models_defined);

    // Read and store Stimuli
    try
    {
        HeartConfig::Instance()->GetStimuli(stimuli_applied, stimulated_areas);
    }
    catch(Exception& e)
    {
        // No stimuli provided
        std::cout << "Warning: No stimuli provided. Simulation will be run anyway." << std::endl;
    }

    // Read and store Cell Heterogeneities
    try
    {
        HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas,
                                                        scale_factor_gks,
                                                        scale_factor_ito,
                                                        scale_factor_gkr);
    }
    catch(Exception& e)
    {
        // No cell heterogeneities provided
    }
}


int main(int argc, char *argv[])
{
    std::cout << "Copyright (C) University of Oxford, 2005-2009 \n\n\
\
Chaste is free software: you can redistribute it and/or modify \n\
it under the terms of the Lesser GNU General Public License as published by \n\
the Free Software Foundation, either version 2.1 of the License, or \n\
(at your option) any later version. \n\n\
\
Chaste is distributed in the hope that it will be useful, \n\
but WITHOUT ANY WARRANTY; without even the implied warranty of \n\
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n\
Lesser GNU General Public License for more details. \n\n\
\
You should have received a copy of the Lesser GNU General Public License \n\
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.\n\n";

    //Compilation information
    std::cout<<"This version of Chaste was compiled on:\n";
    std::cout<<UNAME<<" (uname)\n";
    std::cout<<"from revision number "<<GetChasteVersion()<<" with build type "<<BUILD_TYPE<<".\n\n";

    try
    {
        PETSCEXCEPT(PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL) );

        if (argc<2)
        {
            std::cout  << "Usage: Chaste parameters_file\n";
            return -1;
        }

        parameter_file = std::string(argv[1]);
        ReadParametersFromFile();

        switch(domain)
        {
            case domain_type::Mono :
            {
                switch (HeartConfig::Instance()->GetSpaceDimension())
                {
                    case 3:
                    {
                        ChasteSlabCellFactory<3> cell_factory;
                        MonodomainProblem<3> mono_problem( &cell_factory);

                        mono_problem.ConvertOutputToMeshalyzerFormat(true);
                        mono_problem.Initialise();
                        mono_problem.Solve();

                        break;
                    }

                    case 2:
                    {
                        ChasteSlabCellFactory<2> cell_factory;
                        MonodomainProblem<2> mono_problem( &cell_factory);

                        mono_problem.ConvertOutputToMeshalyzerFormat(true);

                        mono_problem.Initialise();
                        mono_problem.Solve();

                        break;
                    }

                    case 1:
                    {
                        ChasteSlabCellFactory<1> cell_factory;
                        MonodomainProblem<1> mono_problem( &cell_factory);

                        mono_problem.ConvertOutputToMeshalyzerFormat(true);

                        mono_problem.Initialise();
                        mono_problem.Solve();

                        break;
                    }
                    default :
                        EXCEPTION("Space dimension not supported!!!");
                }
                break;
            }

            case domain_type::Bi :
            {
                switch (HeartConfig::Instance()->GetSpaceDimension())
                {
                    case 3:
                    {
                        ChasteSlabCellFactory<3> cell_factory;
                        BidomainProblem<3> bi_problem( &cell_factory);

                        bi_problem.ConvertOutputToMeshalyzerFormat(true);
                        bi_problem.Initialise();
                        bi_problem.Solve();

                        break;
                    }
                    case 2:
                    {
                        ChasteSlabCellFactory<2> cell_factory;
                        BidomainProblem<2> bi_problem( &cell_factory);

                        bi_problem.ConvertOutputToMeshalyzerFormat(true);
                        bi_problem.Initialise();
                        bi_problem.Solve();

                        break;
                    }
                    case 1:
                    {
                        ChasteSlabCellFactory<1> cell_factory;
                        BidomainProblem<1> bi_problem( &cell_factory);

                        bi_problem.ConvertOutputToMeshalyzerFormat(true);
                        bi_problem.Initialise();
                        bi_problem.Solve();

                        break;
                    }
                    default :
                        EXCEPTION("Space dimension not supported!!!");
                }
                break;
            }
            default :
                EXCEPTION("Unknown domain type!!!");
        }
        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
    catch(Exception& e)
    {
        std::cout << e.GetMessage() << std::endl;
    }

    PetscFinalize();

    return 0;
}


