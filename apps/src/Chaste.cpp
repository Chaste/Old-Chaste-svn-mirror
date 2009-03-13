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

#include "TetrahedralMesh.hpp"
#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include "PetscTools.hpp"
#include <petscvec.h>
#include <vector>
#include "AbstractCardiacCellFactory.hpp"
#include "TimeStepper.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "TrianglesMeshWriter.hpp"
#include "MultiStimulus.hpp"
#include <ctime>

#include "ChasteParameters.hpp"
#include <memory>

#include "AbstractStimulusFunction.hpp"
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"

#include "BackwardEulerFoxModel2002Modified.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BackwardEulerLuoRudyIModel1991.hpp"

#include "FoxModel2002Modified.hpp"
#include "FaberRudy2000Version3.cpp"
#include "FaberRudy2000Version3Optimised.hpp"

#include "DiFrancescoNoble1985OdeSystem.hpp"
#include "Mahajan2008OdeSystem.hpp"
#include "TenTusscher2006OdeSystem.hpp"
#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"

#include "OrthotropicConductivityTensors.hpp"

#include "Hdf5ToMeshalyzerConverter.hpp"

// Path to the parameter file
std::string parameter_file;

// User-modifiable parameters.  Real values will be read from a config file.
std::string  output_directory = "/";      // Location to put simulation results
domain_type domain = domain_type::Mono;
ionic_model_type ionic_model = ionic_model_type::LuoRudyI;

std::vector<SimpleStimulus> stimuli_applied;
std::vector<ChasteCuboid> stimuled_areas;

std::vector<double> scale_factor_gks;
std::vector<double> scale_factor_ito;
std::vector<ChasteCuboid> cell_heterogeneity_areas;

std::vector< c_vector<double,3> > specific_conductivities;
std::vector<ChasteCuboid> conductivity_heterogeneity_areas;

class ChasteSlabCellFactory : public AbstractCardiacCellFactory<3>
{
public:
    ChasteSlabCellFactory() : AbstractCardiacCellFactory<3>()
    {
    }


    AbstractCardiacCell* CreateCellWithIntracellularStimulus(AbstractStimulusFunction* intracellularStimulus, unsigned node)
    {
        switch(ionic_model)
        {
            case(ionic_model_type::LuoRudyI):
                return new LuoRudyIModel1991OdeSystem(mpSolver, intracellularStimulus);
                break;

            case(ionic_model_type::BackwardEulerLuoRudyIModel1991):
                return new BackwardEulerLuoRudyIModel1991(intracellularStimulus);
                break;

            case(ionic_model_type::BackwardEulerFoxModel2002Modified):
                return new BackwardEulerFoxModel2002Modified(intracellularStimulus);
                break;
            
            case(ionic_model_type::DifrancescoNoble):
                return new DiFrancescoNoble1985OdeSystem(mpSolver, intracellularStimulus);
                break;
            
            case(ionic_model_type::MahajanShiferaw):
                return new Mahajan2008OdeSystem(mpSolver, intracellularStimulus);
                break;
                
            case(ionic_model_type::tenTusscher2006):
                return new TenTusscher2006OdeSystem(mpSolver, intracellularStimulus);
                break;
            
            case(ionic_model_type::HodgkinHuxley):
                return new HodgkinHuxleySquidAxon1952OriginalOdeSystem(mpSolver, intracellularStimulus);
                break;
                
            case(ionic_model_type::FaberRudy2000):
                {
                    FaberRudy2000Version3*  faber_rudy_instance = new FaberRudy2000Version3(mpSolver, intracellularStimulus);

                    for (unsigned ht_index = 0;
                         ht_index < cell_heterogeneity_areas.size();
                         ++ht_index)
                    {
                        if ( cell_heterogeneity_areas[ht_index].DoesContain(mpMesh->GetNode(node)->GetPoint()) )
                        {
                            faber_rudy_instance->SetScaleFactorGks(scale_factor_gks[ht_index]);
                            faber_rudy_instance->SetScaleFactorIto(scale_factor_ito[ht_index]);
                        }
                    }

                    return faber_rudy_instance;
                    break;
                }

            case(ionic_model_type::FaberRudy2000Version3Optimised):
                return new FaberRudy2000Version3Optimised(mpSolver, intracellularStimulus);
                break;

            default:
                EXCEPTION("Unknown ionic model!!!");
        }

        return NULL;
    }


    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        // Memory leak, these pointers should freed somewhere
        MultiStimulus* node_specific_stimulus = new MultiStimulus();

        // Check which of the defined stimuli contain the current node
        for (unsigned stimulus_index = 0;
             stimulus_index < stimuli_applied.size();
             ++stimulus_index)
        {
            if ( stimuled_areas[stimulus_index].DoesContain(mpMesh->GetNode(node)->GetPoint()) )
            {
                node_specific_stimulus->AddStimulus(&stimuli_applied[stimulus_index]);
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
    HeartConfig::Instance()->SetParametersFile(parameter_file);

    output_directory = HeartConfig::Instance()->GetOutputDirectory();

    domain = HeartConfig::Instance()->GetDomain();
    ionic_model = HeartConfig::Instance()->GetIonicModel();

    // Read and store Stimuli
    HeartConfig::Instance()->GetStimuli(stimuli_applied, stimuled_areas);

    // Read and store Cell Heterogeneities
    try
    {
        HeartConfig::Instance()->GetCellHeterogeneities(cell_heterogeneity_areas,
                                                        scale_factor_gks,
                                                        scale_factor_ito);
    }
    catch(Exception& e)
    {
        // Ignore the exception
    }               
    
    /// \todo Read and store Conductivity Heterogeneities

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
    std::cout<<UNAME<<"\n";
    std::cout<<"from revision number "<<SVN_REV<<" with build type "<<BUILD_TYPE<<".\n\n";
    
    PETSCEXCEPT(PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL) );

    if (argc<2)
    {
        std::cout  << "Usage: Chaste parameters_file\n";
        return -1;
    }

    parameter_file = std::string(argv[1]);
    ReadParametersFromFile();

    ChasteSlabCellFactory cell_factory;

    switch(domain)
    {
        case domain_type::Mono :
        {
            MonodomainProblem<3> mono_problem( &cell_factory);

            mono_problem.ConvertOutputToMeshalyzerFormat(true);

            mono_problem.Initialise();
            mono_problem.Solve();

            break;
        }

        case domain_type::Bi :
        {
            BidomainProblem<3> bi_problem( &cell_factory);

            bi_problem.ConvertOutputToMeshalyzerFormat(true);

            bi_problem.Initialise();
            bi_problem.Solve();

            break;
        }

        default :
            EXCEPTION("Unknown domain type!!!");
    }
            
    HeartEventHandler::Headings();
    HeartEventHandler::Report();

    PetscFinalize();

    return 0;
}


