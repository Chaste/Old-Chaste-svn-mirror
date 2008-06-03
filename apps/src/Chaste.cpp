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

#include "OrthotropicConductivityTensors.hpp"

#include "Hdf5ToMeshalyzerConverter.hpp"

// Path to the parameter file
std::string parameter_file;

// User-modifiable parameters.  Real values will be read from a config file.
double simulation_duration = -1; // ms
double slab_x = -1;          // mm
double slab_y = -1;          // mm
double slab_z = -1;          // mm
double inter_node_space = -1;    // mm

double intra_x_cond;
double intra_y_cond;
double intra_z_cond;

std::string mesh_file_prefix;
anisotropic_type media = anisotropic_type::Orthotropic;

std::string  output_directory = "/";      // Location to put simulation results
std::string  mesh_output_directory = "/"; // Location for generated mesh files
domain_type domain = domain_type::Mono;
ionic_model_type ionic_model = ionic_model_type::LuoRudyIModel1991OdeSystem;

std::vector<InitialStimulus> stimuli_applied;
std::vector<ChasteCuboid> stimuled_areas;

std::vector<double> scale_factor_gks;
std::vector<double> scale_factor_ito;
std::vector<ChasteCuboid> cell_heterogeneity_areas;

std::vector< c_vector<double,3> > specific_conductivities;
std::vector<ChasteCuboid> conductivity_heterogeneity_areas;

bool create_slab;
bool load_mesh;

// Parameters fixed at compile time
const std::string  output_filename_prefix = "Run";
const double ode_time_step = 0.005;     // ms
const double pde_time_step = 0.02;     // ms
const double printing_time_step = 1; // ms

// Scale factor because Chaste code expects lengths in cm, but params use mm.
const double scale_factor = 1/10.0;


class ChasteSlabCellFactory : public AbstractCardiacCellFactory<3>
{
public:
    ChasteSlabCellFactory() : AbstractCardiacCellFactory<3>(ode_time_step)
    {
    }

    
    AbstractCardiacCell* CreateCellWithIntracellularStimulus(AbstractStimulusFunction* intracellularStimulus, unsigned node)
    {
        switch(ionic_model)
        {
            case(ionic_model_type::LuoRudyIModel1991OdeSystem):
                return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, intracellularStimulus, mpZeroStimulus);
                break;   
      
            case(ionic_model_type::BackwardEulerLuoRudyIModel1991):
                return new BackwardEulerLuoRudyIModel1991(mTimeStep, intracellularStimulus, mpZeroStimulus);
                break;                
            
            case(ionic_model_type::BackwardEulerFoxModel2002Modified):
                return new BackwardEulerFoxModel2002Modified(mTimeStep, intracellularStimulus, mpZeroStimulus);
                break;
    
            case(ionic_model_type::FaberRudy2000Version3):
                {
                    FaberRudy2000Version3*  faber_rudy_instance = new FaberRudy2000Version3(mpSolver, mTimeStep, intracellularStimulus, mpZeroStimulus);
                    
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
                return new FaberRudy2000Version3Optimised(mpSolver, mTimeStep, intracellularStimulus, mpZeroStimulus);
                break;
                
            default:
                EXCEPTION("Unknown ionic model!!!");
        }   
        
        return NULL;
    }
    
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {        
        // Memory leak, this pointers should freed somewhere
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
    try
    {
        std::auto_ptr<chaste_parameters_type> p_params(ChasteParameters(parameter_file));
        
        simulation_type simulation_params = p_params->Simulation();
		physiological_type physiological_params = p_params->Physiological();
        
        simulation_duration = simulation_params.SimulationDuration();
        
        create_slab = simulation_params.Mesh().Slab() != NULL;
        load_mesh = simulation_params.Mesh().LoadMesh() != NULL; 
        
        if (create_slab)
        {
            slab_x = simulation_params.Mesh().Slab()->SlabX();     // mm
            slab_y = simulation_params.Mesh().Slab()->SlabY();   // mm
            slab_z = simulation_params.Mesh().Slab()->SlabZ();   // mm
            inter_node_space = simulation_params.Mesh().Slab()->InterNodeSpace(); // mm
        }
        else // (load_mesh)
        {
            mesh_file_prefix = simulation_params.Mesh().LoadMesh()->name();
            media = simulation_params.Mesh().LoadMesh()->media();
        }
        
        intra_x_cond = physiological_params.IntracellularConductivities().longi();
        intra_y_cond = physiological_params.IntracellularConductivities().trans();
        intra_z_cond = physiological_params.IntracellularConductivities().normal();                
        output_directory = simulation_params.OutputDirectory();
        mesh_output_directory = simulation_params.MeshOutputDirectory();
        domain = simulation_params.Domain();
        ionic_model = simulation_params.IonicModel();
        
        // Read and store Stimuli
        simulation_type::Stimulus::container& stimuli = simulation_params.Stimulus();
        for (simulation_type::Stimulus::iterator i = stimuli.begin();
             i != stimuli.end();
             ++i)
        {                     
            stimulus_type stimulus(*i);           
            point_type point_a = stimulus.Location().CornerA();
            point_type point_b = stimulus.Location().CornerB();
            
            ChastePoint<3> chaste_point_a (scale_factor* point_a.x(), 
                                           scale_factor* point_a.y(),
                                           scale_factor* point_a.z());

            ChastePoint<3> chaste_point_b (scale_factor* point_b.x(),
                                           scale_factor* point_b.y(),
                                           scale_factor* point_b.z());
                        
            stimuli_applied.push_back( InitialStimulus(stimulus.Strength(), stimulus.Duration(), stimulus.Delay() ) );
            stimuled_areas.push_back( ChasteCuboid( chaste_point_a, chaste_point_b ) );
        }

        // Read and store Cell Heterogeneities
        simulation_type::CellHeterogeneity::container& hts = simulation_params.CellHeterogeneity();
        for (simulation_type::CellHeterogeneity::iterator i = hts.begin();
             i != hts.end();
             ++i)
        {                     
            cell_heterogeneity_type ht(*i);           
            point_type point_a = ht.Location().CornerA();
            point_type point_b = ht.Location().CornerB();
            
            // method get() should be called for Y and Z since they have been defined optional in the schema
            // {Y,Z}.set() can be called to know if they have been defined
            ChastePoint<3> chaste_point_a (scale_factor* point_a.x(), 
                                           scale_factor* point_a.y(),
                                           scale_factor* point_a.z());

            ChastePoint<3> chaste_point_b (scale_factor* point_b.x(),
                                           scale_factor* point_b.y(),
                                           scale_factor* point_b.z());
                        
            scale_factor_gks.push_back (ht.ScaleFactorGks());
            scale_factor_ito.push_back (ht.ScaleFactorIto());                                    
            cell_heterogeneity_areas.push_back( ChasteCuboid( chaste_point_a, chaste_point_b ) );
        }
        
        // Read and store Conductivity Heterogeneities

    }
    catch (const xml_schema::exception& e)
    {
         std::cerr << e << std::endl;
         EXCEPTION("XML parsing error");
    }
}

template<unsigned PROBLEM_DIM>
void SetupProblem(AbstractCardiacProblem<3, PROBLEM_DIM>& rProblem)
{
    rProblem.SetEndTime(simulation_duration);   // ms
    rProblem.SetPdeTimeStepAndPrintEveryNthTimeStep(pde_time_step, (int) (printing_time_step/pde_time_step));
    rProblem.SetOutputDirectory(output_directory+"/results");
    rProblem.SetOutputFilenamePrefix("Chaste");
    rProblem.SetCallChaste2Meshalyzer(false);
    
    if (load_mesh)
    {
        std::string fibre_file = mesh_file_prefix + ".fibres";
        rProblem.SetFibreOrientation(fibre_file);
    }
    
    rProblem.SetIntracellularConductivities(Create_c_vector(intra_x_cond, intra_y_cond, intra_z_cond));
}


void CreateSlab(ConformingTetrahedralMesh<3,3>* pMesh)
{
    // construct mesh. Note that mesh is measured in cm
    unsigned slab_nodes_x = (unsigned)round(slab_x/inter_node_space);
    unsigned slab_nodes_y = (unsigned)round(slab_y/inter_node_space);
    unsigned slab_nodes_z = (unsigned)round(slab_z/inter_node_space);        
   
    pMesh->ConstructCuboid(slab_nodes_x,
                         slab_nodes_y,
                         slab_nodes_z,
                         true);
    // place at origin
    pMesh->Translate(-(double)slab_nodes_x/2.0,
                   -(double)slab_nodes_y/2.0,
                   -(double)slab_nodes_z/2.0);
    // scale
    double mesh_scale_factor = inter_node_space*scale_factor;
    pMesh->Scale(mesh_scale_factor, mesh_scale_factor, mesh_scale_factor);

    OutputFileHandler handler(output_directory,false);
    std::string output_dir_full_path = handler.GetOutputDirectoryFullPath();


    // write out the mesh that was used if we are the master process
    if (PetscTools::AmMaster())
    {
        // Meshalyzer output format
        //MeshalyzerMeshWriter<3,3> mesh_writer(output_directory+"/mesh", "Slab", false);
        //mesh_writer.WriteFilesUsingMesh(mesh);
        // Triangles output format
        TrianglesMeshWriter<3,3> triangles_writer(output_directory+"/mesh", "Slab", false);
        triangles_writer.WriteFilesUsingMesh(*pMesh);
        
        // copy input parameters file to results directory
        //system(("cp " + parameter_file + " " + output_dir_full_path).c_str());
    }    
}


int main(int argc, char *argv[]) 
{
    std::cout << "Copyright (C) University of Oxford, 2008 \n\n\
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
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.\n\n ";

    try
    {
        PETSCEXCEPT(PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL) );
        
        // solver and preconditioner options
        //PetscOptionsSetValue("-ksp_type", "cg");
        //PetscOptionsSetValue("-pc_type", "bjacobi");
        //PetscOptionsSetValue("-options_table", "");
        
        if (argc<2)
        {
            std::cout  << "Usage: Chaste parameters_file\n";
            return -1;
        }
        
        parameter_file = std::string(argv[1]);
        ReadParametersFromFile();
                
        ChasteSlabCellFactory cell_factory;       
        ConformingTetrahedralMesh<3,3> mesh;
        
        bool orthotropic = (media == anisotropic_type::Orthotropic);
        
        switch(domain)
        {
            case domain_type::Mono :
            {            
                MonodomainProblem<3> mono_problem( &cell_factory, orthotropic );

                SetupProblem(mono_problem);
                                
                if (create_slab)
                {
                    CreateSlab(&mesh);                               
                    mono_problem.SetMesh(&mesh);     
                }
                else // (load_mesh) 
                {
                    mono_problem.SetMeshFilename(mesh_file_prefix);
                }
                
                mono_problem.Initialise();
                mono_problem.Solve();
        
                break;
            }
            
            case domain_type::Bi :
            {
                BidomainProblem<3> bi_problem( &cell_factory, orthotropic );                

                SetupProblem(bi_problem);
                bi_problem.SetExtracellularConductivities(Create_c_vector(6.2, 2.4, 2.4));
                                
                if (create_slab)
                {   
                    CreateSlab(&mesh);                            
                    bi_problem.SetMesh(&mesh);     
                }
                else // (load_mesh) 
                {
                    bi_problem.SetMeshFilename(mesh_file_prefix);
                }
                
                bi_problem.Initialise();
                bi_problem.Solve();

                break;
            }
                
            default :
                EXCEPTION("Unknown domain type!!!");
        }        
        
    }
    catch(Exception& e)
    {
        std::cerr << e.GetMessage() << "\n";
        return 1;
    }
    
    Hdf5ToMeshalyzerConverter converter(output_directory+"/results", "Chaste");
    
    EventHandler::Headings();
    EventHandler::Report();

    PetscFinalize();

    return 0;    
}    


