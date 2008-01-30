#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include "PetscTools.hpp"
#include <petscvec.h>
#include <vector>
#include "AbstractCardiacCellFactory.hpp"
#include "TimeStepper.hpp"
#include "MeshalyzerMeshWriter.cpp"
#include "TrianglesMeshWriter.cpp"
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

// Path to the parameter file
std::string parameter_file;

// User-modifiable parameters.  Real values will be read from a config file.
double simulation_duration = -1; // ms
double slab_width = -1;          // mm
double slab_height= -1;          // mm
double inter_node_space = -1;    // mm
double face_stimulus_width = -1; // mm
double quadrant_stimulus_delay = -1; // ms
std::string  output_directory = "/";      // Location to put simulation results
std::string  mesh_output_directory = "/"; // Location for generated mesh files
domain_type domain = domain_type::Mono;
ionic_model_type ionic_model = ionic_model_type::LuoRudyIModel1991OdeSystem;

std::vector<InitialStimulus> stimuli_applied;
std::vector<ChasteCuboid> stimuled_areas;

// Parameters fixed at compile time
const std::string  output_filename_prefix = "Run";
const double ode_time_step = 0.005;     // ms
const double pde_time_step = 0.02;     // ms
const double printing_time_step = 1; // ms

const double intracellular_cond = 1.75;
const double extracellular_cond = 7.0;

// Scale factor because Chaste code expects lengths in cm, but params use mm.
const double scale_factor = 1/10.0;


class ChasteSlabCellFactory : public AbstractCardiacCellFactory<3>
{
public:
    ChasteSlabCellFactory() : AbstractCardiacCellFactory<3>(ode_time_step)
    {
    }

    
    AbstractCardiacCell* CreateCellWithIntracellularStimulus(AbstractStimulusFunction* intracellularStimulus)
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
                return new FaberRudy2000Version3(mpSolver, mTimeStep, intracellularStimulus, mpZeroStimulus);
                break;

            case(ionic_model_type::FaberRudy2000Version3Optimised):
                return new FaberRudy2000Version3Optimised(mpSolver, mTimeStep, intracellularStimulus, mpZeroStimulus);
                break;
                
            default:
                EXCEPTION("Unknown ionic model!!!");
        }   
        
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
        
        return CreateCellWithIntracellularStimulus(node_specific_stimulus);                
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
        simulation_duration = p_params->SimulationDuration();
        slab_width = p_params->SlabWidth();     // mm
        slab_height = p_params->SlabHeight();   // mm
        inter_node_space = p_params->InterNodeSpace(); // mm
        output_directory = p_params->OutputDirectory();
        mesh_output_directory = p_params->MeshOutputDirectory();
        domain = p_params->Domain();
        ionic_model = p_params->IonicModel();
        
        chaste_parameters_type::Stimulus::container& stimuli = p_params->Stimulus();
        
        for (chaste_parameters_type::Stimulus::iterator i = stimuli.begin();
             i != stimuli.end();
             ++i)
        {                     
            stimulus_type stimulus(*i);           
            point_type point_a = stimulus.Location().CornerA();
            point_type point_b = stimulus.Location().CornerB();
            
            // method get() should be called for Y and Z since they have been defined optional in the schema
            // {Y,Z}.set() can be called to know if they have been defined
            ChastePoint<3> chaste_point_a (scale_factor* point_a.X(), 
                                           scale_factor*point_a.Y().get(),
                                           scale_factor*point_a.Z().get());

            ChastePoint<3> chaste_point_b (scale_factor*point_b.X(),
                                           scale_factor*point_b.Y().get(),
                                           scale_factor*point_b.Z().get());
                        
            stimuli_applied.push_back( InitialStimulus(stimulus.Strength(), stimulus.Duration(), stimulus.Delay() ) );
            stimuled_areas.push_back( ChasteCuboid( chaste_point_a, chaste_point_b ) );
        }
    }
    catch (const xml_schema::exception& e)
    {
         std::cerr << e << std::endl;
         EXCEPTION("XML parsing error");
    }
}

template<unsigned PROBLEM_DIM>
void SetupProblem(AbstractCardiacProblem<3, PROBLEM_DIM>& rProblem,
             ConformingTetrahedralMesh<3,3>& rMesh)
{
    rProblem.SetMesh(&rMesh);
    rProblem.SetEndTime(simulation_duration);   // ms
    rProblem.SetPdeTimeStep(pde_time_step); // ms
    rProblem.SetPrintingTimeStep(printing_time_step); // ms
    rProblem.SetOutputDirectory(output_directory+"/results");
    rProblem.SetOutputFilenamePrefix("Chaste");
    rProblem.SetCallChaste2Meshalyzer(false);  
    
    rProblem.Initialise();
}

int main(int argc, char *argv[]) 
{
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
        
        // construct mesh. Note that mesh is measured in cm
        unsigned slab_nodes_width = (unsigned)round(slab_width/inter_node_space);
        unsigned slab_nodes_height = (unsigned)round(slab_height/inter_node_space);
       
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(slab_nodes_width,
                             slab_nodes_height,
                             slab_nodes_width,
                             true);
        // place at origin
        mesh.Translate(-(double)slab_nodes_width/2.0,
                       -(double)slab_nodes_height/2.0,
                       -(double)slab_nodes_width/2.0);
        // scale
        double mesh_scale_factor = inter_node_space*scale_factor;
        mesh.Scale(mesh_scale_factor, mesh_scale_factor, mesh_scale_factor);

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
            triangles_writer.WriteFilesUsingMesh(mesh);
            
            // copy input parameters file to results directory
            //system(("cp " + parameter_file + " " + output_dir_full_path).c_str());
        }
    
        ChasteSlabCellFactory cell_factory;
        cell_factory.SetMesh( &mesh );
        
        switch(domain)
        {
            case domain_type::Mono :
            {
                MonodomainProblem<3> mono_problem( &cell_factory );
                SetupProblem(mono_problem, mesh);

                mono_problem.GetMonodomainPde()->SetIntracellularConductivityTensor(intracellular_cond*identity_matrix<double>(3));
                mono_problem.Solve();
                break;
            }
            case domain_type::Bi :
            {
                BidomainProblem<3> bi_problem( &cell_factory );
                SetupProblem(bi_problem, mesh);

                bi_problem.GetBidomainPde()->SetIntracellularConductivityTensor(intracellular_cond*identity_matrix<double>(3));
                bi_problem.GetBidomainPde()->SetExtracellularConductivityTensor(extracellular_cond*identity_matrix<double>(3));               
                bi_problem.Solve();            
                break;
            }    
            default:
                EXCEPTION("Unknown domain type!!!");
        }

    }
    catch(Exception& e)
    {
        std::cerr << e.GetMessage() << "\n";
        return 1;
    }
    
    EventHandler::Headings();
    EventHandler::Report();

    return 0;    
}    


