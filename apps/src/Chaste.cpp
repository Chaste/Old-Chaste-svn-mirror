#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include "PetscTools.hpp"
#include <petscvec.h>
#include <vector>
#include "AbstractCardiacCellFactory.hpp"
#include "TimeStepper.hpp"
#include "MeshalyzerMeshWriter.cpp"
#include "TrianglesMeshWriter.cpp"
#include "SumStimulus.hpp"
#include <ctime>

#include "ChasteParameters.hpp"
#include <memory>


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


// Parameters fixed at compile time
const std::string  output_filename_prefix = "Run";
const double ode_time_step = 0.02;     // ms
const double pde_time_step = 0.02;     // ms
const double printing_time_step = 1; // ms

const double intracellular_cond = 1.75;
const double extracellular_cond = 7.0;

// Scale factor because Chaste code expects lengths in cm, but params use mm.
const double scale_factor = 1/10.0;

class ChasteSlabCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
    InitialStimulus *mpStimulus2;
    SumStimulus *mpSumStimulus;
    
public:
    ChasteSlabCellFactory() : AbstractCardiacCellFactory<3>(ode_time_step)
    {
        mpStimulus = new InitialStimulus(-25.5*1000, 0.5);
        mpStimulus2= new InitialStimulus(-25.5*1000, 0.5, quadrant_stimulus_delay);
        mpSumStimulus = new SumStimulus(mpStimulus, mpStimulus2);
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
        double x=mpMesh->GetNode(node)->GetPoint()[0];
        double z=mpMesh->GetNode(node)->GetPoint()[2];
        
        if ( x <= inter_node_space*scale_factor*(face_stimulus_width-slab_width) )
        {
            if (x<=0 && z<=0)
            {
                return CreateCellWithIntracellularStimulus(mpSumStimulus);
            }
            else
            {
                return CreateCellWithIntracellularStimulus(mpStimulus);
            }
        }
        else
        {
            if (x<=0 && z<=0)
            {
                return CreateCellWithIntracellularStimulus(mpStimulus2);
            }
            else
            {
                return CreateCellWithIntracellularStimulus(mpZeroStimulus);
            }
        }
    }
    
    ~ChasteSlabCellFactory(void)
    {
        delete mpStimulus;
        delete mpStimulus2;
        delete mpSumStimulus;
    }
};

void ReadParametersFromFile()
{
    try
    {
        std::auto_ptr<ChasteParameters::type> p_params(ChasteParameters(parameter_file));
        simulation_duration = p_params->SimulationDuration();
        slab_width = p_params->SlabWidth();     // mm
        slab_height = p_params->SlabHeight();   // mm
        inter_node_space = p_params->InterNodeSpace(); // mm
        face_stimulus_width = p_params->FaceStimulusWidth(); // mm
        quadrant_stimulus_delay = p_params->QuadrantStimulusDelay(); // ms
        output_directory = p_params->OutputDirectory();
        mesh_output_directory = p_params->MeshOutputDirectory();
        domain = p_params->Domain();
        ionic_model = p_params->IonicModel();
    }
    catch (const xml_schema::exception& e)
    {
         std::cerr << e << std::endl;
         EXCEPTION("XML parsing error");
    }
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
        
        if (argc!=2)
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
    
                mono_problem.SetMesh(&mesh);
                mono_problem.SetEndTime(simulation_duration);   // ms
                mono_problem.SetPdeTimeStep(pde_time_step); // ms
                mono_problem.SetPrintingTimeStep(printing_time_step); // ms
                mono_problem.SetOutputDirectory(output_directory+"/results");
                mono_problem.SetOutputFilenamePrefix("Chaste");
                mono_problem.SetCallChaste2Meshalyzer(false);  
                
                mono_problem.Initialise();
                mono_problem.GetMonodomainPde()->SetIntracellularConductivityTensor(intracellular_cond*identity_matrix<double>(3));
                mono_problem.Solve();
                break;
            }
            case domain_type::Bi :
            {
                BidomainProblem<3> bi_problem( &cell_factory );
    
                bi_problem.SetMesh(&mesh);
                bi_problem.SetEndTime(simulation_duration);   // ms
                bi_problem.SetPdeTimeStep(pde_time_step); // ms
                bi_problem.SetPrintingTimeStep(printing_time_step); // ms
                bi_problem.SetOutputDirectory(output_directory+"/results");
                bi_problem.SetOutputFilenamePrefix("Chaste");
                bi_problem.SetCallChaste2Meshalyzer(false);  
                
                bi_problem.Initialise();
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


