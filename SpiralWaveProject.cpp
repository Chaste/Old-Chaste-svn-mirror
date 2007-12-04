#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include "AbstractCardiacCellFactory.hpp"
#include "TimeStepper.hpp"
#include "MeshalyzerMeshWriter.cpp"
#include "SumStimulus.hpp"
#include <ctime>

#include "SpiralParameters.hpp"
#include <memory>

//#include "FoxModel2002.hpp"
#include "BackwardEulerFoxModel2002Modified.hpp"

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

// Parameters fixed at compile time
const std::string  cell_var_name = "Ca_i"; // Variable to output to results files
const std::string  output_filename_prefix = "Run";
const double ode_time_step = 0.02;     // ms
const double pde_time_step = 0.02;     // ms
const double printing_time_step = 0.1; // ms

// Scale factor because Chaste code expects lengths in cm, but params use mm.
const double scale_factor = 1/10.0;

class SpiralWaveCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
    InitialStimulus *mpStimulus2;
    SumStimulus *mpSumStimulus;
public:
    SpiralWaveCellFactory() : AbstractCardiacCellFactory<3>(ode_time_step)
    {
        mpStimulus = new InitialStimulus(-600.0*1000, 0.5);
        mpStimulus2= new InitialStimulus(-600.0*1000, 0.5, quadrant_stimulus_delay);
        mpSumStimulus = new SumStimulus(mpStimulus, mpStimulus2);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        double x=mpMesh->GetNode(node)->GetPoint()[0];
        double z=mpMesh->GetNode(node)->GetPoint()[2];
        if ( x <= inter_node_space*scale_factor*(face_stimulus_width-slab_width) )
        {
            if (x<=0 && z<=0)
            {
                return new BackwardEulerFoxModel2002Modified(mTimeStep, mpSumStimulus);
            }
            else
            {
                return new BackwardEulerFoxModel2002Modified(mTimeStep, mpStimulus);
            }
        }
        else
        {
            if (x<=0 && z<=0)
            {
                return new BackwardEulerFoxModel2002Modified(mTimeStep, mpStimulus2);
            }
            else
            {
                return new BackwardEulerFoxModel2002Modified(mTimeStep, mpZeroStimulus);
            }
        }
    }
    
    ~SpiralWaveCellFactory(void)
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
        std::auto_ptr<SpiralParameters::type> p_params(SpiralParameters(parameter_file));
        simulation_duration = p_params->SimulationDuration();
        slab_width = p_params->SlabWidth();     // mm
        slab_height = p_params->SlabHeight();   // mm
        inter_node_space = p_params->InterNodeSpace(); // mm
        face_stimulus_width = p_params->FaceStimulusWidth(); // mm
        quadrant_stimulus_delay = p_params->QuadrantStimulusDelay(); // ms
        output_directory = p_params->OutputDirectory();
        mesh_output_directory = p_params->MeshOutputDirectory();
    }
    catch (const xml_schema::exception& e)
    {
         std::cerr << e << std::endl;
         EXCEPTION("XML parsing error");
    }
}

int main(int argc, char *argv[])
{
    PETSCEXCEPT(PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL) );
    
    // solver and preconditioner options
    //PetscOptionsSetValue("-ksp_type", "cg");
    //PetscOptionsSetValue("-pc_type", "bjacobi");
    //PetscOptionsSetValue("-options_table", "");
    
    if (argc!=2)
    {
        std::cout  << "Usage: SpiralWaveProject parameters_file\n";
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
    
    SpiralWaveCellFactory cell_factory;
    
    cell_factory.SetMesh( &mesh );
    
    MonodomainPde<3> monodomain_pde(&cell_factory);
    
    // Set up assembler
    MonodomainDg0Assembler<3, 3> monodomain_assembler(&mesh, &monodomain_pde);
  
    // initial condition;
    Vec initial_condition= DistributedVector::CreateVec();
    DistributedVector ic(initial_condition);

    // Set a constant initial voltage throughout the mMesh and get the state variable number and name
    for (DistributedVector::Iterator index = DistributedVector::Begin();
         index != DistributedVector::End();
         ++index)
         {
            ic[index] = monodomain_pde.GetCardiacCell(index.Global)->GetVoltage();
         }
    ic.Restore();
    

    AbstractCardiacCell& some_cell =  *(monodomain_pde.GetCardiacCell(DistributedVector::Begin().Global));
    unsigned cell_var_number = some_cell.GetStateVariableNumberByName(cell_var_name);
    std::string cell_var_units = some_cell.GetStateVariableUnitsByNumber(cell_var_number);

    
    ParallelColumnDataWriter writer(output_directory,output_filename_prefix);
        
    writer.DefineFixedDimension("Node", "dimensionless", mesh.GetNumNodes() );
    unsigned time_var_id = writer.DefineUnlimitedDimension("Time","msecs");
    unsigned voltage_var_id = writer.DefineVariable("V","mV");
    unsigned cell_var_id = writer.DefineVariable(cell_var_name,cell_var_units);
    writer.EndDefineMode();
    writer.PutVariable(time_var_id, 0.0);
    writer.PutVector(voltage_var_id, initial_condition); 
    for (DistributedVector::Iterator index = DistributedVector::Begin();
         index != DistributedVector::End();
         ++index)
    {
        writer.PutVariable(cell_var_id,
                           monodomain_pde.GetCardiacCell(index.Global)->GetStateVariableValueByNumber(cell_var_number),
                           index.Global);
    }   
    Vec voltage;
    Vec extra_var = DistributedVector::CreateVec();
    DistributedVector dist_extra_var(extra_var);
    // main solve loop
    TimeStepper stepper(0.0, simulation_duration, printing_time_step);
    while (!stepper.IsTimeAtEnd())
    {
        monodomain_assembler.SetTimes(stepper.GetTime(), stepper.GetNextTime(), pde_time_step);
        monodomain_assembler.SetInitialCondition( initial_condition );
        
       
        try
        {
             voltage=monodomain_assembler.Solve();
        }
        catch (Exception &e)
        {
            writer.Close();
            throw e;
        }
        writer.AdvanceAlongUnlimitedDimension(); // creates a new file
        writer.PutVariable(time_var_id, stepper.GetTime());
        writer.PutVector(voltage_var_id, voltage);
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            writer.PutVariable(cell_var_id,
                               monodomain_pde.GetCardiacCell(index.Global)->GetStateVariableValueByNumber(cell_var_number),
                               index.Global);
        }

        VecDestroy(initial_condition);
        initial_condition=voltage;
        stepper.AdvanceOneTimeStep();
    }
    writer.Close();
    
    // write out the mesh that was used if we are the master process
    if (DistributedVector::IsGlobalIndexLocal(0))
    {
        MeshalyzerMeshWriter<3,3> mesh_writer(mesh_output_directory, "SlabMesh", false);
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
}
