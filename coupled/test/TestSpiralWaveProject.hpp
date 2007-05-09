#ifndef TESTSPIRALWAVEPROJECT_HPP_
#define TESTSPIRALWAVEPROJECT_HPP_

#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "TimeStepper.hpp"
#include "MeshalyzerMeshWriter.cpp"
#include "SumStimulus.hpp"
#include <time.h>

const double simulation_duration = 5.0;    // ms *

const double slab_width = 2;     // mm *
const double slab_height= 1;     // mm *
const double inter_node_space = 0.25;// mm *
const double face_stimulus_width = 0.25; // mm *
const double quadrant_stimulus_delay = 120.0; // ms

const std::string  output_directory="SpiralWave"; // *
const std::string  mesh_output_directory="Slab"; // *
const std::string  cell_var_name="CaI"; // *

const std::string  output_filename_prefix="Run";
const double ode_time_step=0.01;    // ms
const double pde_time_step=0.01;    // ms
const double printing_time_step =0.1;// ms
double scale_factor = inter_node_space/10.0;

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
        if ( x <= scale_factor*(face_stimulus_width-slab_width) )
        {
            if (x<=0 && z<=0)
            {
                return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpSumStimulus);
            }
            else
            {
                return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus);
            }
        }
        else
        {
            if (x<=0 && z<=0)
            {
                return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus2);
            }
            else
            {
                return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
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

class TestSpiralWaveProject : public CxxTest::TestSuite
{
public:

    void TestMonodomain3DSlab() throw (Exception)
    {
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
        double scale_factor = inter_node_space/10.0;
        mesh.Scale(scale_factor, scale_factor, scale_factor);
        
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
        
        // and write out the mesh that was used
        MeshalyzerMeshWriter<3,3> mesh_writer(mesh_output_directory, "SlabMesh" );
        mesh_writer.WriteFilesUsingMesh(mesh);
        
    }
        
};

#endif /*TESTSPIRALWAVEPROJECT_HPP_*/
