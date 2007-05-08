#ifndef _TESTMONODOMAINDG0ASSEMBLERFORCONVERGENCE_HPP_
#define _TESTMONODOMAINDG0ASSEMBLERFORCONVERGENCE_HPP_


#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>

#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

// For chmod()
#include <sys/types.h>
#include <sys/stat.h>


class PointStimulusCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    InitialStimulus* mpStimulus;
    
public:
    PointStimulusCellFactory(double timeStep) : AbstractCardiacCellFactory<1>(timeStep)
    {
        // set the new stimulus
        //mpStimulus = new InitialStimulus(-600, 0.5);
// Note: the above stimulus used to be the original one and the one used in
//       other test series. It, however, doesn't allow for small space steps,
//       hence we have increased its amplitude (but not too much, as otherwise
//       the cell model will blow up).
        mpStimulus = new InitialStimulus(-1500*1000, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (mpMesh->GetNode(node)->GetPoint()[0] == 0.0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
        }
    }
    
    ~PointStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};


class TestMonodomainDg0AssemblerForConvergence : public CxxTest::TestSuite
{
private:
    OutputFileHandler *mpOutputFileHandler;
public:

    void WriteTemporaryMeshFiles(std::string meshFilename, int middleNode)
    {
        int nb_of_eles = 2*middleNode;
        int nb_of_nodes = nb_of_eles+1;
        
        // Nodes file
        out_stream p_node_file = mpOutputFileHandler->OpenOutputFile(meshFilename+".node");
        (*p_node_file) << std::scientific;
        
        (*p_node_file) << nb_of_nodes << "\t1\t0\t1" << std::endl;
        for (int i = 0; i < nb_of_nodes; i++)
        {
            int b = 0;
            if ((i == 0) || (i == nb_of_nodes-1))
            {
                b = 1;
            }
            (*p_node_file) << i << "\t" << 0.08*i/nb_of_eles << "\t" << b << std::endl;
        }
        p_node_file->close();
        
        // Elements file
        out_stream p_ele_file = mpOutputFileHandler->OpenOutputFile(meshFilename+".ele");
        (*p_ele_file) << std::scientific;
        
        (*p_ele_file) << nb_of_eles << "\t2\t0" << std::endl;
        for (int i = 0; i < nb_of_eles; i++)
        {
            (*p_ele_file) << i << "\t" << i << "\t" << i+1 << std::endl;
        }
        p_ele_file->close();
    }
    
    void TestMonodomainDg01DSpaceAndTime()
    {
        const std::string mesh_filename = "1D_0_to_0.8mm";
        OutputFileHandler output_file_handler("MonodomainConvergenceMesh");
        std::string mesh_pathname = output_file_handler.GetTestOutputDirectory()
                                    + mesh_filename;
        mpOutputFileHandler = &output_file_handler;
        
        //Invariant: middle_node*space_step = 0.04 cm
//       double space_step=0.01;   // cm
        //Loop over the space stepif (space_step < 0.0016)
        double prev_voltage_for_space = -999;   // To ensure that the first test fails
        bool converging_in_space = false;
        double space_step;   // cm
        double time_step;   // ms
        double probe_voltage;
        ReplicatableVector voltage_replicated;
        
        int middle_node = 1;
        
        do
//       for (int middle_node=1; middle_node<100; middle_node*=2){
        {
            bool converging_in_time = false;
            double prev_voltage_for_time = -999;   // To ensure that the first test fails
            
            space_step=0.04/middle_node;   // cm
            time_step = 0.04;   // ms
            
            WriteTemporaryMeshFiles(mesh_filename, middle_node);
            
            std::cout<<"================================================================================"<<std::endl;
            std::cout<<"Solving with a space step of "<<space_step<<" cm"<<std::endl;
            
            do
            {
                PointStimulusCellFactory cell_factory(time_step);
                MonodomainProblem<1> monodomain_problem(&cell_factory);
                
                monodomain_problem.SetMeshFilename(mesh_pathname);
                monodomain_problem.SetEndTime(200);   // 200 ms
                
                monodomain_problem.SetPdeTimeStep(time_step);
                monodomain_problem.Initialise();
                
                std::cout<<"   Solving with a time step of "<<time_step<<" ms"<<std::endl;
                
                try
                {
                    monodomain_problem.Solve();
                    
                    Vec voltage=monodomain_problem.GetVoltage();
                    voltage_replicated.ReplicatePetscVector(voltage);
                    
                    probe_voltage=voltage_replicated[middle_node];
                    
                    double relerr = fabs ((probe_voltage - prev_voltage_for_time) / prev_voltage_for_time);
                    std::cout<<"   >>> Convergence test: probe_voltage = "<<probe_voltage<<" mV | prev_voltage_for_time = "<<prev_voltage_for_time
                    <<" mV | relerr = "<<relerr<<std::endl;
                    
                    if (relerr < 1e-2)
                    {
                        converging_in_time = true;
                    }
                    else
                    {
                        // Get ready for the next test by halving the time step
                        
                        time_step *= 0.5;
                    }
                    
                    prev_voltage_for_time = probe_voltage;
                    
                }
                catch (Exception e)
                {
                    // An exception has been caught, meaning that the time step is too big, so halve it
                    std::cout << "   >>> Convergence test: an exception was thrown (" << e.GetMessage() << ")" << std::endl;
                    std::cout << "   >>>                   We assume that the time step was too big" << std::endl;
                    
                    time_step *= 0.5;
                }
            }
            while (!converging_in_time);   //do while: time_step
            
            double relerr = fabs ((probe_voltage - prev_voltage_for_space) / prev_voltage_for_space);
            std::cout<<">>> Convergence test: probe_voltage = "<<probe_voltage<<" mV | prev_voltage_for_space = "<<prev_voltage_for_space
            <<" mV | relerr = "<<relerr<<std::endl;
            
            if (relerr < 1e-2)
            {
                converging_in_space = true;
            }
            else
            {
                // Get ready for the next test by halving the space step (done by doubling the index of the middle node)
                middle_node*=2;
            }
            
            prev_voltage_for_space = probe_voltage;
        }
        while (!converging_in_space);   //do while: space_step
        
        std::cout<<"================================================================================"<<std::endl;
        
        std::cout << "Converged both in space ("<< space_step <<" cm) and time ("<< time_step << " ms)" << std::endl;
        
        TS_ASSERT_DELTA(space_step, 0.01, 0.0);
        TS_ASSERT_DELTA(time_step, 0.02, 0.0);
        TS_ASSERT_DELTA(probe_voltage, -10.3, 1.0);
        // Note: the delta is because of floating point issues (!!)
    }
    
    
    // here we keep pde_time_step constant (=0.01) and half the ode
    // timestep until we get convergence. Turns out we need the ode
    // timestep to be 10 times smaller than the pde timestep. Note
    // that ComputeExceptVoltage() method in the cell classes now
    // forces the voltage to remain constant using the mSetVoltageDerivativeToZero
    // flag - without it this test would fail and exception would even
    // be thrown.
    void TestMonodomainForConvergenceInOdeTimestep()
    {
        double pde_time_step = 0.01;
        double ode_time_step = 2*pde_time_step; // this is twice pde_time_step because ode_time_step
        // is halved at beginning to 'do' loop
        
        double rel_error = 999;
        double previous_voltage = -999;
        double probe_voltage = -999;
        
        double TOL = 1e-2;
        std::cout << "Tolerance for relative error is " << TOL << "\n\n";
        while ( rel_error > TOL )
        {
            ode_time_step *= 0.5;
            std::cout << "Ode timestep is " << ode_time_step << "\n\n";
            
            PointStimulusCellFactory cell_factory(ode_time_step);
            MonodomainProblem<1> monodomain_problem(&cell_factory);
            
            monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
            int end_node = 100; //101st node
            
            monodomain_problem.SetEndTime(16);   // ms - just long enough for AP to reach end node
            monodomain_problem.SetPdeTimeStep(pde_time_step);
            monodomain_problem.Initialise();
            
            TS_ASSERT_THROWS_NOTHING(monodomain_problem.Solve());
            
            Vec voltage = monodomain_problem.GetVoltage();
            ReplicatableVector voltage_replicated;
            voltage_replicated.ReplicatePetscVector(voltage);
            
            probe_voltage = voltage_replicated[end_node];
            rel_error = fabs ((probe_voltage - previous_voltage) / previous_voltage);
            std::cout << "voltage at end node is: " << probe_voltage    << "\n"
            << "previous voltage was  : " << previous_voltage << "\n"
            << "relative error is     : " << rel_error << "\n\n";
            
            previous_voltage = probe_voltage;
            
            
            // force quit if ode_time_step gets very small, to avoid infinite loop
            if (ode_time_step<1e-4)
            {
                TS_FAIL("");
            }
        }
        
        /* !!-- note that 0.001 (ie pde_timestep/10) would also have worked --!! */
        TS_ASSERT_DELTA(ode_time_step, 0.000625, 0.0)
        
        // finally, verify the AP has reaching the end node
        TS_ASSERT_LESS_THAN(0.0, probe_voltage);
    }
};

#endif //_TESTMONODOMAINDG0ASSEMBLERFORCONVERGENCE_HPP_
