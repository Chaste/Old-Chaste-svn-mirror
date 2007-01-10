#ifndef TEST3DBIDOMAINFORCONVERGENCE_HPP_
#define TEST3DBIDOMAINFORCONVERGENCE_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "OutputFileHandler.hpp"


class PointStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    // define a new stimulus
    InitialStimulus* mpStimulus;
    
public:
    PointStimulusCellFactory(double timeStep, double num_elements) : AbstractCardiacCellFactory<3>(timeStep)
    {
    	// scale stimulus depending on space_step of elements
        mpStimulus = new InitialStimulus(-100000* pow(num_elements/12.0, 1.0/3.0), 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        double x = mpMesh->GetNode(node)->GetPoint()[0];
        //double y = mpMesh->GetNode(node)->GetPoint()[1];
        //double z = mpMesh->GetNode(node)->GetPoint()[2];
        
        if (x*x<=1e-10)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus, mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus, mpZeroStimulus);
        }
    }
    
    ~PointStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};


class Test3dBidomainForConvergence : public CxxTest::TestSuite
{
public:

    void Test3dBidomainSpaceAndTime()
    {
    	double num_elements[5] =           {  12,   152,   1016,   7790 ,  61687 };
        std::string file_num_elements[5] = { "12", "152", "1016", "7790", "61687"}; //, "488702"};
        double approx_space_steps[5] = {0.043,0.018,0.01,0.005,0.0025}; //, 0.00125};
        int opposite_corner_node = 6; // the node at (0.2,0.2,0.2)
        int num_meshes=5;                
        // To ensure that the first test fails        
        double prev_voltage_for_space = -999;   
        bool converging_in_space = false;
        bool failed_to_converge_in_space = false;

        double time_step;   // ms

        double probe_voltage;
        ReplicatableVector voltage_replicated;
        
        int current_file_num = 0;                    
        
        do //do while: space_step
        {
            bool converging_in_time = false;
            // To ensure that the first test fails
            double prev_voltage_for_time = -999;   
            
            time_step = 0.04;  // ms 
            
            
            std::string mesh_pathname = "mesh/test/data/cube_2mm_" + file_num_elements[current_file_num] + "_elements";
         

            std::cout<<"================================================================================"<<std::endl  << std::flush;
            std::cout<<"Solving with an (approx) space step of "<< approx_space_steps[current_file_num] <<" cm - mesh " << current_file_num <<std::endl  << std::flush;
            
            do  //do while: time_step
            {
                PointStimulusCellFactory cell_factory(time_step, num_elements[current_file_num]);
                BidomainProblem<3> bidomain_problem(&cell_factory);
                
                bidomain_problem.SetMeshFilename(mesh_pathname);
                bidomain_problem.SetEndTime(3.52);   // ms        

                
/*              bidomain_problem.SetOutputDirectory("bitemp");
                bidomain_problem.SetOutputFilenamePrefix("bitemp");
                bidomain_problem.SetPrintingTimeStep(0.1);
                bidomain_problem.SetWriteInfo();
*/

                bidomain_problem.SetPdeTimeStep(time_step);
                bidomain_problem.Initialise();

                std::cout<<"   Solving with a time step of "<<time_step<<" ms"<<std::endl  << std::flush;
                
                try
                {
                    bidomain_problem.Solve();
                    Vec voltage=bidomain_problem.GetVoltage();
                    voltage_replicated.ReplicatePetscVector(voltage);
                    
                    probe_voltage=voltage_replicated[opposite_corner_node];
                    
                    double relerr = fabs ((pro
    {
    	double num_elements[5] =           {  12,   152,   1016,   7790 ,  61687 };
        std::string file_num_elements[5] = { "12", "152", "1016", "7790", "61687"}; //, "488702"};
        double approx_space_steps[5] = {0.043,0.018,0.01,0.005,0.0025}; //, 0.00125};
        int opposite_corner_node = 6; // the node at (0.2,0.2,0.2)
        int num_meshes=3;                
        // Tobe_voltage - prev_voltage_for_time) / prev_voltage_for_time);
                    std::cout<<"   >>> Convergence test: probe_voltage = "<<probe_voltage<<" mV | prev_voltage_for_time = "<<prev_voltage_for_time
                    <<" mV | relerr = "<<relerr<<std::endl  << std::flush;
                    
                    if (relerr < 1e-2)
                    {
                        converging_in_time = true;
                    }
                    else 
                    {
                        // Get ready for the next test by halving the time step
                        time_step *= 0.5;
                    }
                    
                    if(time_step < 1e-4)
                    {
                        std::cout << "**** NO TIMESTEP GREATER THAN 1e-4 FOUND WHICH WORKS, MOVING ONTO NEXT MESH...****\n" << std::flush;
                        converging_in_time = true;
                    }
                    
                    prev_voltage_for_time = probe_voltage;
                    
                }
                catch (Exception e)
                {
                    // An exception has been caught, meaning that the time step is too big, so halve it
                    std::cout << "   >>> Convergence test: an exception was thrown (" << e.GetMessage() << ")" << std::endl  << std::flush;
                    std::cout << "   >>>                   We assume that the time step was too big" << std::endl << std::flush;
                    
                    time_step *= 0.5;
                }
            }
            while (!converging_in_time);   //do while: time_step
            
            
            double relerr = fabs ((probe_voltage - prev_voltage_for_space) / prev_voltage_for_space);
            std::cout<<">>> Convergence test: probe_voltage = "<<probe_voltage<<" mV | prev_voltage_for_space = "<<prev_voltage_for_space
            <<" mV | relerr = "<<relerr<<std::endl << std::flush;
            
            if (relerr < 1e-2)
            {
                converging_in_space = true;
            }
            else
            {
                // Use the next mesh next time 
                current_file_num++;
                if(current_file_num==num_meshes)
                {
                    TS_FAIL("Could not converge for any of the meshes used");
                    failed_to_converge_in_space = true;
                }
            }
            
            prev_voltage_for_space = probe_voltage;
        }
        while (!converging_in_space && !failed_to_converge_in_space);   //do while: space_step
        
        if (converging_in_space)
        {
	        std::cout<<"================================================================================"<<std::endl << std::flush;
        
    	    std::cout << "Converged both in space ("<< approx_space_steps[current_file_num] <<" cm) and time ("<< time_step << " ms)" << std::endl << std::flush;
        }    
        
       // TS_ASSERT_DELTA(approx_space_steps[current_file_num], 0.005, 0.0);
       // TS_ASSERT_DELTA(time_step, 0.005, 0.0);
       // TS_ASSERT_DELTA(probe_voltage, -10.3432, 0.0001);
        // Note: the delta is because of floating point issues (!!)
    }
};




#endif /*TEST3DBIDOMAINFORCONVERGENCE_HPP_*/
