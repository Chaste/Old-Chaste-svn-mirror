#ifndef TESTBIDOMAINFORCONVERGENCE1D_HPP_
#define TESTBIDOMAINFORCONVERGENCE1D_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>

#include "ConformingTetrahedralMesh.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "OutputFileHandler.hpp"


class PointStimulusCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    InitialStimulus* mpStimulus;
    
public:
    PointStimulusCellFactory(double timeStep) : AbstractCardiacCellFactory<1>(timeStep)
    {
        mpStimulus = new InitialStimulus(-1500000, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
    {
        if (mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0)
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


class TestBidomainForConvergence1d : public CxxTest::TestSuite 
{
private:
    OutputFileHandler *mpOutputFileHandler;
public:

    void WriteTemporaryMeshFiles(std::string meshFilename, int middleNode)
    {
        int nb_of_eles = 2*middleNode;
        int nb_of_nodes = nb_of_eles+1;

        // Nodes file
        std::ofstream *p_node_file = mpOutputFileHandler->OpenOutputFile(meshFilename+".node");
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
        delete p_node_file;

        // Elements file
        std::ofstream *p_ele_file = mpOutputFileHandler->OpenOutputFile(meshFilename+".ele");
        (*p_ele_file) << std::scientific;

        (*p_ele_file) << nb_of_eles << "\t2\t0" << std::endl;
        for (int i = 0; i < nb_of_eles; i++)
        {
            (*p_ele_file) << i << "\t" << i << "\t" << i+1 << std::endl;
        }
        p_ele_file->close();
        delete p_ele_file;
    }
    
    
    void TestBidomain1DSpaceAndTime()
    {
       const std::string mesh_filename = "1D_0_to_0.8mm";
       OutputFileHandler output_file_handler("BidomainConvergenceMesh");
       mpOutputFileHandler = &output_file_handler;

       //Invariant: middle_node*space_step = 0.04 cm
//       double space_step=0.01;   // cm
       //Loop over the space stepif (space_step < 0.0016)
       double prev_voltage_for_space = -999;   // To ensure that the first test fails
       bool converging_in_space = false;
       double space_step;   // cm
       double time_step;   // ms
       double probe_voltage;

       int middle_node = 4;

       do
//       for (int middle_node=1; middle_node<100; middle_node*=2){
       {
            bool converging_in_time = false;
            double prev_voltage_for_time = -999;   // To ensure that the first test fails

            space_step=0.04/middle_node;   // cm
            time_step = 0.04;   // ms

            WriteTemporaryMeshFiles(mesh_filename, middle_node);
            

            std::cout<<"================================================================================"<<std::endl  << std::flush;
            std::cout<<"Solving with a space step of "<<space_step<<" cm"<<std::endl  << std::flush;

            do
            {
                PointStimulusCellFactory cell_factory(time_step);
                BidomainProblem<1> bidomain_problem(&cell_factory);
        
                bidomain_problem.SetMeshFilename(mesh_filename);
                bidomain_problem.SetEndTime(200);   // 200 ms
                
                bidomain_problem.SetPdeTimeStep(time_step);
                bidomain_problem.Initialise();
                
                std::cout<<"   Solving with a time step of "<<time_step<<" ms"<<std::endl  << std::flush;

                try
                {                
                    bidomain_problem.Solve();
                    
 
    
                    double* p_voltage_array;
                    int lo, hi;
                    bidomain_problem.GetVoltageArray(&p_voltage_array, lo, hi); 
            
            
                    double my_voltage = 0.0;
                    if ((lo <= middle_node) && (middle_node < hi))
                    {
                        int local_index = middle_node-lo;
                        my_voltage = p_voltage_array[local_index];
                    }
                    MPI_Allreduce(&my_voltage, &probe_voltage, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        
                    double relerr = fabs ((probe_voltage - prev_voltage_for_time) / prev_voltage_for_time);
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
        
                    prev_voltage_for_time = probe_voltage;
            
                    bidomain_problem.RestoreVoltageArray(&p_voltage_array);
                }
                catch (Exception e)
                {
                    // An exception has been caught, meaning that the time step is too big, so halve it
                    std::cout << "   >>> Convergence test: an exception was thrown (" << e.GetMessage() << ")" << std::endl  << std::flush;
                    std::cout << "   >>>                   We assume that the time step was too big" << std::endl << std::flush;
                
                    time_step *= 0.5;
                }               
            } while(!converging_in_time);    //do while: time_step   

            double relerr = fabs ((probe_voltage - prev_voltage_for_space) / prev_voltage_for_space);
            std::cout<<">>> Convergence test: probe_voltage = "<<probe_voltage<<" mV | prev_voltage_for_space = "<<prev_voltage_for_space
                     <<" mV | relerr = "<<relerr<<std::endl << std::flush;

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
        } while (!converging_in_space);   //do while: space_step
        
        std::cout<<"================================================================================"<<std::endl << std::flush;

        std::cout << "Converged both in space ("<< space_step <<" cm) and time ("<< time_step << " ms)" << std::endl << std::flush;

        TS_ASSERT_DELTA(space_step, 0.005, 0.0);
        TS_ASSERT_DELTA(time_step, 0.005, 0.0);
        TS_ASSERT_DELTA(probe_voltage, -10.343, 0.0001);
        // Note: the delta is because of floating point issues (!!)
    }    
};


#endif /*TESTBIDOMAINFORCONVERGENCE1D_HPP_*/
