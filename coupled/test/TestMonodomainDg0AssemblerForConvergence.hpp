#ifndef _TESTMONODOMAINDG0ASSEMBLERFORCONVERGENCE_HPP_
#define _TESTMONODOMAINDG0ASSEMBLERFORCONVERGENCE_HPP_

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
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

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
        mpStimulus = new InitialStimulus(-1500, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
    {
        if (mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0)
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
public:

    void WriteTemporaryMeshFiles(std::string meshFilename, int middleNode)
    {
            int nb_of_eles = 2*middleNode;
            int nb_of_nodes = nb_of_eles+1;

            // Nodes file

            std::ofstream node_file((meshFilename+".node").c_str());
            node_file << std::scientific;

            node_file << nb_of_nodes << "\t1\t0\t1" << std::endl;
            
            for (int i = 0; i < nb_of_nodes; i++)
            {
                int b = 0;
                if ((i == 0) || (i == nb_of_nodes-1))
                {
                    b = 1;
                }

                node_file << i << "\t" << 0.08*i/nb_of_eles << "\t" << b << std::endl;
            }

            node_file.close();

            // Elements file

            std::ofstream ele_file((meshFilename+".ele").c_str());

            ele_file << std::scientific;

            ele_file << nb_of_eles << "\t2\t0" << std::endl;
            
            for (int i = 0; i < nb_of_eles; i++)
            {
                ele_file << i << "\t" << i << "\t" << i+1 << std::endl;
            }

            ele_file.close();
        }
    
    void TestMonodomainDg01DSpaceAndTime()
    {
//       const std::string mesh_filename = "testoutput/1D_0_to_1mm";
       const std::string mesh_filename = "/tmp/1D_0_to_0.8mm";
       //Invariant: middle_node*space_step = 0.04 cm
//       double space_step=0.01;   // cm
       //Loop over the space stepif (space_step < 0.0016)
       double prev_voltage_for_space = -999;   // To ensure that the first test fails
       bool converging_in_space = false;
       double space_step;   // cm
       double time_step;   // ms
       double probe_voltage;

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
        
                monodomain_problem.SetMeshFilename(mesh_filename);
                monodomain_problem.SetEndTime(200);   // 200 ms
//                monodomain_problem.SetEndTime(10);   // 10 ms  < --- Delete this!
                
                monodomain_problem.SetOutputDirectory("/tmp/MonoDg01D");
                monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_1d");
                monodomain_problem.SetPdeTimeStep(time_step);
                monodomain_problem.Initialise();

                std::cout<<"   Solving with a time step of "<<time_step<<" ms"<<std::endl;

                try
                {                
                    monodomain_problem.Solve();
    
                    double* p_voltage_array;
                    int lo, hi;
                    monodomain_problem.GetVoltageArray(&p_voltage_array, lo, hi); 
            
            
                    double my_voltage = 0.0;
                    if ((lo <= middle_node) && (middle_node < hi))
                    {
                        int local_index = middle_node-lo;
                        my_voltage = p_voltage_array[local_index];
                    }
                    MPI_Allreduce(&my_voltage, &probe_voltage, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        
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
            
                    monodomain_problem.RestoreVoltageArray(&p_voltage_array);
                }
                catch (Exception e)
                {
                    // An exception has been caught, meaning that the time step is too big, so half it
                
                    std::cout << "   >>> Convergence test: the time step is too big and, as a result, an exception has been thrown" << std::endl;
                
                    time_step *= 0.5;
                }               
              } while(!converging_in_time);    //do while: time_step   

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
        } while(!converging_in_space);   //do while: space_step
        
        std::cout<<"================================================================================"<<std::endl;

        std::cout << "Converged both in space ("<< space_step <<" cm) and time ("<< time_step << " ms)" << std::endl;

        TS_ASSERT_DELTA(space_step, 0.005, 0.0);
        TS_ASSERT_DELTA(time_step, 0.005, 0.0);
        TS_ASSERT_DELTA(probe_voltage, -10.339, 0.0001);
        // Note: the delta is because of floating point issues (!!)
    }    
};

#endif //_TESTMONODOMAINDG0ASSEMBLERFORCONVERGENCE_HPP_
