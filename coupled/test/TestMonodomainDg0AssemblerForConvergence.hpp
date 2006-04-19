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
#include "MonodomainProblemIteration7.hpp"
#include "AbstractCardiacCellFactory.hpp"


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
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus, mTimeStep);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus, mTimeStep);
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

    // Solve on a 1D string of cells, 1mm long with a space step of 0.1mm and a
    // time step that starts at 0.01ms and gets halved until we get convergence.
    void TestMonodomainDg01D()
    {
        bool converging = false;
        double prev_voltage = -999;   // To ensure that the first test fails
        double time_step = 0.01;   // ms

        // Determine a converging time step

        do
        {
            PointStimulusCellFactory cell_factory(time_step);
            MonodomainProblemIteration7<1> monodomain_problem(&cell_factory);
    
            monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
            monodomain_problem.SetEndTime(2);   // 2 ms
            monodomain_problem.SetOutputDirectory("testoutput/MonoDg01D");
            monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_1d");
            monodomain_problem.SetPdeTimeStep(time_step);
            monodomain_problem.Initialise();
            monodomain_problem.Solve();
    
            double* voltage_array;
            int lo, hi;
            monodomain_problem.GetVoltageArray(&voltage_array, lo, hi); 
    
            // test whether voltages and gating variables are in correct ranges
            for(int global_index=lo; global_index<hi; global_index++) 
            {
                // assuming LR model has Ena = 54.4 and Ek = -77
                double Ena   =  54.4;   // mV
                double Ek    = -77.0;   // mV
    
                TS_ASSERT_LESS_THAN_EQUALS(   voltage_array[global_index-monodomain_problem.mLo] , Ena +  30);
                TS_ASSERT_LESS_THAN_EQUALS(  -voltage_array[global_index-monodomain_problem.mLo] + (Ek-30), 0);
    
                std::vector<double> odeVars = monodomain_problem.GetCardiacCell(global_index)->rGetStateVariables();
                for(int j=0; j<8; j++)
                {
                    // if not voltage or calcium ion conc, test whether between 0 and 1
                    if((j!=4) && (j!=3))
                    {
                        TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);
                        TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);
                    }
                }
            }

            double my_voltage = 0.0;
            if ((monodomain_problem.mLo <= 5) && (monodomain_problem.mHi > 5))
            {
                my_voltage = voltage_array[5-monodomain_problem.mLo];
            }
            double probe_voltage;
            MPI_Allreduce(&my_voltage, &probe_voltage, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

            double relerr = fabs ((probe_voltage - prev_voltage) / prev_voltage);
// Note: the first time around this iteration, the relative error will be wrong
//       because of the initial value of prev_voltage
//            std::cout << "Testing for dt = " << monodomain_problem.time_step
//                << "ms (dx = 0.01cm) err = " << relerr << " V = " << voltage_array[5] << "mV"
//                << std::endl << std::flush;

            if (relerr < 1e-2)
            {
                converging = true;
//                std::cout << "Converged at dt = " << monodomain_problem.time_step << "ms" << std::endl;
            }
            else
            {
                // Get ready for the next test by halving the time step
    
                time_step *= 0.5;
            }

            prev_voltage = probe_voltage;
    
            monodomain_problem.RestoreVoltageArray(&voltage_array);
          } while(!converging);
                
        // Do we end up with the expected time step?
        
        TS_ASSERT_DELTA(time_step, 0.00125, 0.0);

        // Determine a converging space step, using the converging time step 
        // above

        int middle_node = 5;
        double space_step = 0.01;   // cm
        converging = false;
        prev_voltage = -999;   // To ensure that the first test fails
        const std::string mesh_filename = "testoutput/1D_0_to_1mm";
        
        PointStimulusCellFactory cell_factory(time_step);

        do
        {
            int nb_of_eles = 2*middle_node;
            int nb_of_nodes = nb_of_eles+1;

            // Nodes file

            std::ofstream node_file((mesh_filename+".node").c_str());

            node_file << std::scientific;

            node_file << nb_of_nodes << "\t1\t0\t1" << std::endl;
            
            for (int i = 0; i < nb_of_nodes; i++)
            {
                int b = 0;
                if ((i == 0) || (i == nb_of_nodes-1))
                {
                    b = 1;
                }

                node_file << i << "\t" << 0.1*i/nb_of_eles << "\t" << b << std::endl;
            }

            node_file.close();

            // Elements file

            std::ofstream ele_file((mesh_filename+".ele").c_str());

            ele_file << std::scientific;

            ele_file << nb_of_eles << "\t2\t0" << std::endl;
            
            for (int i = 0; i < nb_of_eles; i++)
            {
                ele_file << i << "\t" << i << "\t" << i+1 << std::endl;
            }

            ele_file.close();

            // Solve the problem itself
           
            MonodomainProblemIteration7<1> monodomain_problem(&cell_factory);
    
            monodomain_problem.SetMeshFilename(mesh_filename);
            monodomain_problem.SetEndTime(2);   // 2 ms
            monodomain_problem.SetOutputDirectory("testoutput/MonoDg01D");
            monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_1d");
            monodomain_problem.SetPdeTimeStep(time_step);
            monodomain_problem.Initialise();
            monodomain_problem.Solve();
    
            double* voltage_array;
            int lo, hi;
            monodomain_problem.GetVoltageArray(&voltage_array, lo, hi); 
    
            // test whether voltages and gating variables are in correct ranges
            for(int global_index=lo; global_index<hi; global_index++)
            {
                // assuming LR model has Ena = 54.4 and Ek = -77
                double Ena   =  54.4;   // mV
                double Ek    = -77.0;   // mV
    
                TS_ASSERT_LESS_THAN_EQUALS(   voltage_array[global_index-monodomain_problem.mLo] , Ena +  30);
                TS_ASSERT_LESS_THAN_EQUALS(  -voltage_array[global_index-monodomain_problem.mLo] + (Ek-30), 0);
    
                std::vector<double> odeVars = monodomain_problem.GetCardiacCell(global_index)->rGetStateVariables();
                for(int j=0; j<8; j++)
                {
                    // if not voltage or calcium ion conc, test whether between 0 and 1
                    if((j!=4) && (j!=3))
                    {
                        TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);
                        TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);
                    }
                }
            }

            double my_voltage = 0.0;
            if ((monodomain_problem.mLo <= middle_node) && (monodomain_problem.mHi > middle_node))
            {
                my_voltage = voltage_array[middle_node-monodomain_problem.mLo];
            }
            double probe_voltage;
            MPI_Allreduce(&my_voltage, &probe_voltage, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

            double relerr = fabs ((probe_voltage - prev_voltage) / prev_voltage);
// Note: the first time around this iteration, the relative error will be wrong
//       because of the initial value of prev_voltage
//            std::cout << "Testing for dx = " << space_step
//                << "cm (dt = " << time_step << "ms) err = " << relerr << " V = " 
//                << voltage_array[middle_node] << "mV" << std::endl << std::flush;
            
            prev_voltage = probe_voltage;

            if (relerr < 1e-2)
            {
                converging = true;
//                std::cout << "Converged at dx = " << space_step << "cm" << std::endl;
            }
            else
            {
                // Get ready for the next test by halving the space step and 
                // doubling the middle node
    
                space_step *= 0.5;
                middle_node *= 2;
            }
    
             monodomain_problem.RestoreVoltageArray(&voltage_array);
        } while(!converging);
        
        // Conclusion
        
//        std::cout << "In conclusion, the 'NewMonodomainLR91_1d' model converges for dx = " << space_step << "cm and dt = " << time_step << "ms" << std::endl;
        
        // Do we end up with the expected space step?
        
        TS_ASSERT_DELTA(space_step, 0.0025, 0.0);
    }
};


#endif //_TESTMONODOMAINDG0ASSEMBLERFORCONVERGENCE_HPP_
