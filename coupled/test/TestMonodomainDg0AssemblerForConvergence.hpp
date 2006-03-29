#ifndef _TESTMONODOMAINDG0ASSEMBLERFORCONVERGENCE_HPP_
#define _TESTMONODOMAINDG0ASSEMBLERFORCONVERGENCE_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"

#include <cxxtest/TestSuite.h>
#include <petsc.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>

#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Node.hpp"
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "ColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"
#include "PropagationPropertiesCalculator.hpp"

#include "MonodomainPde.hpp"
#include "MockEulerIvpOdeSolver.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "AbstractMonodomainProblemStimulus.hpp"

class PointStimulus1D: public AbstractMonodomainProblemStimulus<1>
{
public:
    virtual void Apply(MonodomainPde<1> *pPde, 
                       ConformingTetrahedralMesh<1,1> *pMesh)
    {
//        static InitialStimulus stimulus(-600.0, 0.5);
// Note: the above stimulus used to be the original one and the one used in 
//       other test series. It, however, doesn't allow for small space steps,
//       hence we have increased its amplitude (but not too much, as otherwise
//       the cell model will blow up).
        static InitialStimulus stimulus(-1500.0, 0.5);

        for (int i = 0; i < pMesh->GetNumNodes(); i++)
        {
            if (pMesh->GetNodeAt(i)->GetPoint()[0] == 0.0)
            {
                pPde->SetStimulusFunctionAtNode(i, &stimulus);
            }
        }
    }
};

class TestMonodomainDg0Assembler : public CxxTest::TestSuite 
{   
private:
    /**
     * Refactor code to set up a PETSc vector holding the initial condition.
     */
    Vec CreateInitialConditionVec(int size)
    {
        Vec initial_condition;
        VecCreate(PETSC_COMM_WORLD, &initial_condition);
        VecSetSizes(initial_condition, PETSC_DECIDE, size);
        VecSetFromOptions(initial_condition);
        return initial_condition;
    }
    
public:

    // Solve on a 1D string of cells, 1mm long with a space step of 0.1mm and a
    // time step that starts at 0.01ms and gets halved until we get convergence.
    void TestMonodomainDg01D()
    {
        PointStimulus1D point_stimulus_1D;

        bool converging = false;
        double prev_voltage = -999;   // To ensure that the first test fails
        double time_step = 0.01;   // ms

        // Determine a converging time step

        do
        {
            MonodomainProblem<1> monodomainProblem("mesh/test/data/1D_0_to_1mm_10_elements",
                                                   2, // ms
                                                   "testoutput/MonoDg01D",
                                                   "NewMonodomainLR91_1d",
                                                   &point_stimulus_1D);
                                                   
            monodomainProblem.SetTimeSteps(time_step, time_step);
            monodomainProblem.Solve();
    
            double* voltage_array;
    
            // test whether voltages and gating variables are in correct ranges
    
            int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &voltage_array);
    
            for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
            {
                // assuming LR model has Ena = 54.4 and Ek = -77
                double Ena   =  54.4;   // mV
                double Ek    = -77.0;   // mV
    
                TS_ASSERT_LESS_THAN_EQUALS(   voltage_array[global_index-monodomainProblem.mLo] , Ena +  30);
                TS_ASSERT_LESS_THAN_EQUALS(  -voltage_array[global_index-monodomainProblem.mLo] + (Ek-30), 0);
    
                std::vector<double> odeVars = monodomainProblem.mMonodomainPde->GetOdeVarsAtNode(global_index);
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
            if ((monodomainProblem.mLo <= 5) && (monodomainProblem.mHi > 5))
            {
                my_voltage = voltage_array[5-monodomainProblem.mLo];
            }
            double probe_voltage;
            MPI_Allreduce(&my_voltage, &probe_voltage, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

            double relerr = fabs ((probe_voltage - prev_voltage) / prev_voltage);
// Note: the first time around this iteration, the relative error will be wrong
//       because of the initial value of prev_voltage
//            std::cout << "Testing for dt = " << monodomainProblem.time_step
//                << "ms (dx = 0.01cm) err = " << relerr << " V = " << voltage_array[5] << "mV"
//                << std::endl << std::flush;

            if (relerr < 1e-2)
            {
                converging = true;
//                std::cout << "Converged at dt = " << monodomainProblem.time_step << "ms" << std::endl;
            }
            else
            {
                // Get ready for the next test by halving the time step
    
                time_step *= 0.5;
            }

            prev_voltage = probe_voltage;
    
            VecRestoreArray(monodomainProblem.mCurrentVoltage, &voltage_array);      
            VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
            VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
            VecDestroy(monodomainProblem.mCurrentVoltage);
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
           
            MonodomainProblem<1> monodomainProblem(mesh_filename,
                                                   2, // ms
                                                   "testoutput/MonoDg01D",
                                                   "NewMonodomainLR91_1d",
                                                   &point_stimulus_1D);
                                                   
            monodomainProblem.SetTimeSteps( time_step, time_step);
            monodomainProblem.Solve();
    
            double* voltage_array;
    
            // test whether voltages and gating variables are in correct ranges
    
            int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &voltage_array);
    
            for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
            {
                // assuming LR model has Ena = 54.4 and Ek = -77
                double Ena   =  54.4;   // mV
                double Ek    = -77.0;   // mV
    
                TS_ASSERT_LESS_THAN_EQUALS(   voltage_array[global_index-monodomainProblem.mLo] , Ena +  30);
                TS_ASSERT_LESS_THAN_EQUALS(  -voltage_array[global_index-monodomainProblem.mLo] + (Ek-30), 0);
    
                std::vector<double> odeVars = monodomainProblem.mMonodomainPde->GetOdeVarsAtNode(global_index);
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
            if ((monodomainProblem.mLo <= middle_node) && (monodomainProblem.mHi > middle_node))
            {
                my_voltage = voltage_array[middle_node-monodomainProblem.mLo];
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
    
            VecRestoreArray(monodomainProblem.mCurrentVoltage, &voltage_array);      
            VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
            VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
            VecDestroy(monodomainProblem.mCurrentVoltage);
        } while(!converging);
        
        // Conclusion
        
//        std::cout << "In conclusion, the 'NewMonodomainLR91_1d' model converges for dx = " << space_step << "cm and dt = " << time_step << "ms" << std::endl;
        
        // Do we end up with the expected space step?
        
        TS_ASSERT_DELTA(space_step, 0.0025, 0.0);
    }
};

#endif //_TESTMONODOMAINDG0ASSEMBLERFORCONVERGENCE_HPP_
